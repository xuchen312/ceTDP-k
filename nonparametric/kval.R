#!/usr/bin/env Rscript

# =============================================================
# Nonparametric Cluster Extent Threshold Computation
#
# Example:
# Rscript kval.R \
#   --dirpath=/data/imgs \
#   --mask=/data/mask.nii.gz \
#   --conn=18 \
#   --nperm=1000 \
#   --alpha=0.05 \
#   --kval=14 \
#   --cft=3.1 \
#   --out=results.txt
# =============================================================

# -----------------------------
# load packages
# -----------------------------
if (!requireNamespace("optparse", quietly=TRUE)) 
  install.packages("optparse", repos="https://cloud.r-project.org")
if (!requireNamespace("RNifti", quietly=TRUE)) 
  install.packages("RNifti", repos="https://cloud.r-project.org")
if (!requireNamespace("Rcpp", quietly=TRUE)) 
  install.packages("Rcpp", repos="https://cloud.r-project.org")
library(optparse)
library(RNifti)
library(Rcpp)

# -----------------------------
# detect script directory
# -----------------------------
args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_full[grep("--file=", args_full)])
script_dir <- dirname(normalizePath(script_path))

# -----------------------------
# options
# -----------------------------
option_list <- list(
  make_option(c("-d","--dirpath"), type="character", help="Directory containing cope images"),
  make_option(c("-m","--mask"), type="character", help="Mask NIfTI file"),
  make_option(c("-c","--conn"), type="integer", default=18, help="3D connectivity criterion (6/18/26)"),
  make_option(c("-n","--nperm"), type="integer", default=1000, help="Number of permutations"),
  make_option(c("-a","--alpha"), type="double", default=0.05, help="Significance level"),
  make_option(c("-k","--kval"), type="integer", help="Cluster extent threshold (optional)"),
  make_option(c("-t","--cft"), type="double", default=3.1, help="Cluster forming threshold (optional)"),
  make_option(c("-o","--out"), type="character", default="results.txt", help="Output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

# -----------------------------
# checking
# -----------------------------
if (is.null(opt$dirpath)) stop("--dirpath is required")
if (is.null(opt$mask)) stop("--mask is required")

if (!opt$conn %in% c(6,18,26)) stop("--conn must be 6 / 18 / 26")
if (is.null(opt$k) && is.null(opt$cft)) stop("You must provide either --kval or --cft (or both)")

# -----------------------------
# load cpp (same folder as script)
# -----------------------------
cppfile <- file.path(script_dir, "kval.cpp")
if (!file.exists(cppfile)) stop("Cannot find kval.cpp in script directory")
cat("Compiling C++:", cppfile, "\n")
Rcpp::sourceCpp(cppfile)

# -----------------------------
# file loading
# -----------------------------
files <- list.files(path=opt$dirpath, pattern="\\.nii(\\.gz)?$", recursive=TRUE, full.names=TRUE)
if (length(files) == 0) stop("No cope images found")
mask <- RNifti::readNifti(opt$mask)

n <- length(files)
dims <- RNifti::niftiHeader(files[1])$dim[2:4]

cat("Subjects:", n, "\n")
cat("Dims:", paste(dims, collapse="x"), "\n")

# -----------------------------
# load data
# -----------------------------
cat("Loading data...\n")

cope <- array(0, c(dims, n))
for (i in seq_along(files)) {
  cope[,,,i] <- RNifti::readNifti(files[i])
}

# -----------------------------
# compute stats
# -----------------------------
V <- array(cope, dim=c(prod(dims), n))  # vectorized cope

mu  <- rowMeans(V)
sdv <- sqrt(rowMeans(V^2) - mu^2)

valid <- (sdv > 0) & (as.vector(mask) != 0)

tvals <- numeric(length(mu))
tvals[valid] <- sqrt(n) * mu[valid] / sdv[valid]

pvals <- rep(1, length(mu))
pvals[valid] <- pt(tvals[valid], df=n-1, lower.tail=FALSE)

zvals <- rep(0, length(mu))
zvals[valid] <- qnorm(pvals[valid], lower.tail=FALSE)

mask <- array(0, dims)
pmap <- array(1, dims)
zmap <- array(0, dims)

mask[valid] <- 1
pmap[valid] <- pvals[valid]
zmap[valid] <- zvals[valid]

# ------------------------------------
# cluster prep
# ------------------------------------

# mask indices (only compute once)
mask_idx <- which(mask != 0)
m <- length(mask_idx)

# voxel-wise vectors
z <- zmap[mask_idx]
p <- pmap[mask_idx]

# sorting once
ord <- order(p)
sortp <- p[ord]

# rank mapping
ranks <- integer(m)
ranks[ord] <- seq_len(m)

# index mapping (0-based for C++) for unsorted p-values
indexp <- as.integer(mask_idx - 1)

# mask rank volume
maskI <- integer(prod(dims))
maskI[mask_idx] <- seq_len(m)

cat("Building adjacency...\n")
adj <- findAdjList(maskI, indexp, dims, m, opt$conn)

# -----------------------------
# observed
# -----------------------------
nperm <- opt$nperm
alpha <- opt$alpha

zval <- NULL
maxsize <- NULL
allsize <- NULL

if (!is.null(opt$k)) {
  zval <- numeric(nperm + 1)
  id <- FindIndex(opt$k, m, adj, indexp, ord, ranks)
  if (id > 0) zval[nperm + 1] <- zmap[id]
}

if (!is.null(opt$cft)) {
  maxsize <- numeric(nperm + 1)
  allsize <- numeric(nperm + 1)
  
  rnk <- which(sortp <= pnorm(opt$cft, lower.tail=FALSE))
  if (length(rnk) > 0) {
    rnk <- rnk[length(rnk)]
    maxsize[nperm + 1] <- FindSize(rnk, m, adj, ord, ranks)
  }
  allsize[nperm + 1] <- sum(z >= opt$cft)
}

# -------------------------------------------
# permutations
# -------------------------------------------
# create permutation matrix (sign flipping)
set.seed(nperm)
perms <- matrix(sample(c(-1,1), nperm*n, replace=TRUE), nrow=nperm, ncol=n)

# run permutations
cat("Running permutations:\n")

for (np in 1:nperm) {
  
  cat(np, "\n")
  
  sign_vec <- perms[np, ]
  
  # apply permutation only on in-mask voxels
  data_perm <- V[mask_idx,, drop=FALSE]
  data_perm <- data_perm * sign_vec
  
  mu <- rowMeans(data_perm)
  sdv <- sqrt(rowMeans(data_perm^2) - mu^2)
  
  valid <- sdv > 0
  
  t_perm <- rep(0, m)
  t_perm[valid] <- sqrt(n) * mu[valid] / sdv[valid]
  
  p_perm <- rep(1, m)
  p_perm[valid] <- pt(t_perm[valid], df=n-1, lower.tail=FALSE)
  
  z_perm <- rep(0, m)
  z_perm[valid] <- qnorm(p_perm[valid], lower.tail=FALSE)
  
  zmap_perm <- array(0, dim=dims)
  zmap_perm[mask_idx] <- z_perm
  
  # ranking once
  ord <- order(p_perm)
  sortp <- p_perm[ord]
  
  ranks <- integer(m)
  ranks[ord] <- 1:m
  
  # ---------------- FindIndex ----------------
  if (!is.null(opt$k)) {
    id <- FindIndex(opt$k, m, adj, indexp, ord, ranks)
    if (id > 0) zval[np] <- zmap_perm[id]
  }
  
  # ---------------- FindSize ----------------
  if (!is.null(opt$cft)) {
    rnk <- which(sortp <= pnorm(opt$cft, lower.tail=FALSE))
    if (length(rnk) > 0) {
      rnk <- rnk[length(rnk)]
      maxsize[np] <- FindSize(rnk, m, adj, ord, ranks)
    }
    allsize[np] <- sum(z_perm >= opt$cft)
  }
}

# -----------------------------
# output
# -----------------------------
out <- opt$out
if (!grepl("\\.txt$", out)) out <- paste0(out, ".txt")
cat("", file=out)  # clear file

cut <- ceiling((nperm + 1) * (1 - alpha))

cat("RESULTS\n", file=out, append=TRUE)
cat("========================\n", file=out, append=TRUE)

if (!is.null(opt$k)) {
  zthres <- sort(zval)[cut]
  cat("METHOD: K-based\n", file=out, append=TRUE)
  cat("k =", opt$k, "\n", file=out, append=TRUE)
  cat("alpha =", alpha, "\n", file=out, append=TRUE)
  cat("z-threshold =", zthres, "\n", file=out, append=TRUE)
  cat("\n", file=out, append=TRUE)
}

if (!is.null(opt$cft)) {
  sizethres <- sort(maxsize)[cut]
  allsthres <- sort(allsize)[cut]
  cat("METHOD: Cluster-forming\n", file=out, append=TRUE)
  cat("cft =", opt$cft, "\n", file=out, append=TRUE)
  cat("cluster-size threshold =", sizethres, "\n", file=out, append=TRUE)
  cat("all cluster-size threshold =", allsthres, "\n", file=out, append=TRUE)
}

cat("Done.\n")
