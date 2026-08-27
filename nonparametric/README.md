# ceTDP-k Nonparametric Implementation

This directory contains the **nonparametric implementation of ceTDP-k**, the cluster extent threshold computation component of Stage 1 in the ceTDP framework.

It estimates the minimum significant cluster extent threshold (`k`) using permutation-based inference without relying on Random Field Theory (RFT). The implementation also supports the reverse thresholding operation, in which the cluster-forming threshold (CFT) is determined for a given cluster extent threshold (`k`).

## Method Overview

This implementation uses a permutation-based cluster inference framework:

- Voxel-wise test statistics are computed from the input images.
- A null distribution is generated using sign-flipping permutations.
- Clusters are identified using voxel adjacency graphs implemented in C++.
- Cluster-level thresholds are estimated from the empirical permutation distribution.

Unlike the parametric implementations, this approach does not rely on RFT assumptions.

## Thresholding Operations

The implementation supports two complementary thresholding operations.

### 1. CFT-to-k

Given a cluster-forming threshold (`cft`) and a significance level (`alpha`), the method estimates the minimum significant cluster extent threshold from the empirical permutation distribution.

- **Input:** CFT (`cft`)
- **Output:** Minimum significant cluster extent

### 2. k-to-CFT

Given a cluster extent threshold (`kval`) and a significance level (`alpha`), the method determines the corresponding cluster-forming threshold from the empirical permutation distribution.

- **Input:** Cluster extent threshold (`kval`)
- **Output:** CFT

The two operations can be run independently or together in a single execution by providing the corresponding arguments.

## Requirements

- R (>= 4.0 recommended)
- Packages:
  - `optparse`
  - `RNifti`
  - `Rcpp`

Required packages are installed automatically if missing.

## Usage

The main script is:

```bash
Rscript kval.R
```

Both short and long command-line options are supported.

### Mode 1: CFT-to-k

To estimate the minimum significant cluster extent threshold for a given CFT:

**Long-option form**

```bash
Rscript kval.R \
  --dirpath=/data/imgs \
  --mask=/data/mask.nii.gz \
  --conn=18 \
  --nperm=1000 \
  --alpha=0.05 \
  --cft=3.1 \
  --out=results.txt
```

**Short-option form**

```bash
Rscript kval.R \
  -d /data/imgs \
  -m /data/mask.nii.gz \
  -c 18 \
  -n 1000 \
  -a 0.05 \
  -t 3.1 \
  -o results.txt
```

The script reports the minimum significant cluster extent threshold corresponding to the specified CFT.

### Mode 2: k-to-CFT

To determine the CFT corresponding to a given cluster extent threshold:

**Long-option form**

```bash
Rscript kval.R \
  --dirpath=/data/imgs \
  --mask=/data/mask.nii.gz \
  --conn=18 \
  --nperm=1000 \
  --alpha=0.05 \
  --kval=14 \
  --out=results.txt
```

**Short-option form**

```bash
Rscript kval.R \
  -d /data/imgs \
  -m /data/mask.nii.gz \
  -c 18 \
  -n 1000 \
  -a 0.05 \
  -k 14 \
  -o results.txt
```

The script reports the CFT corresponding to the specified cluster extent threshold.

### Running both operations

Both `cft` and `kval` can be provided in the same execution. In this case, the script computes both sets of results:

```bash
Rscript kval.R \
  -d /data/imgs \
  -m /data/mask.nii.gz \
  -c 18 \
  -n 1000 \
  -a 0.05 \
  -k 14 \
  -t 3.1 \
  -o results.txt
```

The output file will contain the results for both the `k`-based and CFT-based calculations.

## Command-Line Arguments

| Short | Long | Description |
|------|------|-------------|
| `-d` | `--dirpath` | Directory containing the input images (e.g., COPE or statistic maps) |
| `-m` | `--mask` | Brain mask in NIfTI format |
| `-c` | `--conn` | Voxel connectivity used for cluster formation: 6, 18, or 26 |
| `-n` | `--nperm` | Number of sign-flipping permutations |
| `-a` | `--alpha` | Significance level |
| `-k` | `--kval` | Cluster extent threshold used for k-to-CFT calculation |
| `-t` | `--cft` | Cluster-forming threshold used for CFT-to-k calculation |
| `-o` | `--out` | Output file; `.txt` is added automatically if omitted |

### Required arguments

The following arguments are required:

- `--dirpath` / `-d`
- `--mask` / `-m`

The following threshold arguments are optional individually:

- `--cft` / `-t`
- `--kval` / `-k`

At least one of `cft` or `kval` must be provided.

If both are provided, both calculations are performed.

## Inputs

### Input images

The `dirpath` argument should point to a directory containing the subject-level input images.

For example:

```text
/data/imgs/
├── cope01.nii.gz
├── cope02.nii.gz
├── cope03.nii.gz
└── ...
```

The images should have matching spatial dimensions and correspond to the same analysis.

### Brain mask

The `mask` argument specifies a NIfTI brain mask defining the analysis volume.

### Voxel connectivity

The `conn` argument specifies the voxel adjacency criterion used to define connected clusters:

- `6`
- `18`
- `26`

The default is `18`.

### Permutations

The `nperm` argument controls the number of sign-flipping permutations used to construct the empirical null distribution.

A larger number of permutations generally provides a more stable estimate of the empirical tail distribution but increases computation time.

### Significance level

The `alpha` argument specifies the significance level used to determine the empirical threshold.

The default is:

```text
alpha = 0.05
```

## Outputs

Results are written to the specified output text file.

### CFT-to-k

When `cft` is provided, the output includes:

```text
METHOD: Cluster-forming
cft = ...
cluster-size threshold = ...
all cluster-size threshold = ...
```

The primary output for ceTDP-k is:

> **cluster-size threshold (`k`)**

where `k` denotes the minimum cluster size required for cluster-level significance under the specified CFT and significance level.

The total supra-threshold voxel threshold is an additional permutation-based summary statistic and is not used as the ceTDP-k threshold `k`.

### k-to-CFT

When `kval` is provided, the output includes:

```text
METHOD: K-based
k = ...
alpha = ...
z-threshold = ...
```

Here, `k` specifies the cluster extent threshold, and the resulting z-threshold is the corresponding CFT.

### Both calculations

If both `cft` and `kval` are provided, the output file contains both result sections.

## Integration with ceTDP

This implementation is the **nonparametric branch of ceTDP-k**, providing the cluster extent threshold computation component of Stage 1.

The resulting cluster extent threshold can be used in the subsequent Stage 1 workflow for region specification and cluster information preparation, and is required as an input to Stage 2 (ceTDP-cons / ceTDP-appr).

The cluster information required by Stage 2 can be prepared using the `cluster-export` utilities provided in the main ceTDP software suite.

## Notes

- No RFT assumptions are required.
- The null distribution is generated using sign-flipping permutations.
- Results depend on the number of permutations (`nperm`).
- Voxel connectivity (`conn`) determines how clusters are formed.
- In CFT-to-k calculation, the choice of CFT affects the resulting cluster extent threshold.
- In k-to-CFT calculation, the specified cluster extent threshold determines the corresponding CFT.
- The input images must have compatible spatial dimensions.
- The analysis mask defines the search volume used for cluster inference.
- Computation time increases with the number of permutations and the size of the analysis volume.
