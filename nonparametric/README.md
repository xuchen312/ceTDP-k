# ceTDP-k Nonparametric Implementation

This directory contains the **nonparametric implementation of ceTDP-k**, the cluster extent threshold computation component of Stage 1 in the ceTDP framework.

It estimates the minimum significant cluster extent threshold (`k`) using permutation-based inference without relying on Random Field Theory (RFT). The implementation also supports a reverse thresholding mode in which the cluster-forming threshold (CFT) is determined for a given cluster extent threshold (`k`).

## Method Overview

This implementation uses a permutation-based cluster inference framework:

- Voxel-wise test statistics are computed from the input images.
- A null distribution is generated using sign-flipping permutations.
- Clusters are identified using disjoint-set forest implemented in C++.
- Cluster-level thresholds are estimated from the empirical permutation distribution.

Unlike the parametric implementations, this approach does not rely on RFT assumptions.

## Inference Modes

The implementation supports two complementary thresholding modes.

### Mode 1: CFT-to-`k`

Given a cluster-forming threshold (`cft`) and a significance level (`alpha`), the method estimates the minimum significant cluster extent threshold (`k`) from the empirical permutation distribution.

- **Input:** CFT (`cft`)
- **Output:** Minimum significant cluster extent

### Mode 2: `k`-to-CFT

Given a cluster extent threshold (`kval`) and a significance level (`alpha`), the method determines the corresponding cluster-forming threshold (`cft`) from the empirical permutation distribution.

- **Input:** Cluster extent threshold (`kval`)
- **Output:** CFT

Thus, the implementation supports both **CFT-to-`k`** and **`k`-to-CFT** thresholding modes.

## Requirements

- R (>=4.0 recommended)
- Packages:
  - `optparse`
  - `RNifti`
  - `Rcpp`

Required packages are installed automatically if missing.

## Usage

### Mode 1: CFT-to-`k`

To estimate the minimum significant cluster extent threshold for a given CFT:

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

The script computes the minimum significant cluster extent threshold (`k`) corresponding to the specified CFT.

### Mode 2: `k`-to-CFT

To determine the CFT corresponding to a given cluster extent threshold:

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

The script estimates the CFT (`cft`) corresponding to the specified cluster extent threshold.

## Inputs

| Argument | Description |
|----------|-------------|
| `dirpath` | Directory containing the input images (e.g., COPE or statistic maps) |
| `mask` | Brain mask in NIfTI format |
| `conn` | Voxel connectivity used for cluster formation (6, 18, or 26) |
| `nperm` | Number of permutations |
| `alpha` | Significance level |
| `cft` | Cluster-forming threshold. Used in **CFT-to-`k`** mode |
| `kval` | Cluster extent threshold. Used in **`k`-to-CFT** mode |
| `out` | Output file |

The inference mode is determined by the threshold parameter provided:

- Providing `cft` selects **CFT-to-`k`** mode.
- Providing `kval` selects **`k`-to-CFT** mode.

## Outputs

### Mode 1: CFT-to-`k`

The primary output is:

> **k (minimum significant cluster extent threshold, in voxels)**

### Mode 2: `k`-to-CFT

The primary output is:

> **cft (cluster-forming threshold)**

Additional permutation-based statistics may also be reported in the output file.

## Integration with ceTDP

This implementation is the **nonparametric branch of ceTDP-k**, providing the cluster extent threshold computation component of Stage 1.

The resulting cluster extent threshold (`k`) can be used in the subsequent Stage 1 workflow for region specification and cluster information preparation, and is required as an input to Stage 2 (ceTDP-cons / ceTDP-appr).

## Notes

- No RFT assumptions are required.
- Results depend on the number of permutations (`nperm`).
- Voxel connectivity (`conn`) determines how clusters are formed.
- In **CFT-to-`k`** mode, the choice of CFT affects the resulting cluster extent threshold (`k`).
- In **`k`-to-CFT** mode, the specified cluster extent threshold (`kval`) determines the corresponding CFT.
