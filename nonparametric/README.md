# ceTDP-k Nonparametric Implementation (Stage 1)

This directory provides a **nonparametric implementation** of Stage 1 (ceTDP-k) using permutation-based inference.

It estimates minimum significant cluster extent threshold (`k`) without relying on Random Field Theory (RFT).

## Method Overview

This implementation uses a **permutation-based cluster inference framework**:

- Voxel-wise test statistics are computed from input images
- A null distribution is generated using sign-flipping permutations
- Cluster formation is performed using voxel adjacency graphs (C++ backend)
- Cluster extent thresholds (`k`) are estimated from empirical null distributions

## Requirements

- R ($\ge$ 4.0 recommended)
- Packages:
  - `optparse`
  - `RNifti`
  - `Rcpp`

Dependencies are automatically installed if missing.

## Usage

Run the main script:

```bash
Rscript kval.R \
  --dirpath=/data/imgs \
  --mask=/data/mask.nii.gz \
  --conn=18 \
  --nperm=1000 \
  --alpha=0.05 \
  --kval=14 \
  --cft=3.1 \
  --out=results.txt
```

## Inputs

- `dirpath`: Input image directory (COPE or statistic maps)
- `mask`: Brain mask (NIfTI)
- `conn`: Voxel connectivity (6 / 18 / 26)
- `nperm`: Number of permutations
- `alpha`: Significance level
- `kval`: Optional cluster extent threshold mode
- `cft`: Cluster-forming threshold
- `out`: Output file

## Outputs

The script outputs:

- Cluster extent threshold (`k`)
- Optional z-threshold or cluster-size statistics

## Notes

- This implementation corresponds to the nonparametric branch of Stage 1 (ceTDP-k).
- The resulting cluster extent threshold (`k`) is used as input for Stage 2 (ceTDP-cons & ceTDP-appr).
- No RFT assumptions are used.
- Results depend on number of permutations (`nperm`).
