# ceTDP-k FSL Implementation

This directory contains the **FSL-based implementation of ceTDP-k**, the cluster extent threshold computation component of Stage 1 in the ceTDP framework.

It relies on standard FSL tools for estimating image smoothness and performing Random Field Theory (RFT)-based cluster-extent inference.

## Requirements

- FSL installed and configured
- GLM results (e.g., z-stat maps)
- Brain mask
- Residual images (recommended for smoothness estimation)

## Workflow Overview

The procedure consists of two main steps:

### 1. Smoothness estimation

Smoothness can be estimated from either residual images (recommended) or z-stat maps.

#### Using residuals (recommended)

```bash
smoothest -r <residual_image.nii.gz> -m <brain_mask.nii.gz>
```

#### Using z-stat maps (less optimal)

```bash
smoothest -z <z_map.nii.gz> -m <brain_mask.nii.gz>
```

The `smoothest` command returns:

- **DLH** — determinant of the Lambda matrix used in RFT
- **VOLUME** — number of voxels within the analysis mask
- **RESELS** — number of resolution elements (smoothness-adjusted search volume)

### 2. Cluster-level inference

Use FSL `cluster` to obtain the minimum significant cluster extent threshold (`k`).

**Single-line command**

```bash
cluster -i <zmap.nii.gz> -t <CDT> -p <alpha> --dlh=<DLH> --volume=<VOLUME> --minclustersize --no_table
```
**Multi-line equivalent**

```bash
cluster \
  -i <zmap.nii.gz> \
  -t <CDT> \
  -p <alpha> \
  --dlh=<DLH> \
  --volume=<VOLUME> \
  --minclustersize \
  --no_table
```

## Inputs

- Z-statistic map
- Brain mask
- Cluster-defining threshold (CDT)
- Significance level (alpha)
- Smoothness estimates (DLH, RESELS and VOLUME)

## Output

> **k (minimum significant cluster extent threshold, in voxels)**

The resulting `k` is used for subsequent region specification in Stage 1 and as an input to Stage 2 of ceTDP.

## Notes

- This implementation is a wrapper around standard FSL RFT-based cluster inference.
- It does not modify FSL statistical models or inference procedures.
- Results depend on accurate smoothness estimation.
