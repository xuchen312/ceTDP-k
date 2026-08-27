# ceTDP-k FSL Implementation

This directory contains the **FSL-based implementation of ceTDP-k**, the cluster extent threshold computation component of Stage 1 in the ceTDP framework.

It uses standard FSL tools to estimate image smoothness and perform Random Field Theory (RFT)-based cluster-level inference.

## Requirements

- FSL installed and configured
- Z-statistic map
- Brain mask
- Residual image (recommended for smoothness estimation)

## Workflow Overview

The procedure consists of two main steps:

### 1. Smoothness estimation

Image smoothness can be estimated from either residual images (recommended) or Z-statistic maps.

#### Using residuals (recommended)

```bash
smoothest -r <residual_image.nii.gz> -m <brain_mask.nii.gz>
```

#### Using a Z-statistic map (less optimal)

```bash
smoothest -z <z_map.nii.gz> -m <brain_mask.nii.gz>
```

The `smoothest` command reports several quantities, including:

- **DLH** — determinant of the Lambda matrix used for RFT-based inference
- **VOLUME** — number of voxels within the analysis mask
- **RESELS** — number of resolution elements (smoothness-adjusted search volume)

For the `cluster` command below, **DLH** and **VOLUME** are used as inputs.

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

The `--minclustersize` option returns the minimum cluster extent threshold corresponding to the specified cluster-defining threshold (`CDT`) and significance level (`alpha`).

## Inputs

The procedure requires:

- Z-statistic map
- Brain mask
- Cluster-defining threshold (`CDT`)
- Significance level (`alpha`)
- Smoothness estimates (`DLH` and `VOLUME`)

A residual image is recommended for estimating smoothness.

## Output

> **k (minimum significant cluster extent threshold, in voxels)**

The resulting `k` is used in the subsequent Stage 1 workflow for region specification and cluster information preparation, and is required as an input to Stage 2 of ceTDP.

## Notes

- This implementation uses standard FSL tools for RFT-based cluster-level inference.
- It does not modify FSL statistical models or inference procedures.
- Accurate smoothness estimation is important for valid RFT-based inference.
- Results depend on the statistical map, analysis mask, cluster-defining threshold, significance level, and estimated smoothness.
