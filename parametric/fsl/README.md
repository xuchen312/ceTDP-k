# ceTDP-k-FSL

This folder provides a **FSL-based implementation** of the Stage 1 **ceTDP-k** procedure for computing cluster extent thresholds (k) using Random Field Theory (RFT).

It relies on standard FSL tools for estimating smoothness and performing cluster-level inference.

## Requirements

- FSL installed and configured
- Standard neuroimaging preprocessing results (e.g., GLM output, z-stat maps, or residual images)

## Workflow Overview

The procedure consists of two main steps:

### 1. Smoothness estimation

Smoothness can be estimated using either residual images (recommended) or z-stat maps.

#### Using residuals (recommended)

```bash
smoothest -r <residual_image.nii.gz> -m <brain_mask.nii.gz>
```

#### Using z-stat maps (less optimal)

```bash
smoothest -z <z_map.nii.gz> -m <brain_mask.nii.gz>
```

#### Output from `smoothest`

- `DLH`: determinant of the Lambda matrix (used in RFT)
- `VOLUME`: number of voxels in the mask 
- `RESELS`: number of resolution elements (smoothness-adjusted volume)

### 2. Cluster-level inference

Use FSL `cluster` to compute cluster extent thresholds.

```bash
cluster -i <z-map> -t <CDT> -p <alpha> --dlh=<DLH> --volume=<VOLUME> --minclustersize --no_table
cluster \
  -i <z_map.nii.gz> \
  -t <CDT> \
  -p <alpha> \
  --dlh=<DLH> \
  --volume=<VOLUME> \
  --minclustersize \
  --no_table
```

## Inputs

- z-map (statistical map from GLM)
- brain mask
- CDT (cluster-defining threshold)
- alpha (significance level)
- smoothness estimates (DLH, RESELS, VOLUME)

## Output

- Cluster extent threshold (k) used in ceTDP-k Stage 1

This threshold is passed to Stage 2 for TDP estimation.

## Notes

- This is a pipeline wrapper around standard FSL tools
- It does not modify FSL statistical models
- Results depend on accurate smoothness estimation
