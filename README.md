# ceTDP-k: Cluster Extent Threshold Computation (Stage 1) of ceTDP

**ceTDP-k** is the Stage 1 component of the ceTDP framework for cluster-extent inference in neuroimaging.

It computes the minimum significant cluster extent thresholds (*k*) that are required as input for Stage 2, which estimates the True Discovery Proportion (TDP) within chosen regions.

It can also help users to prepare the input cluster information files (CSV) for Stage 2.

## Overview

Cluster-extent inference using Random Field Theory (RFT) is widely used in fMRI analysis, but typically provides limited information about the spatial distribution of activation within clusters.

The ceTDP framework addresses this by:

- **Stage 1**: computing minimal cluster extent thresholds (k)
- **Stage 2**: estimating TDP bounds within any regions

This repository implements Stage 1.

## Methods

Stage 1 supports two main approaches:

### Parametric approach

Based on Random Field Theory (RFT), implemented using:
- SPM (MATLAB-based implementation)
- FSL-based pipeline (wrapper using FSL tools)

### Nonparametric approach

Based on nonparametric inference, implemented in R and C++

## Repository Structure

```text
ceTDP-k/
│
├── parametric/
│   ├── spm/        # MATLAB implementation (SPM-based)
│   ├── fsl/        # FSL pipeline (workflow-based, see README.md)
│
├── nonparametric/  # R & C++ implementation
|
└── README.md
```

