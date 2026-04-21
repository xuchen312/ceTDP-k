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

#### SPM (MATLAB-based implementation)

The SPM implementation consists of a single function:

`spm_ceTDP_minsize.m`

This function is designed for **interactive use within the SPM GUI environment**.

Typical usage:

- Run the function without inputs
- Select an SPM results file via the SPM interface
- Choose the desired thresholding options within the GUI

The function then computes the minimum significant cluster extent threshold (k) based on SPM's RFT-based cluster inference.

The primary output is:

> **k (cluster extent threshold in voxels)**

which is used for Stage 2 (ceTDP-cons & ceTDP-appr).

The function also prints standard SPM thresholding information (e.g., T/P thresholds), which are provided for reporting purposes only.

```matlab
====================================================
 Cluster extent threshold (FWE corrected)
====================================================

 Threshold:
   P = 0.001  (T = 3.3962)
   alpha = 0.050

 Result:
   k = 129 voxels
   corrected P = 0.0497
```

#### FSL-based implementation

Follows a pipeline using standard FSL tools (see parametric/fsl/README.md).

### Nonparametric approach

Based on nonparametric inference, implemented in R and C++

## Repository Structure

```text
ceTDP-k/
│
├── parametric/
│   ├── spm/        # MATLAB SPM-based single-function implementation (GUI + console output)
│   ├── fsl/        # FSL-based pipeline implementation (see README.md)
│
├── nonparametric/  # R & C++ implementation (nonparametric inference)
|
└── README.md
```

