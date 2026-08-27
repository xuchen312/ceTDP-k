# ceTDP-k: Cluster Extent Threshold Computation in ceTDP Stage 1

**ceTDP-k** is the **cluster extent threshold computation component of Stage 1** in the ceTDP framework for cluster-extent inference in neuroimaging.

It computes the minimum significant cluster extent threshold (`k`) used in subsequent region specification and Stage 2 of ceTDP, where `k` is used for **True Discovery Proportion (TDP)** inference within selected regions.

This repository provides parametric and nonparametric implementations for estimating `k` using Random Field Theory (RFT)-based and permutation-based approaches.

## Overview

Cluster-extent inference using RFT is widely used in fMRI analysis, but typically provides limited information about the spatial localization of activation within significant clusters.

The ceTDP framework addresses this by:

- **Stage 1**: computing the minimum significant cluster extent threshold (`k`) and specifying regions for subsequent inference
- **Stage 2**: estimating TDP bounds within selected regions

This repository implements **the cluster extent threshold computation component of Stage 1 (ceTDP-k)**. Region specification is a separate part of Stage 1 and is not implemented in this repository.

## Methods

This repository provides two approaches for cluster extent threshold computation:

---

### 1. Parametric approach

The parametric approach is based on **Gaussian Random Field Theory** and provides implementations using SPM and FSL.

#### SPM (MATLAB-based implementation)

The SPM implementation consists of a single function:

`spm_ceTDP_minsize.m`

This function is designed for interactive use within the SPM GUI environment.

It computes the minimum significant cluster extent threshold (`k`) based on SPM's RFT-based cluster-level inference.

The primary output is:

> **k (cluster extent threshold in voxels)**

This value is used in subsequent Stage 1 processing and as an input to Stage 2 (ceTDP-cons and ceTDP-appr).

See `parametric/spm/README.md` for implementation details.

#### FSL-based implementation

The FSL implementation uses standard FSL tools for RFT-based cluster inference and obtains the corresponding minimum significant cluster extent threshold (`k`).

See `parametric/fsl/README.md` for implementation details.

---

### 2. Nonparametric approach

The nonparametric implementation estimates the minimum significant cluster extent threshold (`k`) using permutation-based inference without relying on RFT assumptions.

It is implemented in R and C++.

See `nonparametric/README.md` for implementation details.

## Output

The primary output of all implementations is:

> **k (minimum significant cluster extent threshold in voxels)**

The resulting `k` is used by the subsequent ceTDP workflow for region specification and Stage 2 TDP inference.

## Repository Structure

```text
ceTDP-k/
│
├── parametric/
│   ├── spm/        # MATLAB SPM-based implementation
│   └── fsl/        # FSL-based pipeline implementation
│
├── nonparametric/  # R & C++ permutation-based implementation
│
└── README.md
```
