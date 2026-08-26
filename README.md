# ceTDP-k: Cluster Extent Threshold Computation in ceTDP Stage 1

**ceTDP-k** is the **cluster extent threshold computation component of Stage 1** in the ceTDP framework for cluster-extent inference in neuroimaging.

It computes the minimum significant cluster extent threshold (`k`) required for region specification and Stage 2 of ceTDP, where `k` is used for **True Discovery Proportion (TDP)** inference within selected regions.

This repository provides parametric and nonparametric implementations for estimating `k` using Random Field Theory (RFT)-based and permutation-based approaches.

## Overview

Cluster-extent inference using RFT is widely used in fMRI analysis, but typically provides limited information about the spatial localization of activation within significant clusters.

The ceTDP framework addresses this by:

- **Stage 1**: computing the minimum significant cluster extent threshold (`k`) and specifying regions for Stage 2 inference
- **Stage 2**: estimating TDP bounds within selected regions

This repository implements the cluster extent threshold computation component of Stage 1.

## Methods

Stage 1 supports two main approaches:

---

### 1. Parametric approach

Based on RFT, implemented using:

#### SPM (MATLAB-based implementation)

The SPM implementation consists of a single function:

`spm_ceTDP_minsize.m`

This function is designed for interactive use within the SPM GUI environment.

It computes the minimum significant cluster extent threshold (`k`) based on SPM's RFT-based inference.

The primary output is:

> **k (cluster extent threshold in voxels)**

which is used as an input for Stage 2 (ceTDP-cons & ceTDP-appr).

See `parametric/spm/README.md` for implementation details.

#### FSL-based implementation

The FSL implementation uses a standard FSL pipeline for RFT-based cluster inference and obtains the corresponding minimum significant cluster extent threshold (`k`).

See `parametric/fsl/README.md`.

---

### 2. Nonparametric approach

The nonparametric implementation estimates the minimum significant cluster extent threshold (`k`) using permutation-based inference, without relying on RFT assumptions.

It is implemented in R and C++.

See `nonparametric/README.md`.

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

