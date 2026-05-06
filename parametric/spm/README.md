# ceTDP-k SPM Implementation (Stage 1)

This directory contains the **SPM-based implementation** of Stage 1 (ceTDP-k), which computes the cluster extent threshold using Random Field Theory (RFT).

## Method Overview

This implementation is a MATLAB wrapper built on top of **SPM12 cluster inference**.

It derives the minimum significant cluster extent threshold (`k`) based on a user-defined cluster-forming threshold.

Key characteristics:

- Based on Gaussian RFT
- Fully integrated with SPM12 results structure (`xSPM`)
- Interactive GUI-based workflow
- No standalone statistical model outside SPM

## Main Function

The implementation consists of a single MATLAB function:

```matlab
**spm_ceTDP_minsize.m**
```

## Usage

### Recommended usage (interactive mode)

Run the function without inputs:

```matlab
spm_ceTDP_minsize
```

Typical workflow:

- Select an SPM results structure (`SPM`) via the SPM GUI
- Choose the desired thresholding options within the GUI
- Run the computation

## Output

The function computes:

**k (cluster extent threshold in voxels)**

This is the minimum cluster size required for significance under the selected SPM configuration.

## Additional Output (Reporting Only)

The function also displays standard SPM statistical information, including:

- Cluster-forming threshold (T / P)
- FWE-corrected significance level
- Cluster-level corrected p-values

These values are provided for reporting purposes only and are **not used as Stage 2 inputs**.

Example outputs:

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

## Integration with ceTDP-k

This implementation is part of the **parametric branch of Stage 1 (ceTDP-k)**.

The resulting cluster extent threshold (**k**) is used as input for Stage 2 (ceTDP-cons / ceTDP-appr).

## Dependencies

- MATLAB
- SPM12 (must be added to MATLAB path)

Example:

```matlab
addpath(genpath('/path/to/spm12'));
```

## Notes

- This implementation does not replace SPM inference
- It operates fully within SPM’s statistical framework
- Results depend on SPM model specification and thresholding choices
