# ceTDP-k SPM Implementation

This directory contains the **SPM-based implementation** of ceTDP-k, the cluster extent threshold computation component of Stage 1.

## Method Overview

This implementation is a MATLAB function built on top of **SPM12 cluster-level inference**.

It computes the minimum significant cluster extent threshold (`k`) based on the selected SPM cluster-forming threshold and cluster-level inference settings.

Key characteristics:

- Based on Gaussian Random Field Theory (RFT)
- Integrated with SPM12 statistical results (`SPM` / `xSPM`)
- Interactive SPM GUI workflow
- No standalone statistical model outside SPM

## Main Function

The main function is:

```matlab
spm_ceTDP_minsize.m
```

## Usage

### Recommended usage (interactive mode)

Run the function without inputs:

```matlab
spm_ceTDP_minsize
```

Typical workflow:

1. Select an SPM results structure (`SPM`) via the SPM GUI.
2. Specify the desired thresholding and cluster-level inference options.
3. Run the computation.

## Output

The function computes:

> **k (cluster extent threshold in voxels)**

This is the minimum cluster extent, in voxels, required for cluster-level significance under the selected SPM configuration.

## Additional Output (Reporting Only)

The function also displays standard SPM statistical information, including:

- Cluster-forming threshold (T / P)
- FWE-corrected significance level
- Cluster-level corrected p-values

These values are provided for reporting purposes and are **not required as Stage 2 inputs**.

Example output:

```text
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

## Integration with ceTDP

This implementation provides the **SPM-based parametric computation of ceTDP-k**.

The resulting cluster extent threshold (`k`) is used in the subsequent Stage 1 workflow for region specification and cluster information preparation, and is required as an input for Stage 2 (ceTDP-cons / ceTDP-appr).

## Dependencies

- MATLAB
- SPM12 (must be added to the MATLAB path)

Example:

```matlab
addpath(genpath('/path/to/spm12'));
```

## Notes

- This implementation does not replace SPM cluster-level inference.
- It operates within SPM's statistical framework.
- Results depend on SPM model specification and thresholding choices.
