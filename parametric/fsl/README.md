# ceTDP-k-FSL

Require the installation of FSL.

1. Estimate smoothness

   - Estimate smoothness using a 4D residual image (e.g., `res4d.nii.gz`)
     ```bash
     smoothest -r <residual image> -m <mask>
     ```

   - Estimate smoothness from z-stat image (less ideal than using residuals)
     ```bash
     smoothest -z <z-map> -m <mask>
     ```

  Output:

  - `DLH`: determinant of the Lambda matrix (used in random field theory) 
  - `VOLUME`: number of voxels in mask 
  - `RESELS`: resolution elements (volume adjusted for smoothness)

2. Run cluster to compute threshold
   ```bash
   cluster -i <z-map> -t <CDT> -p <alpha> --dlh=<DLH> --volume=<VOLUME> --minclustersize --no_table
   ```
