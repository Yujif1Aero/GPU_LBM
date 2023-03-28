# GPU_LBM
This is LBM programs of which the development and research are progress.

March 28 2023
--Launching this github
    Lid-driven cavity problem benchmark with half-way bounce back scheme and bounce back scheme boundary condition.
    
    You can calculate this LBM with half-way bounce back scheme( Default )
    $./bin/lbm-sim data/lbm.dat -cpu
    If you can calculate this LBM with bounce back scheme, switch a name directory of src_baunce_20230328 to src
    Initial condition is in ./data/lbm.data
    Output data are in ./img

<img width="476" alt="image" src="https://user-images.githubusercontent.com/116667889/228243975-9bdc9a1b-8b14-4c70-b39d-f5860cf077aa.png">
