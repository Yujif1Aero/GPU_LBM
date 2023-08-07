#ifndef _MAIN_C_
#define _MAIN_C_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
using namespace std;

#include "boundary.hpp"
#include "collision.hpp"
#include "initialization.hpp"
#include "lbm_model.hpp"
#include "stopwatch.h"
#include "streaming.hpp"
#include "utils.hpp"
#include "visualization.hpp"
#define BCINITIAL_CHK

int main(int argc, char* argv[]) {
    float *collide_field = NULL, *stream_field = NULL, *swap = NULL, tau, nu, mlups_sum;
    int num_cells;

    /* Boundary condition parameter */
    float bc_velocity[D_LBM];
    float bc_pressure;
    int t = 0, timesteps, timesteps_per_plotting, gpu_enabled;

    /*GRID parameter*/
    char meshdir[250];
    int xmax, ymax, zmax, xstart, ystart, zstart, xend, yend, zend;
    float LengthRef;


    vector<double> xd, yd, zd;
    vector<int> bcd;

    /* INITIAL parameter*/
    float RhoInt, VelInt[D_LBM], PresInt, Feq[Q_LBM];

    clock_t mlups_time;
    size_t field_size;

    /* process parameters */
    ReadParameters(&LengthRef, &tau, &nu, VelInt, &RhoInt, &PresInt, bc_velocity, &bc_pressure, &timesteps, &timesteps_per_plotting, meshdir, argc, argv, &gpu_enabled);

    /* check if provided parameters are legitimate */
    printf("valida \n");
    ValidateModel(bc_velocity, LengthRef, tau, RhoInt);

    /* 	 grid generation */
    printf("Grid%s \n", meshdir);
    InitialiseGrid(meshdir, xmax, ymax, zmax, xstart, ystart, zstart, xend, yend, zend, xd, yd, zd, bcd);

    /* initializing fields */
    printf("Array initialisation \n");
    num_cells = xmax * ymax * zmax;
    field_size = Q_LBM * num_cells * sizeof(float);
    collide_field = (float*)malloc(field_size);
    stream_field = (float*)malloc(field_size);

    printf("Initialisation \n");
    InitialiseFields(collide_field, stream_field, xmax, ymax, zmax, gpu_enabled, Feq, RhoInt, VelInt, bcd);
    //WriteFluidVtkOutput(collide_field, "../Initial/initfluidfiled", t, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax, xd, yd, zd, bcd);
    WriteAllVtkOutput(collide_field, "./Initial/initallfiled", t, xmax, ymax, zmax, xd, yd, zd);
    WritePhysics(collide_field, "./Initial/physics-filed", t, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax, gpu_enabled);
    writebcd(bcd, "./Initial/init", xmax, ymax, zmax, gpu_enabled);
  

#ifdef BCINITIAL_CHK
    /* Treat boundaries */
    printf("CHK Boundary condition \n");
    TreatBoundary(collide_field, bc_velocity, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax, bcd, RhoInt);
    printf("END CHK BC \n");
        // WriteFluidVtkOutput(collide_field, "../Initial/initfluidfiled_bc", t, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax, xd, yd, zd, bcd);
    printf("Output Vtk for CHK BC \n");
    WriteAllVtkOutput(collide_field, "./Initial/initallfiled_bc", t, xmax, ymax, zmax, xd, yd, zd);
    WritePhysics(collide_field, "./Initial/physics-filed_bc", t, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax, gpu_enabled);
    printf("END output Vtk for CHK BC \n");
    //exit(0);
#endif
    mlups_sum = 0.f;
    for (t = 0; t <= timesteps; t++) {
        printf("Time step: #%d\r", t);
        // printf("Time step: #%d\n", t);
        Stopwatch sw;
        sw.start();
        /* Compute post collision distributions */
        DoCollision(collide_field, tau, xmax, ymax, zmax, xstart, ystart, zstart, xend, yend, zend, bcd);
        /* Copy pdfs from neighbouring cells into collide field */
        DoStreaming(collide_field, stream_field, xmax, ymax, zmax, xstart, ystart, zstart, xend, yend, zend, bcd);
        // /* Perform the swapping of collide and stream fields */
        swap = collide_field;
        collide_field = stream_field;
        stream_field = swap;
        // /* Treat boundaries */
        //cout << " Treat boundary t=" << t << endl;
        TreatBoundary(collide_field, bc_velocity, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax, bcd, RhoInt);

        sw.stop();
        mlups_time = sw.getMs();
        /* Print out the MLUPS value */
        mlups_sum += num_cells / (MLUPS_EXPONENT * (float)mlups_time);
        // if (VERBOSE)
        // printf("MLUPS: %f\n", num_cells / (MLUPS_EXPONENT * (float)mlups_time));
        /* Print out vtk output if needed */
        if (!(t % timesteps_per_plotting)) {
            WriteAllVtkOutput(collide_field, "./img/all-fielda", t, xmax, ymax, zmax, xd, yd, zd);
           // WriteFluidVtkOutput(collide_field, "../img/fluid-field", t, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax, xd, yd, zd, bcd);
        }
    }

    printf("Average MLUPS: %f\n", mlups_sum / (t + 1));

    // if (VERBOSE) {
    // 	WriteField(collide_field, "../Final/disfunc-field", t, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax,
    // 		gpu_enabled);
    // 	WritePhysics(collide_field, "../Final/physics-filed", t, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax, gpu_enabled);

    // }

    /* Free memory */
    free(collide_field);
    free(stream_field);
    xd.clear();
    yd.clear();
    zd.clear();
    bcd.clear();

    // FreeDeviceFields(&collide_field_d, &stream_field_d, &flag_field_d);

    printf("Simulation complete.\n");
    return 0;
}

#endif
