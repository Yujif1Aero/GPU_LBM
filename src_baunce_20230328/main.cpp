#ifndef _MAIN_C_
#define _MAIN_C_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>

using namespace std;
#include "boundary.hpp"
#include "collision.hpp"
#include "initialization.hpp"
#include "lbm_model.hpp"
#include "streaming.hpp"
#include "utils.hpp"
#include "visualization.hpp"

int main(int argc, char* argv[]) {
    float *collide_field = NULL, *stream_field = NULL, *swap = NULL, tau, wall_velocity[D_LBM], num_cells, mlups_sum;
    int xlength, t = 0, timesteps, timesteps_per_plotting, gpu_enabled;
    /*GRID parameter*/
    int xstep, ystep, zstep, xstart, ystart, zstart, xend, yend, zend;
    /* INITIAL parameter*/
    float RhoInt = 1.0, VelInt = 0.0, Feq[Q_LBM];

    clock_t mlups_time;
    size_t field_size;

    /* process parameters */
    ReadParameters(&xlength, &tau, VelInt, wall_velocity, &timesteps,
                   &timesteps_per_plotting, argc, argv, &gpu_enabled);

    /* check if provided parameters are legitimate */
    ValidateModel(wall_velocity, xlength, tau, VelInt, RhoInt);

    // grid generation //
    InitialiseGrid(xlength, xstep, ystep, zstep, xstart, ystart, zstart, xend, yend, zend);

    /* initializing fields */
    // num_cells = pow(xlength + 2, D_LBM);
    num_cells = xstep * ystep * zstep;
    field_size = Q_LBM * num_cells * sizeof(float);
    collide_field = (float*)malloc(field_size);
    stream_field = (float*)malloc(field_size);
    // flag_field = (int*)malloc(num_cells * sizeof(int));
    InitialiseFields(collide_field, stream_field, xstep, ystep, zstep, gpu_enabled, Feq, RhoInt, wall_velocity);
    // writeFlagField(flag_field, "./Initial/init", xstep, ystep, zstep, gpu_enabled);
    // TreatBoundary(collide_field, wall_velocity, xstart, ystart,
    //               zstart, xend, yend, zend, xstep, ystep, zstep);

    WriteFluidVtkOutput(collide_field, "./Initial/initfluidfiled", t, xstart, ystart,
                        zstart, xend, yend, zend, xstep, ystep, zstep);
    //WriteAllVtkOutput(collide_field, "./Initial/initallfiled", t, xstep, ystep, zstep);

    WriteField(collide_field, "./Initial/field", t, xstart, ystart, zstart, xend, yend, zend, xstep, ystep, zstep,
               gpu_enabled);
    WritePhysics(collide_field, "Initial/physics-filed", t, xstart, ystart,
                 zstart, xend, yend, zend, xstep, ystep, zstep, gpu_enabled);
    for (t = 0; t <= timesteps; t++) {
        printf("Time step: #%d\n", t);
        mlups_time = clock();
        /* Compute post collision distributions */
        DoCollision(collide_field, tau, xstep, ystep, zstep, xstart, ystart, zstart, xend, yend, zend);
        /* Copy pdfs from neighbouring cells into collide field */
        DoStreaming(collide_field, stream_field, xstep, ystep, zstep, xstart, ystart, zstart, xend, yend, zend);
        /* Perform the swapping of collide and stream fields */
        swap = collide_field;
        collide_field = stream_field;
        stream_field = swap;
        /* Treat boundaries */
        TreatBoundary(collide_field, wall_velocity, xstart, ystart, zstart, xend, yend, zend, xstep, ystep, zstep);
        mlups_time = clock() - mlups_time;
        /* Print out the MLUPS value */
        mlups_sum += num_cells / (MLUPS_EXPONENT * (float)mlups_time / CLOCKS_PER_SEC);
        //if (VERBOSE)
            printf("MLUPS: %f\n", num_cells / (MLUPS_EXPONENT * (float)mlups_time / CLOCKS_PER_SEC));
        /* Print out vtk output if needed */
        if (!(t % timesteps_per_plotting)) {
            WriteFluidVtkOutput(collide_field, "img/fluid-field", t, xstart + 1, ystart + 1, zstart, xend - 1, yend - 1, zend, xstep, ystep, zstep);
        }
    }

    printf("Average MLUPS: %f\n", mlups_sum / (t + 1));

    if (VERBOSE) {
        WriteField(collide_field, "Final/disfunc-field", t, xstart, ystart, zstart, xend, yend, zend, xstep, ystep, zstep,
                   gpu_enabled);
        WritePhysics(collide_field, "Final/physics-filed", t, xstart + 1, ystart + 1, zstart, xend - 1, yend - 1, zend, xstep, ystep, zstep, gpu_enabled);

    }

    /* Free memory */
    free(collide_field);
    free(stream_field);

    // FreeDeviceFields(&collide_field_d, &stream_field_d, &flag_field_d);

    printf("Simulation complete.\n");
    return 0;
}

#endif
