#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <unistd.h>
using namespace std;

#include "initialization.hpp"
#include "lbm_model.hpp"
#include "utils.hpp"

/**
 * Prints help message that includes program usage instructions and control
 * flags.
 */
void PrintHelpMessage() {
    printf("List of control flags:\n");
    printf("\t -gpu             all computations are to be performed on gpu "
           "(now unable) \n");
    printf("\t -cpu             all computations are to be performed on cpu\n");
    printf("\t -help            prints this help message\n");
    printf("NOTE: Control flags are mutually exclusive and only one flag at a "
           "time is allowed\n");
    printf("Example program usage:\n");
    printf("\t./lbm-sim ./data/lbm.dat -cpu\n");
    exit(1);
}

void ReadParameters(int* xlength, float* tau, float& VelInt,
                    float* velocity_wall, int* timesteps,
                    int* timesteps_per_plotting, int argc, char* argv[],
                    int* gpu_enabled) {
    float *velocity_wall_1, *velocity_wall_2, *velocity_wall_3, re, nu;
    /* printing out help message */
    if (argc < 3)
        PrintHelpMessage();
    if (!strcmp(argv[1], "-help") || !strcmp(argv[2], "-help"))
        PrintHelpMessage();

    /* checking parameters */
    if (access(argv[1], R_OK) != 0)
        ERROR("Provided configuration file path either doesn't exist or can "
              "not be read.");
    *gpu_enabled = 0;

    /* reading parameters */
    // READ_FLOAT(argv[1], *tau);
    READ_FLOAT(argv[1], re);

    velocity_wall_1 = &velocity_wall[0];
    velocity_wall_2 = &velocity_wall[1];
    velocity_wall_3 = &velocity_wall[2];

    READ_FLOAT(argv[1], *velocity_wall_1);
    READ_FLOAT(argv[1], *velocity_wall_2);
    READ_FLOAT(argv[1], *velocity_wall_3);

    READ_INT(argv[1], *xlength);
    READ_INT(argv[1], *timesteps);
    READ_INT(argv[1], *timesteps_per_plotting);
    VelInt = sqrt((*velocity_wall_1) * (*velocity_wall_1) +
                  (*velocity_wall_2) * (*velocity_wall_2) +
                  (*velocity_wall_3) * (*velocity_wall_3));
    nu = VelInt * (float)(*xlength) / re;
    *tau = C_S_POW2_INV * nu / (C * C * DT) + 1.0 / 2.0;
}

void InitialiseFields(float* collide_field, float* stream_field, int xstep, int ystep, int zstep, int gpu_enabled, float Feq[Q_LBM], float RhoInt, float wall_velocity[D_LBM]) {
    int x, y, z, i;
    float dot_prod_cu, dot_prod_cu2, dot_prod_uu;
    float VelInt = 0.0; // only cavity
    float Velzero = 0.0;
    /* NOTE: We use y=ystep-1 as the moving wall */
    /* Initialization including dummy wall*/
    for (z = 0; z < zstep; z++) {
        for (y = 0; y < ystep; y++) {
            for (x = 0; x < xstep; x++) {

                for (i = 0; i < Q_LBM; i++) {

                    /* Initializing condition */
                    dot_prod_cu = LATTICE_VELOCITIES[i][0] * Velzero +
                                  LATTICE_VELOCITIES[i][1] * Velzero +
                                  LATTICE_VELOCITIES[i][2] * Velzero;
                    dot_prod_cu2 = dot_prod_cu * dot_prod_cu;
                    dot_prod_uu = Velzero * Velzero;

                    Feq[i] = LATTICE_WEIGHTS[i] * (RhoInt) *
                             (1.0 + dot_prod_cu * C_S_POW2_INV +
                              dot_prod_cu2 * C_S_POW4_INV / 2.0 -
                              dot_prod_uu * C_S_POW2_INV / 2.0);

                    /* Probability distribution function can not be less
                     * than 0 */
                    if (Feq[i] < 0)
                        ERROR("Probability distribution function can not "
                              "be negative.");
                    stream_field[Q_LBM *
                                     (x + y * ystep + z * zstep * zstep) +
                                 i] = Feq[i];
                    collide_field[Q_LBM *
                                      (x + y * ystep + z * zstep * zstep) +
                                  i] = Feq[i];
                }
                //}
            }
        }
    }
}

void InitialiseGrid(int xlength, int& xstep, int& ystep, int& zstep, int& xstart, int& ystart, int& zstart, int& xend, int& yend, int& zend) {
#ifdef D3Q19
    xstep = xlength + 2;
    ystep = xlength + 2;
    zstep = xlength + 2;
#endif // D3Q19
#ifdef D2Q9
    xstep = xlength + 2;
    ystep = xlength + 2;
    zstep = 1;
    // x = 0  is  dummy wall
    xstart = 1; //fluid region //wall is 0 + 1/2
    // x = xstep -1 is dummy wall
    xend = xstep - 2; // fluid region //wall is xstep - 1 - 1/2
    // y = 0 is dummy wall
    ystart = 1; // fluid region // wall is 0 + 1/2
    // y = ystep -1 // is dummy wall
    yend = ystep - 2; //fluid region // wall is ystep -1 -1/2 
    // z = 0
    zstart = 0;
    // z = zstep - 1
    zend = 0;
#endif // D2Q9
}
