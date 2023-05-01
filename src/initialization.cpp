#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
// #include <unistd.h>
using namespace std;

#include "initialization.hpp"
#include "lbm_model.hpp"
#include "utils.hpp"

#define WORKAROUND 1

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

void ReadParameters(int* xlength, float* tau, float& VelInt, float* velocity_wall, int* timesteps, int* timesteps_per_plotting, int argc, char* argv[], int* gpu_enabled) {
    /*  cavity flow */
    float *velocity_wall_1, *velocity_wall_2, *velocity_wall_3, re, nu;
    /* printing out help message */
    if (argc < 3)
        PrintHelpMessage();
    if (!strcmp(argv[1], "-help") || !strcmp(argv[2], "-help"))
        PrintHelpMessage();

        /* checking parameters */
#if !defined(WORKAROUND)
    if (access(argv[1], R_OK) != 0)
        ERROR("Provided configuration file path either doesn't exist or can "
              "not be read.");
#endif
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

void ReadParameters(float* LengthRef, float* tau, float* nu, float* VelInt, float* RhoInt, float* PresInt, float* velocity_bc, float* pressure_bc, int* timesteps, int* timesteps_per_plotting, int argc, char* argv[], int* gpu_enabled) {
    /* NOTE ALL VARIABLES are NON-DIMENTIONAL
    L = L* / dx*, t = t* /dt*, rho = rho* /rho0*, Re = u*L* / nu* = uL/nu
    * means dimensional */

    /*  open boundary */
    float re, abs_vel_bc;
    float *velocity_bcx, *velocity_bcy, *velocity_bcz;
    float *VelIntx, *VelInty, *VelIntz;
    /* printing out help message */
    if (argc < 3)
        PrintHelpMessage();
    if (!strcmp(argv[1], "-help") || !strcmp(argv[2], "-help"))
        PrintHelpMessage();

        /* checking parameters */
#if !defined(WORKAROUND)
    if (access(argv[1], R_OK) != 0)
        ERROR("Provided configuration file path either doesn't exist or can "
              "not be read.");
#endif
    *gpu_enabled = 0;

    /* reading parameters */
    // READ_FLOAT(argv[1], *tau);
    READ_FLOAT(argv[1], re);

    velocity_bcx = &velocity_bc[0];
    velocity_bcy = &velocity_bc[1];
    velocity_bcz = &velocity_bc[2];

    VelIntx = &VelInt[0];
    VelInty = &VelInt[1];
    VelIntz = &VelInt[2];

    // READ_FLOAT(argv[1], *RhoInt);
    // READ_FLOAT(argv[1], *PresInt);
    *RhoInt = 1.0;
    *PresInt = C_S_POW2 * *RhoInt;
    READ_FLOAT(argv[1], *VelIntx);
    READ_FLOAT(argv[1], *VelInty);
    READ_FLOAT(argv[1], *VelIntz);
    READ_FLOAT(argv[1], *velocity_bcx);
    READ_FLOAT(argv[1], *velocity_bcy);
    READ_FLOAT(argv[1], *velocity_bcz);
    // READ_FLOAT(argv[1], *pressure_bc);
    *pressure_bc = C_S_POW2 * *RhoInt;

    READ_FLOAT(argv[1], *LengthRef);
    READ_INT(argv[1], *timesteps);
    READ_INT(argv[1], *timesteps_per_plotting);
    abs_vel_bc = sqrt(*velocity_bcx * *velocity_bcx +
                      *velocity_bcy * *velocity_bcy +
                      *velocity_bcz * *velocity_bcz);
    *nu = abs_vel_bc * (*LengthRef) / re;
    *tau = C_S_POW2_INV * *nu + 0.5;
}

// void InitialiseFields(float* collide_field, float* stream_field, int xmax, int ymax, int zmax, int gpu_enabled, float Feq[Q_LBM], float RhoInt, float wall_velocity[D_LBM]) {
//     /* only cavity*/
//     int x, y, z, i;
//     float dot_prod_cu, dot_prod_cu2, dot_prod_uu;
//     float VelInt = 0.0;
//     float Velzero = 0.0;
//     /* NOTE: We use y=ymax-1 as the moving wall */
//     /* Initialization including dummy wall*/
//     for (z = 0; z < zmax; z++) {
//         for (y = 0; y < ymax; y++) {
//             for (x = 0; x < xmax; x++) {

//                 for (i = 0; i < Q_LBM; i++) {

//                     /* Initializing condition */
//                     dot_prod_cu = LATTICE_VELOCITIES[i][0] * Velzero +
//                                   LATTICE_VELOCITIES[i][1] * Velzero +
//                                   LATTICE_VELOCITIES[i][2] * Velzero;
//                     dot_prod_cu2 = dot_prod_cu * dot_prod_cu;
//                     dot_prod_uu = Velzero * Velzero;

//                     Feq[i] = LATTICE_WEIGHTS[i] * (RhoInt) *
//                              (1.0 + dot_prod_cu * C_S_POW2_INV +
//                               dot_prod_cu2 * C_S_POW4_INV / 2.0 -
//                               dot_prod_uu * C_S_POW2_INV / 2.0);

//                     /* Probability distribution function can not be less
//                      * than 0 */
//                     if (Feq[i] < 0)
//                         ERROR("Probability distribution function can not "
//                               "be negative.");
//                     stream_field[Q_LBM * (x + y * ymax + z * zmax * zmax) + i] = Feq[i];
//                     collide_field[Q_LBM * (x + y * ymax + z * zmax * zmax) + i] = Feq[i];
//                 }
//             }
//         }
//     }
// }

void InitialiseFields(float* collide_field, float* stream_field, int xmax, int ymax, int zmax, int gpu_enabled, float Feq[Q_LBM], float RhoInt, float* VelInt) {
    /* open boundary */
    float dot_prod_cu, dot_prod_cu2, dot_prod_uu;
    for (int z = 0; z < zmax; z++) {
        for (int y = 0; y < ymax; y++) {
            for (int x = 0; x < xmax; x++) {
                for (int i = 0; i < Q_LBM; i++) {

                    /* Initializing condition */
                    dot_prod_cu = LATTICE_VELOCITIES[i][0] * VelInt[0] +
                                  LATTICE_VELOCITIES[i][1] * VelInt[1] +
                                  LATTICE_VELOCITIES[i][2] * VelInt[2];
                    dot_prod_cu2 = dot_prod_cu * dot_prod_cu;
                    dot_prod_uu = VelInt[0] * VelInt[0] + VelInt[1] * VelInt[2] * VelInt[2];

                    Feq[i] = LATTICE_WEIGHTS[i] * (RhoInt) *
                             (1.0 + dot_prod_cu * C_S_POW2_INV +
                              dot_prod_cu2 * C_S_POW4_INV / 2.0 -
                              dot_prod_uu * C_S_POW2_INV / 2.0);

                    /* Probability distribution function can not be less
                     * than 0 */
                    if (Feq[i] < 0)
                        ERROR("Probability distribution function can not "
                              "be negative.");
                    stream_field[Q_LBM * (x + y * ymax + z * zmax * zmax) + i] = Feq[i];
                    collide_field[Q_LBM * (x + y * ymax + z * zmax * zmax) + i] = Feq[i];
                }
            }
        }
    }
}

void InitialiseGrid(int xlength, int& xmax, int& ymax, int& zmax, int& xstart, int& ystart, int& zstart, int& xend, int& yend, int& zend) {
#ifdef D3Q19
    xmax = xlength + 2;
    ymax = xlength + 2;
    zmax = xlength + 2;
#endif // D3Q19
#ifdef D2Q9
    xmax = xlength + 2;
    ymax = xlength + 2;
    zmax = 1;
    // x = 0  is  dummy wall
    xstart = 1; // fluid region //wall is 0 + 1/2
    // x = xmax -1 is dummy wall
    xend = xmax - 2; // fluid region //wall is xmax - 1 - 1/2
    // y = 0 is dummy wall
    ystart = 1; // fluid region // wall is 0 + 1/2
    // y = ymax -1 // is dummy wall
    yend = ymax - 2; // fluid region // wall is ymax -1 -1/2
    // z = 0
    zstart = 0;
    // z = zmax - 1
    zend = 0;
#endif // D2Q9
}

void InitialiseGrid(int &xmax, int &ymax, int &zmax, int& xstart, int& ystart, int& zstart, int& xend, int& yend, int& zend, vector<double>& xd, vector<double>& yd, vector<double>& zd) {
#ifdef D2Q9
    ifstream mesh("./mesh/grid.dat");
    string line;

    int idomein = 0;
    int i = 0, index = 0;
    //int jmax = 0, kmax = 0, lmax = 0;

    while (getline(mesh, line)) {
        if (i == 0) {
            printf("skip the TITLE line \n");
            cout << line.substr() << endl;
        } else if (i == 1) {
            printf("skip the VARIABLE line \n");
            cout << line.substr() << endl;
        } else if (line.substr(0, 4) == "ZONE") {
            /* read the grid dimensions from the ZONE line */
            
            sscanf(line.c_str(), "ZONE T = \"%*[^\"]\", I=%d, J=%d, K=%d, F=%*s", &xmax, &ymax, &zmax);
            cout << line.substr() << endl;
            xstart = 0; 
            ystart = 0;
            zstart = 0;
            xend = xmax - 1;
            yend = ymax - 1;
            zend = zmax;

            idomein = xmax * ymax * zmax;
            printf("jmax, kmax, lmax, jmax*kmax*lmax are %d, %d, %d, %d\n", xmax, ymax, zmax, idomein);
            /* resize the arrays to the correct size */
            xd.resize(idomein);
            yd.resize(idomein);
            zd.resize(idomein);

        } else {
            /* read the x, y, z values from the line*/
            double xi = 0.0, yi = 0.0, zi = 0.0;
            sscanf(line.c_str(), "%lf, %lf, %lf", &xi, &yi, &zi);
            /* store the vb*/
            index = i - 3;
            xd[index] = xi;
            yd[index] = yi;
            zd[index] = zi;
        }
        i++;
    }

#endif // D2Q9
}