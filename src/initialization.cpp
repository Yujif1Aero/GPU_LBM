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

void ReadParameters(float* LengthRef, float* tau, float* nu, float* VelInt, float* RhoInt, float* PresInt, float* velocity_bc, float* pressure_bc, int* timesteps, int* timesteps_per_plotting, char* filegrid, int argc, char* argv[], int* gpu_enabled) {
    /* NOTE ALL VARIABLES are NON-DIMENTIONAL
    L = L* / dx*, t = t* /dt*, rho = rho* /rho0*, Re = u*L* / nu* = uL/nu
    * means dimensional */

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
    READ_FLOAT(argv[1], re);

    velocity_bcx = &velocity_bc[0];
    velocity_bcy = &velocity_bc[1];
    velocity_bcz = &velocity_bc[2];

    VelIntx = &VelInt[0];
    VelInty = &VelInt[1];
    VelIntz = &VelInt[2];

    *RhoInt = 1.0;
    *PresInt = C_S_POW2 * *RhoInt;
    READ_FLOAT(argv[1], *VelIntx);
    READ_FLOAT(argv[1], *VelInty);
    READ_FLOAT(argv[1], *VelIntz);
    READ_FLOAT(argv[1], *velocity_bcx);
    READ_FLOAT(argv[1], *velocity_bcy);
    READ_FLOAT(argv[1], *velocity_bcz);

    *pressure_bc = C_S_POW2 * *RhoInt;

    READ_FLOAT(argv[1], *LengthRef);
    READ_INT(argv[1], *timesteps);
    READ_INT(argv[1], *timesteps_per_plotting);
    abs_vel_bc = sqrt(*velocity_bcx * *velocity_bcx +
                      *velocity_bcy * *velocity_bcy +
                      *velocity_bcz * *velocity_bcz);
    *nu = abs_vel_bc * (*LengthRef) / re;
    *tau = C_S_POW2_INV * *nu + 0.5;

    READ_STRING(argv[1], filegrid);
}

void InitialiseFields(float* collide_field, float* stream_field, int xmax, int ymax, int zmax, int gpu_enabled, float Feq[Q_LBM], float RhoInt, float* VelInt, vector<int>& bcd) {

    float dot_prod_cu, dot_prod_cu2, dot_prod_uu, intvel = 0.0;
    for (int z = 0; z < zmax; z++) {
        for (int y = 0; y < ymax; y++) {
            for (int x = 0; x < xmax; x++) {
                if (bcd[x + y * xmax + z * xmax * ymax] == FLUID) {

                    for (int i = 0; i < Q_LBM; i++) {

                        /* Initializing condition */
                        dot_prod_cu = LATTICE_VELOCITIES[i][0] * VelInt[0] +
                                      LATTICE_VELOCITIES[i][1] * VelInt[1] +
                                      LATTICE_VELOCITIES[i][2] * VelInt[2];
                        dot_prod_cu2 = dot_prod_cu * dot_prod_cu;
                        dot_prod_uu = VelInt[0] * VelInt[0] + VelInt[1] * VelInt[1] + VelInt[2] * VelInt[2];

                        Feq[i] = LATTICE_WEIGHTS[i] * (RhoInt) *
                                 (1.0 + dot_prod_cu * C_S_POW2_INV / pow(C, 2.0) +
                                  dot_prod_cu2 * C_S_POW4_INV / (2.0 * pow(C, 4.0)) -
                                  dot_prod_uu * C_S_POW2_INV / (2.0 * pow(C, 2.0)));

                        /* Probability distribution function can not be less
                         * than 0 */
                        if (Feq[i] < 0)
                            ERROR("Probability distribution function can not "
                                  "be negative.");
                        stream_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i] = Feq[i];
                        collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i] = Feq[i];
                    }
                } else {
                    for (int i = 0; i < Q_LBM; i++) {

                        /* Initializing condition */
                        dot_prod_cu = LATTICE_VELOCITIES[i][0] * intvel +
                                      LATTICE_VELOCITIES[i][1] * intvel +
                                      LATTICE_VELOCITIES[i][2] * intvel;
                        dot_prod_cu2 = dot_prod_cu * dot_prod_cu;
                        dot_prod_uu = intvel * intvel + intvel * intvel + intvel * intvel;
                        Feq[i] = LATTICE_WEIGHTS[i] * (RhoInt) *
                                 (1.0 + dot_prod_cu * C_S_POW2_INV +
                                  dot_prod_cu2 * C_S_POW4_INV / (2.0 * pow(C, 4.0)) -
                                  dot_prod_uu * C_S_POW2_INV / (2.0 * pow(C, 2.0)));

                        /* Probability distribution function can not be less
                         * than 0 */
                        if (Feq[i] < 0)
                            ERROR("Probability distribution function can not "
                                  "be negative.");
                        stream_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i] = Feq[i];
                        collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i] = Feq[i];
                    }
                }
            }
        }
    }
}

void InitialiseGrid(char* dirmesh, int& xmax, int& ymax, int& zmax, int& xstart, int& ystart, int& zstart, int& xend, int& yend, int& zend, vector<double>& xd, vector<double>& yd, vector<double>& zd, vector<int>& bcd) {
#ifdef D2Q9
    ifstream mesh;
    string line;
    mesh.open(dirmesh, ios::in);
    if (!mesh.is_open()) {
        cerr << "Failed to open the file!" << endl;
    }
    int idomein = 0;
    int i = 0, index = 0;

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
            xstart = 1;
            ystart = 1;
            zstart = 0;
            xend = xmax - 2;
            yend = ymax - 2;
            zend = 0;
            printf("xstart, ystart, zstart, xend, yend, zend are %d, %d, %d, %d, %d, %d\n", xstart, ystart, zstart, xend, yend, zend);
            idomein = xmax * ymax * zmax;
            printf("xmax, ymax, zmax, xmax*ymax*zmax are %d, %d, %d, %d\n", xmax, ymax, zmax, idomein);
            /* resize the arrays to the correct size */
            xd.resize(idomein);
            yd.resize(idomein);
            zd.resize(idomein);
            bcd.resize(idomein);

        } else {
            /* read the x, y, z values from the line*/
            double xi = 0.0, yi = 0.0, zi = 0.0;
            int bcdi = 0;
            sscanf(line.c_str(), "%lf  %lf  %lf  %d", &xi, &yi, &zi, &bcdi);
            /* store the vb*/
            index = i - 3;
            xd[index] = xi;
            yd[index] = yi;
            zd[index] = zi;
            bcd[index] = bcdi;
        }
        i++;
    }

#endif // D2Q9
}