#include <cmath>
#include <iostream>
#include <vector>
using namespace std;
#include "boundary.hpp"
#include "cell_computation.hpp"
#include "lbm_model.hpp"

#define SWICH_OUTLETBC
/**
 * Finds an inverse probability distribution for the provided lattice index.
 */
int inv(int i) {
    return (Q_LBM - 1) - i;
}

void TreatBoundary(float* collide_field, float* wall_velocity, int xstart, int ystart, int zstart, int xend, int yend, int zend, int xmax, int ymax, int zmax) {
    int x, nx, y, ny, z, nz, i;
    // int invC = int(1 / C);
    float dot_prod, density;
    int imxy = 7, imx0 = 3, imxmy = 0;

    // halfway bounce back scheme
    for (z = zstart; z <= zend; z++) {
        for (y = ystart; y <= yend; y++) {
            // non-slip boundary left wall
            // nx = 0; //dummy wall
            nx = xstart; // fluid
            i = imxy;
            /* Assign the boudary cell value */
            collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * ((nx - 1) + (y + 1) * xmax + z * xmax * ymax) + i];

            i = imx0;
            /* Assign the boudary cell value */
            collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * ((nx - 1) + y * xmax + z * xmax * ymax) + i];

            i = imxmy;
            /* Assign the boudary cell value */
            collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * ((nx - 1) + (y - 1) * xmax + z * xmax * ymax) + i];
        }
        for (y = ystart; y <= yend; y++) {
            // non-slip boundary right  wall
            // nx = xmax -1; //dummy wall
            nx = xend; // fluid
            i = inv(imxy);
            /* Assign the boudary cell value */
            collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * ((nx + 1) + (y - 1) * xmax + z * xmax * ymax) + i];

            i = inv(imx0);
            /* Assign the boudary cell value */
            collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * ((nx + 1) + y * xmax + z * xmax * ymax) + i];

            i = inv(imxmy);
            /* Assign the boudary cell value */
            collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * ((nx + 1) + (y + 1) * xmax + z * xmax * ymax) + i];
        }
    }

    for (z = zstart; z <= zend; z++) {
        for (x = xstart; x <= xend; x++) {
            // halfway bounce back scheme
            // no-slip boundary bottom wall
            // ny = 0; // dummy wall
            ny = ystart; // fluid
            for (i = 0; i < 3; i++) {

                if (i == 0)
                    collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * ((x - 1) + (ny - 1) * xmax + z * xmax * ymax) + i];
                if (i == 1)
                    collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * ((x + 1) + (ny - 1) * xmax + z * xmax * ymax) + i];
                if (i == 2)
                    collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * (x + (ny - 1) * xmax + z * xmax * ymax) + i];
            }
            // halfway bounce back scheme
            // moving-slip boundary top wall
            // ny = yend +1 ; // dummy wall
            ny = yend; // fluid
            for (i = 6; i < Q_LBM; i++) {
                dot_prod = LATTICE_VELOCITIES[i][0] * wall_velocity[0] + LATTICE_VELOCITIES[i][1] * wall_velocity[1] + LATTICE_VELOCITIES[i][2] * wall_velocity[2];
                /* Compute density in the neighbour cell */
                ComputeDensity(&collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax)], &density);
                /* Assign the boudary cell value */
                if (i == 6)
                    collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax) + inv(i)] =
                        collide_field[Q_LBM * (x + (ny + 1) * xmax + z * xmax * ymax) + i] - 6.0 * density * LATTICE_WEIGHTS[i] * dot_prod;
                if (i == 7)
                    collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax) + inv(i)] =
                        collide_field[Q_LBM * ((x - 1) + (ny + 1) * xmax + z * xmax * ymax) + i] - 6.0 * density * LATTICE_WEIGHTS[i] * dot_prod;
                if (i == 8)
                    collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax) + inv(i)] =
                        collide_field[Q_LBM * ((x + 1) + (ny + 1) * xmax + z * xmax * ymax) + i] - 6.0 * density * LATTICE_WEIGHTS[i] * dot_prod;
            }
        }
    }
}

void TreatBoundary(float* collide_field, float* velocity_bc, int xstart, int ystart, int zstart, int xend, int yend, int zend, int xmax, int ymax, int zmax, vector<int>& bcd, float RhoInt) {

    float dot_prod, density, velocity[3];
    float dot_prod_cu, dot_prod_cu2, dot_prod_uu, dot_prod_cu2_cs, dot_prod_uu2, dot_prod_uu2_cs;
    float Feq[Q_LBM];
    int nx, ny, nz;
    int nyperimax, nyperimin;
    for (int z = 0; z < zmax; z++) {
        for (int y = 0; y < ymax; y++) {
            for (int x = 0; x < xmax; x++) {
                if (bcd[x + y * xmax + z * xmax * ymax] == FLUID) {
                    for (int i = 0; i < Q_LBM; i++) {
                        nx = x + int(LATTICE_VELOCITIES[i][0]);
                        ny = y + int(LATTICE_VELOCITIES[i][1]);
                        nz = z + int(LATTICE_VELOCITIES[i][2]);

                        // cout << "nx, ny and nz are" << nx << "," << ny << "," << nz << endl;
                        if (bcd[nx + ny * xmax + nz * xmax * ymax] == NO_SLIP) {
                            /* bounce back ref Seta Book*/
                            collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i];
                            // collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i];
                        } else if (bcd[nx + ny * xmax + nz * xmax * ymax] == INLET) {
                            /* Dirichlet boundary condition ref book Seta and Springer LBM at PP 180*/
                            dot_prod = LATTICE_VELOCITIES[i][0] * velocity_bc[0] + LATTICE_VELOCITIES[i][1] * velocity_bc[1] + LATTICE_VELOCITIES[i][2] * velocity_bc[2];

                            ComputeDensity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density);
                            //density = RhoInt;

                            collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] - 2.0 * C_S_POW2_INV * density * LATTICE_WEIGHTS[i] * dot_prod;
                        }
                    }
                }
            }
        }
    }

    for (int z = 0; z < zmax; z++) {
        for (int y = 0; y < ymax; y++) {
            for (int x = 0; x < xmax; x++) {
                if (bcd[x + y * xmax + z * xmax * ymax] == OUTLET) {
                    if (x == 0) {
                        nx = 1;
                    } else if (x == xmax - 1) {
                        nx = xmax - 2;

                        for (int i = 0; i < Q_LBM; i++) {
// cout << "outlet" << bcd[nx + ny * xmax + nz * xmax * ymax] << endl;
#ifdef SWICH_OUTLETBC
                            /* Springer LBM at PP200*/
                            ComputeDensity(&collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax)], &density);
                            ComputeVelocity(&collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax)], &density, velocity);

                            dot_prod_cu = LATTICE_VELOCITIES[i][0] * velocity[0] +
                                          LATTICE_VELOCITIES[i][1] * velocity[1] +
                                          LATTICE_VELOCITIES[i][2] * velocity[2];
                            dot_prod_cu2 = dot_prod_cu * dot_prod_cu;

                            dot_prod_uu = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];

                            Feq[i] = LATTICE_WEIGHTS[i] * density * (1.0 + dot_prod_cu * C_S_POW2_INV / pow(C, 2.0) + dot_prod_cu2 * C_S_POW4_INV / (2.0 * pow(C, 4.0)) - dot_prod_uu * C_S_POW2_INV / (2.0 * pow(C, 2.0)));

                            collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i] = Feq[i];
                            // collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i] = Feq[i];
#else
                           // cout << "nx, x, y =" << nx << " " << x << " " << y << " " << endl;
                            collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax) + i];

#endif
                        }
                    }
                    if (y == 0) {
                        ny = 1;
                    } else if (y == ymax - 1) {
                        ny = ymax - 2;

                        for (int i = 0; i < Q_LBM; i++) {
// cout << "outlet" << bcd[nx + ny * xmax + nz * xmax * ymax] << endl;
#ifdef SWICH_OUTLETBC
                            /* Springer LBM at PP200*/
                            ComputeDensity(&collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax)], &density);
                            ComputeVelocity(&collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax)], &density, velocity);

                            dot_prod_cu = LATTICE_VELOCITIES[i][0] * velocity[0] +
                                          LATTICE_VELOCITIES[i][1] * velocity[1] +
                                          LATTICE_VELOCITIES[i][2] * velocity[2];
                            dot_prod_cu2 = dot_prod_cu * dot_prod_cu;

                            dot_prod_uu = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];

                            Feq[i] = LATTICE_WEIGHTS[i] * density * (1.0 + dot_prod_cu * C_S_POW2_INV / pow(C, 2.0) + dot_prod_cu2 * C_S_POW4_INV / (2.0 * pow(C, 4.0)) - dot_prod_uu * C_S_POW2_INV / (2.0 * pow(C, 2.0)));

                            collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i] = Feq[i];
#else
                            collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax) + i];

#endif
                        }
                    }
                }
            }
        }
    }

    for (int z = 0; z < zmax; z++) {
        for (int y = 0; y < ymax; y++) {
            for (int x = 0; x < xmax; x++) {
                if (bcd[x + y * xmax + z * xmax * ymax] == PERIODIC) {
                    for (int i = 0; i < Q_LBM; i++) {
                        if (y == ymax - 1) {
                            // cout << "periodic " << ny <<"x, y, z ="<< nx << " " << ny << " " << nz << endl;
                            nyperimin = 2;
                            collide_field[Q_LBM * (x + (y - 1) * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (x + nyperimin * xmax + z * xmax * ymax) + i];

                        } else if (y == 0) {
                            // cout << "periodic " << ny <<"x, y, z ="<< nx << " " << ny << " " << nz << endl;
                            nyperimax = ymax - 3;
                            collide_field[Q_LBM * (x + (y + 1) * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (x + nyperimax * xmax + z * xmax * ymax) + i];
                        }
                    }
                }
            }
        }
    }
}