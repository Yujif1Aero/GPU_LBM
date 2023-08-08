#include <cmath>
#include <iostream>
#include <vector>
using namespace std;
#include "boundary.hpp"
#include "cell_computation.hpp"
#include "lbm_model.hpp"

// #define SWICH_OUTLETBC
//  #define SWICH_OUTLETBC1
// #define SWICH_OUTLETBC_PRESSUREOUT_20230803_1
#define SWICH_OUTLETBC_PRESSUREOUT_20230803_2
#define SWICH_OUTLETBC_NEUMANN_20230803
/**
 * Finds an inverse probability distribution for the provided lattice index.
 */
int inv(int i) {
    return (Q_LBM - 1) - i;
}

void TreatBoundary(float* collide_field, float* velocity_bc, int xstart, int ystart, int zstart, int xend, int yend, int zend, int xmax, int ymax, int zmax, vector<int>& bcd, float RhoInt) {

    float density, pressure, density_b, density_b1, velocity[3], velocity_b[3], velocity_b1[3];
    float dot_prod, dot_prod_cu, dot_prod_cu2, dot_prod_uu, dot_prod_cu2_cs, dot_prod_uu2, dot_prod_uu2_cs;
    float RHS[Q_LBM];
    float work;
    int nx, ny, nz;
    int nx_b1, ny_b1, nz_b1 = 0;
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

                            /* Dirichlet boundary condition ref book Seta and Springer LBM at PP 180 cs is at PP 92*/
                            dot_prod = LATTICE_VELOCITIES[i][0] * velocity_bc[0] + LATTICE_VELOCITIES[i][1] * velocity_bc[1] + LATTICE_VELOCITIES[i][2] * velocity_bc[2];

                            ComputeDensity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density);
                            // density = RhoInt;

                            collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] - 2.0 * C_S_POW2_INV * density * LATTICE_WEIGHTS[i] * dot_prod;
#ifdef SWICH_OUTLETBC_PRESSUREOUT_20230803_2
                        } else if (bcd[nx + ny * xmax + nz * xmax * ymax] == OUTLET) {

                            nx_b1 = x - int(LATTICE_VELOCITIES[i][0]);
                            ny_b1 = y - int(LATTICE_VELOCITIES[i][1]);
                            nz_b1 = z - int(LATTICE_VELOCITIES[i][2]);
                            // if (i == 8) {
                            //     cout << "nx_b1, ny_b1, nz_b1 = " << nx_b1 << " " << ny_b1 << " " << nz_b1 << " " << endl;
                            // }
                            if ((nx == xmax - 1 && ny != ymax - 1) || (nx == xmax - 1 && ny != 0) || (nx == 0 && ny != ymax - 1) || (nx == 0 && ny != 0)) {
                                // if (nx == xmax - 1 || nx == 0) {
                                ComputeDensity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b);
                                ComputeVelocity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b, velocity_b);
                                ComputeDensity(&collide_field[Q_LBM * (nx_b1 + y * xmax + z * xmax * ymax)], &density_b1);
                                ComputeVelocity(&collide_field[Q_LBM * (nx_b1 + y * xmax + z * xmax * ymax)], &density_b1, velocity_b1);
                                for (int j = 0; j < 3; j++) {
                                    velocity[j] = velocity_b[j] + 0.5 * (velocity_b[j] - velocity_b1[j]);
                                }
                                density = RhoInt;
                                // work = velocity_bc[0] * velocity_bc[0] + velocity_bc[1] * velocity_bc[1] + velocity_bc[2] * velocity_bc[2];
                                // pressure = RhoInt * C_S_POW2 + 0.5 * RhoInt * work;
                                // density = pressure * C_S_POW2_INV;
                            } else if ((ny == ymax - 1 && nx != xmax - 1) || (ny == ymax - 1 && nx != 0) || (ny == 0 && nx != xmax - 1) || (ny == 0 && nx != 0)) {
                                //  } else if (ny == ymax - 1 || ny == 0) {
                                ComputeDensity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b);
                                ComputeVelocity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b, velocity_b);
                                ComputeDensity(&collide_field[Q_LBM * (x + ny_b1 * xmax + z * xmax * ymax)], &density_b1);
                                ComputeVelocity(&collide_field[Q_LBM * (x + ny_b1 * xmax + z * xmax * ymax)], &density_b1, velocity_b1);
                                for (int j = 0; j < 3; j++) {
                                    velocity[j] = velocity_b[j] + 0.5 * (velocity_b[j] - velocity_b1[j]);
                                }
                                // work = velocity_bc[0] * velocity_bc[0] + velocity_bc[1] * velocity_bc[1] + velocity_bc[2] * velocity_bc[2];
                                // pressure = RhoInt * C_S_POW2 + 0.5 * RhoInt * work;
                                // density = pressure * C_S_POW2_INV;
                                density = RhoInt;
                                // cout << density << endl;
                            } else if ((nx == xmax - 1 && ny == ymax - 1) || (nx == 0 && ny == 0) || (nx == xmax - 1 && ny == 0) || (nx == 0 && ny == ymax - 1)) {
                                // conner
                                ComputeDensity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b);
                                ComputeVelocity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b, velocity_b);
                                ComputeDensity(&collide_field[Q_LBM * (nx_b1 + ny_b1 * xmax + z * xmax * ymax)], &density_b1);
                                ComputeVelocity(&collide_field[Q_LBM * (nx_b1 + ny_b1 * xmax + z * xmax * ymax)], &density_b1, velocity_b1);
                                for (int j = 0; j < 3; j++) {
                                    velocity[j] = velocity_b[j] + 0.5 * (velocity_b[j] - velocity_b1[j]);
                                    // velocity[j] = velocity_b[j];
                                }
                                // density = 0.5 * (density_b + density_b1);
                                //  density = density_b;
                                density = RhoInt;
                            }
                            //     for (int j = 0; j < 3; j++) {
                            //         if ((nx == xmax - 1 && ny == ymax - 1) || (nx == 0 && ny == 0) || (nx == xmax - 1 && ny == 0) || (nx == 0 && ny == ymax - 1)) {
                            //             velocity[j] = velocity_b[j];
                            //         } else {
                            //             velocity[j] = velocity_b[j] + 0.5 * (velocity_b[j] - velocity_b1[j]);
                            //         }
                            //     }
                            //   density = density_b;

                            dot_prod_cu = LATTICE_VELOCITIES[i][0] * velocity[0] +
                                          LATTICE_VELOCITIES[i][1] * velocity[1] +
                                          LATTICE_VELOCITIES[i][2] * velocity[2];
                            dot_prod_cu2 = dot_prod_cu * dot_prod_cu;

                            dot_prod_uu = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];

                            RHS[i] = -collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] + 2. * LATTICE_WEIGHTS[i] * density * (1. + 0.5 * dot_prod_cu2 * C_S_POW4_INV - 0.5 * dot_prod_uu * C_S_POW2_INV);

                            // cout << "dencity= " << density << endl;
                            // dot_prod = LATTICE_VELOCITIES[i][0] * velocity[0] + LATTICE_VELOCITIES[i][1] * velocity[1] + LATTICE_VELOCITIES[i][2] * velocity[2];

                            // collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] - 2.0 * C_S_POW2_INV * density * LATTICE_WEIGHTS[i] * dot_prod;
                            collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = RHS[i];
#endif
#ifdef SWICH_OUTLETBC_PRESSUREOUT_20230803_1
                        } else if (bcd[nx + ny * xmax + nz * xmax * ymax] == OUTLET) {
                            /* Dirichlet boundary condition ref book Springer LBM at PP 200*/
                            if (nx == 0) {
                                nx_b1 = x + 1;
                                ComputeDensity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b);
                                ComputeVelocity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b, velocity_b);
                                ComputeDensity(&collide_field[Q_LBM * (nx_b1 + y * xmax + z * xmax * ymax)], &density_b1);
                                ComputeVelocity(&collide_field[Q_LBM * (nx_b1 + y * xmax + z * xmax * ymax)], &density_b1, velocity_b1);
                                for (int j = 0; j < 3; j++) {
                                    velocity[j] = velocity_b[j] + 0.5 * (velocity_b[j] - velocity_b1[j]);
                                }
                                density = RhoInt;

                                dot_prod_cu = LATTICE_VELOCITIES[i][0] * velocity[0] +
                                              LATTICE_VELOCITIES[i][1] * velocity[1] +
                                              LATTICE_VELOCITIES[i][2] * velocity[2];
                                dot_prod_cu2 = dot_prod_cu * dot_prod_cu;

                                dot_prod_uu = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];

                                // dot_prod_uu2 *= dot_prod_uu;

                                RHS[i] = -collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] + 2. * LATTICE_WEIGHTS[i] * density * (1. + 0.5 * dot_prod_cu2 * C_S_POW4_INV - 0.5 * dot_prod_uu * C_S_POW2_INV);

                                // cout << "dencity= " << density << endl;
                                // dot_prod = LATTICE_VELOCITIES[i][0] * velocity[0] + LATTICE_VELOCITIES[i][1] * velocity[1] + LATTICE_VELOCITIES[i][2] * velocity[2];

                                // collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] - 2.0 * C_S_POW2_INV * density * LATTICE_WEIGHTS[i] * dot_prod;
                                collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = RHS[i];
                            } else if (nx == xmax - 1) {
                                nx_b1 = x - 1;
                                ComputeDensity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b);
                                ComputeVelocity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b, velocity_b);
                                ComputeDensity(&collide_field[Q_LBM * (nx_b1 + y * xmax + z * xmax * ymax)], &density_b1);
                                ComputeVelocity(&collide_field[Q_LBM * (nx_b1 + y * xmax + z * xmax * ymax)], &density_b1, velocity_b1);
                                for (int j = 0; j < 3; j++) {
                                    velocity[j] = velocity_b[j] + 0.5 * (velocity_b[j] - velocity_b1[j]);
                                }
                                // density = RhoInt;
                                work = velocity_bc[0] * velocity_bc[0] + velocity_bc[1] * velocity_bc[1] + velocity_bc[2] * velocity_bc[2];
                                pressure = RhoInt * C_S_POW2 + 0.5 * RhoInt * work;
                                density = pressure * C_S_POW2_INV;

                                dot_prod_cu = LATTICE_VELOCITIES[i][0] * velocity[0] +
                                              LATTICE_VELOCITIES[i][1] * velocity[1] +
                                              LATTICE_VELOCITIES[i][2] * velocity[2];
                                dot_prod_cu2 = dot_prod_cu * dot_prod_cu;

                                dot_prod_uu = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];

                                // dot_prod_uu2 *= dot_prod_uu;

                                RHS[i] = -collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] + 2. * LATTICE_WEIGHTS[i] * density * (1. + 0.5 * dot_prod_cu2 * C_S_POW4_INV - 0.5 * dot_prod_uu * C_S_POW2_INV);

                                // cout << "dencity= " << density << endl;
                                // dot_prod = LATTICE_VELOCITIES[i][0] * velocity[0] + LATTICE_VELOCITIES[i][1] * velocity[1] + LATTICE_VELOCITIES[i][2] * velocity[2];

                                // collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] - 2.0 * C_S_POW2_INV * density * LATTICE_WEIGHTS[i] * dot_prod;
                                collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = RHS[i];

                            } else if (ny == 0) {
                                ny_b1 = y + 1;
                                ComputeDensity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b);
                                ComputeVelocity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b, velocity_b);
                                ComputeDensity(&collide_field[Q_LBM * (x + ny_b1 * xmax + z * xmax * ymax)], &density_b1);
                                ComputeVelocity(&collide_field[Q_LBM * (x + ny_b1 * xmax + z * xmax * ymax)], &density_b1, velocity_b1);
                                for (int j = 0; j < 3; j++) {
                                    velocity[j] = velocity_b[j] + 0.5 * (velocity_b[j] - velocity_b1[j]);
                                }
                                // density = density_b;

                                dot_prod_cu = LATTICE_VELOCITIES[i][0] * velocity[0] +
                                              LATTICE_VELOCITIES[i][1] * velocity[1] +
                                              LATTICE_VELOCITIES[i][2] * velocity[2];
                                dot_prod_cu2 = dot_prod_cu * dot_prod_cu;

                                dot_prod_uu = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];

                                // dot_prod_uu2 *= dot_prod_uu;
                                work = velocity_bc[0] * velocity_bc[0] + velocity_bc[1] * velocity_bc[1] + velocity_bc[2] * velocity_bc[2];
                                pressure = RhoInt * C_S_POW2 + 0.5 * RhoInt * work;
                                density = pressure * C_S_POW2_INV;

                                RHS[i] = -collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] + 2. * LATTICE_WEIGHTS[i] * density * (1. + 0.5 * dot_prod_cu2 * C_S_POW4_INV - 0.5 * dot_prod_uu * C_S_POW2_INV);

                                // cout << "dencity= " << density << endl;
                                // dot_prod = LATTICE_VELOCITIES[i][0] * velocity[0] + LATTICE_VELOCITIES[i][1] * velocity[1] + LATTICE_VELOCITIES[i][2] * velocity[2];

                                // collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] - 2.0 * C_S_POW2_INV * density * LATTICE_WEIGHTS[i] * dot_prod;
                                collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = RHS[i];

                            } else if (ny == ymax - 1) {
                                ny_b1 = y - 1;
                                // cout << "y-1 = " << ny_b1 << endl;
                                ComputeDensity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b);
                                ComputeVelocity(&collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)], &density_b, velocity_b);
                                ComputeDensity(&collide_field[Q_LBM * (x + ny_b1 * xmax + z * xmax * ymax)], &density_b1);
                                ComputeVelocity(&collide_field[Q_LBM * (x + ny_b1 * xmax + z * xmax * ymax)], &density_b1, velocity_b1);
                                for (int j = 0; j < 3; j++) {
                                    velocity[j] = velocity_b[j] + 0.5 * (velocity_b[j] - velocity_b1[j]);
                                }
                                work = velocity_bc[0] * velocity_bc[0] + velocity_bc[1] * velocity_bc[1] + velocity_bc[2] * velocity_bc[2];
                                pressure = RhoInt * C_S_POW2 + 0.5 * RhoInt * work;
                                density = pressure * C_S_POW2_INV;

                                dot_prod_cu = LATTICE_VELOCITIES[i][0] * velocity[0] +
                                              LATTICE_VELOCITIES[i][1] * velocity[1] +
                                              LATTICE_VELOCITIES[i][2] * velocity[2];
                                dot_prod_cu2 = dot_prod_cu * dot_prod_cu;

                                dot_prod_uu = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];

                                // dot_prod_uu2 *= dot_prod_uu;

                                RHS[i] = -collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] + 2. * LATTICE_WEIGHTS[i] * density * (1. + 0.5 * dot_prod_cu2 * C_S_POW4_INV - 0.5 * dot_prod_uu * C_S_POW2_INV);

                                // cout << "dencity= " << density << endl;
                                // dot_prod = LATTICE_VELOCITIES[i][0] * velocity[0] + LATTICE_VELOCITIES[i][1] * velocity[1] + LATTICE_VELOCITIES[i][2] * velocity[2];

                                // collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = collide_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] - 2.0 * C_S_POW2_INV * density * LATTICE_WEIGHTS[i] * dot_prod;
                                collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + inv(i)] = RHS[i];
                            }

#endif
                        }
                    }
                }
#ifdef SWICH_OUTLETBC_NEUMANN_20230803
                if (bcd[x + y * xmax + z * xmax * ymax] == NEUMANN) {
                    //cout << "NEUMANN OUTLET  x, y, z =" << x << " " << y << " " << z << " " << endl;
                    if (x == 0) {
                        nx = 1;
                        nx_b1 = 2;

                        for (int i = 0; i < Q_LBM; i++) {
                            collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (nx_b1 + y * xmax + z * xmax * ymax) + i];
                        }
                    } else if (x == xmax - 1) {
                        nx = xmax - 2;
                        nx_b1 = xmax - 3;

                        for (int i = 0; i < Q_LBM; i++) {
                            collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (nx_b1 + y * xmax + z * xmax * ymax) + i];
                        }
                    } else if (y == 0) {
                        ny = 1;
                        ny_b1 = 2;

                        for (int i = 0; i < Q_LBM; i++) {
                            collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (x + ny_b1 * xmax + z * xmax * ymax) + i];
                        }
                    } else if (y == ymax - 1) {
                        ny = ymax - 2;
                        ny_b1 = ymax - 3;

                        for (int i = 0; i < Q_LBM; i++) {
                            collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (x + ny_b1 * xmax + z * xmax * ymax) + i];
                        }
                    }
                }

#endif
            }
        }
    }

#ifdef SWICH_OUTLETBC
    for (int z = 0; z < zmax; z++) {
        for (int y = 0; y < ymax; y++) {
            for (int x = 0; x < xmax; x++) {
                if (bcd[x + y * xmax + z * xmax * ymax] == OUTLET) {
                    if (x == 0) {
                        nx = 1;
                        nx_b1 = 2;
                        for (int i = 0; i < Q_LBM; i++) {
#ifdef SWICH_OUTLETBC1
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
                            // cout << "OUTLET x nx, nx_b1, x, y =" << x << " " << nx << " " << nx_b1 << " " << x << " " << y << " " << endl;
                            //  int nx_b1 = x - int(LATTICE_VELOCITIES[i][0]);
                            //  int ny_b1 = y - int(LATTICE_VELOCITIES[i][1]);
                            //  int nz_b1 = z - int(LATTICE_VELOCITIES[i][2]);
                            collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (nx_b1 + y * xmax + z * xmax * ymax) + i];

#endif
                        }
                    } else if (x == xmax - 1) {
                        nx = xmax - 2;
                        nx_b1 = xmax - 3;
                        for (int i = 0; i < Q_LBM; i++) {

#ifdef SWICH_OUTLETBC1
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
                            // cout << "OUTLET x nx, nx_b1, x, y =" << x << " " << nx << " " << nx_b1 << " " << x << " " << y << " " << endl;
                            //  int nx_b1 = x - int(LATTICE_VELOCITIES[i][0]);
                            //  int ny_b1 = y - int(LATTICE_VELOCITIES[i][1]);
                            //  int nz_b1 = z - int(LATTICE_VELOCITIES[i][2]);
                            collide_field[Q_LBM * (nx + y * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (nx_b1 + y * xmax + z * xmax * ymax) + i];

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
                if (bcd[x + y * xmax + z * xmax * ymax] == OUTLET) {
                    if (y == 0) {
                        ny = 1;
                        ny_b1 = 2;
                        for (int i = 0; i < Q_LBM; i++) {

#ifdef SWICH_OUTLETBC1
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
                            // int nx_b1 = x - int(LATTICE_VELOCITIES[i][0]);
                            // int ny_b1 = y - int(LATTICE_VELOCITIES[i][1]);
                            // int nz_b1 = z - int(LATTICE_VELOCITIES[i][2]);
                            // cout << "OUTLET y ny, ny_b1, x, y =" << y << " " << ny << " " << ny_b1 << " " << x << " " << y << " " << endl;
                            collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (x + ny_b1 * xmax + z * xmax * ymax) + i];

#endif
                        }
                    } else if (y == ymax - 1) {
                        ny = ymax - 2;
                        ny_b1 = ymax - 3;

                        for (int i = 0; i < Q_LBM; i++) {

#ifdef SWICH_OUTLETBC1
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
                            // int nx_b1 = x - int(LATTICE_VELOCITIES[i][0]);
                            // int ny_b1 = y - int(LATTICE_VELOCITIES[i][1]);
                            // int nz_b1 = z - int(LATTICE_VELOCITIES[i][2]);
                            // cout << "OUTLET y ny, ny_b1, x, y =" << y << " " << ny << " " << ny_b1 << " " << x << " " << y << " " << endl;
                            collide_field[Q_LBM * (x + ny * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (x + ny_b1 * xmax + z * xmax * ymax) + i];

#endif
                        }
                    }
                }
            }
        }
    }
#endif

    for (int z = 0; z < zmax; z++) {
        for (int y = 0; y < ymax; y++) {
            for (int x = 0; x < xmax; x++) {
                if (bcd[x + y * xmax + z * xmax * ymax] == PERIODIC) {

                    for (int i = 0; i < Q_LBM; i++) {
                        if (y == ymax - 1) {
                            nyperimin = 2;
                            collide_field[Q_LBM * (x + (y - 1) * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (x + nyperimin * xmax + z * xmax * ymax) + i];

                        } else if (y == 0) {
                            nyperimax = ymax - 3;
                            collide_field[Q_LBM * (x + (y + 1) * xmax + z * xmax * ymax) + i] = collide_field[Q_LBM * (x + nyperimax * xmax + z * xmax * ymax) + i];
                        }
                    }
                }
            }
        }
    }
}