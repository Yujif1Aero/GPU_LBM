#include "boundary.hpp"
#include "cell_computation.hpp"
#include "lbm_model.hpp"
#include <iostream>
using namespace std;
/**
 * Finds an inverse probability distribution for the provided lattice index.
 */
int inv(int i) {
    return (Q_LBM - 1) - i;
}

void TreatBoundary(float* collide_field, float* wall_velocity, int xstart, int ystart, int zstart, int xend, int yend, int zend, int xstep, int ystep, int zstep) {
    int x, nx, y, ny, z, nz, i;
    //int invC = int(1 / C);
    float dot_prod;
    int imxy = 7, imx0 = 2, imxmy = 0;

    //  bounce back scheme
    for (z = zstart; z <= zend; z++) {
        for (y = ystart; y <= yend; y++) {
            // non-slip boundary left wall
            // nx = 0; //dummy wall
            nx = xstart; // wall
            i = imxy;
            /* Assign the boudary cell value */
            collide_field[Q_LBM * (nx + y * ystep + z * zstep * zstep) + inv(i)] = collide_field[Q_LBM * (nx + y * ystep + z * zstep * zstep) + i];

            i = imx0;
            /* Assign the boudary cell value */
            collide_field[Q_LBM * (nx + y * ystep + z * zstep * zstep) + inv(i)] = collide_field[Q_LBM * (nx + y * ystep + z * zstep * zstep) + i];

            i = imxmy;
            /* Assign the boudary cell value */
            collide_field[Q_LBM * (nx + y * ystep + z * zstep * zstep) + inv(i)] = collide_field[Q_LBM * (nx + y * ystep + z * zstep * zstep) + i];
        }
        for (y = ystart; y <= yend; y++) {
            // non-slip boundary right  wall
            // nx = xstep -1; //dummy wall
            nx = xend; // wall
            i = inv(imxy);
            /* Assign the boudary cell value */
            collide_field[Q_LBM * (nx + y * ystep + z * zstep * zstep) + inv(i)] = collide_field[Q_LBM * (nx + y * ystep + z * zstep * zstep) + i];

            i = inv(imx0);
            /* Assign the boudary cell value */
            collide_field[Q_LBM * (nx + y * ystep + z * zstep * zstep) + inv(i)] = collide_field[Q_LBM * (nx + y * ystep + z * zstep * zstep) + i];

            i = inv(imxmy);
            /* Assign the boudary cell value */
            collide_field[Q_LBM * (nx + y * ystep + z * zstep * zstep) + inv(i)] = collide_field[Q_LBM * (nx + y * ystep + z * zstep * zstep) + i];
        }
    }

    for (z = zstart; z <= zend; z++) {
        for (x = xstart; x <= xend; x++) {
 
            // no-slip boundary bottom wall
            // ny = 0; //dummy wall
            ny = ystart; // wall
            for (i = 0; i < 4; i++) {

                collide_field[Q_LBM * (x + ny * ystep + z * zstep * zstep) + inv(i)] = collide_field[Q_LBM * (x + ny * ystep + z * zstep * zstep) + i];
            }
 
            // moving-slip boundary top wall
            // ny = ystep -1; // dummy wall
            ny = yend; // wall
            for (i = 6; i < Q_LBM; i++) {
                dot_prod = LATTICE_VELOCITIES[i][0] * wall_velocity[0] + LATTICE_VELOCITIES[i][1] * wall_velocity[1] + LATTICE_VELOCITIES[i][2] * wall_velocity[2];
                /* Assign the boudary cell value */
                collide_field[Q_LBM * (x + ny * ystep + z * zstep * zstep) + inv(i)] =
                    collide_field[Q_LBM * (x + ny * ystep + z * zstep * zstep) + i] - 6.0 * LATTICE_WEIGHTS[i] * dot_prod;
            }
        }
    }
}