#include "streaming.hpp"
#include "lbm_model.hpp"
#include <iostream>
using namespace std;

void DoStreaming(float* collide_field, float* stream_field, const int xstep, const int ystep, const int zstep, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend) {
    int x, nx, y, ny, z, nz, i;

    for (z = zstart; z <= zend; z++) {
        for (y = ystart; y <= yend; y++) {
            for (x = xstart; x <= xend; x++) {

                for (i = 0; i < Q_LBM; i++) {
                    // nx = x -  LATTICE_VELOCITIES[i][0];
                    // ny = y -  LATTICE_VELOCITIES[i][1];
                    // nz = z -  LATTICE_VELOCITIES[i][2];

                    nx = x + LATTICE_VELOCITIES[i][0];
                    ny = y + LATTICE_VELOCITIES[i][1];
                    nz = z + LATTICE_VELOCITIES[i][2];

                    stream_field[Q_LBM * (nx + ny * ystep + nz * zstep * zstep) + i] =
                        collide_field[Q_LBM * (x + y * ystep + z * zstep * zstep) + i];
                }
            }
        }
    }
}
