#include <iostream>
#include <vector>
using namespace std;
#include "lbm_model.hpp"
#include "streaming.hpp"

void DoStreaming(float* collide_field, float* stream_field, const int xmax, const int ymax, const int zmax, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend, vector<int>& bcd) {
    int x, nx, y, ny, z, nz, i;

    for (z = zstart; z <= zend; z++) {
        for (y = ystart; y <= yend; y++) {
            for (x = xstart; x <= xend; x++) {
                if (bcd[x + y * xmax + z * xmax * ymax] == FLUID) {
                    for (i = 0; i < Q_LBM; i++) {
                        // nx = x -  LATTICE_VELOCITIES[i][0];
                        // ny = y -  LATTICE_VELOCITIES[i][1];
                        // nz = z -  LATTICE_VELOCITIES[i][2];

                        nx = x + LATTICE_VELOCITIES[i][0];
                        ny = y + LATTICE_VELOCITIES[i][1];
                        nz = z + LATTICE_VELOCITIES[i][2];

                        stream_field[Q_LBM * (nx + ny * xmax + nz * xmax * ymax) + i] =
                            collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i];
                    }
                }
            }
        }
    }
}
