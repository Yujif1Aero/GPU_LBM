#include "collision.hpp"
#include "cell_computation.hpp"
#include "lbm_model.hpp"
#include "utils.hpp"
#include <iostream>
using namespace std;

/**
 * Computes the post-collision distribution functions according to the BGK
 * update rule and stores the results again at the same position.
 */
void ComputePostCollisionDistributions(float* current_cell, float tau, const float* const feq) {
    int i;
    for (i = 0; i < Q_LBM; i++) {
        current_cell[i] = current_cell[i] - (current_cell[i] - feq[i]) / tau;

        /* Probability distribution function can not be less than 0 */
        if (current_cell[i] < 0)
            ERROR("Probability distribution function can not be negative.");
    }
}

void DoCollision(float* collide_field, float tau, const int xstep, const int ystep, const int zstep, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend) {
    float density, velocity[3], feq[Q_LBM], *currentCell;
    int x, y, z;
    for (z = zstart; z <= zend; z++) {
        for (y = ystart + 1; y <= yend - 1; y++) {
            for (x = xstart + 1; x <= xend - 1; x++) {
                currentCell = &collide_field[Q_LBM * (x + y * ystep + z * zstep * zstep)];

                ComputeDensity(currentCell, &density);
                ComputeVelocity(currentCell, &density, velocity);
                ComputeFeq(&density, velocity, feq);
                ComputePostCollisionDistributions(currentCell, tau, feq);
            }
        }
    }
}
