#include <iostream>
#include <vector>
using namespace std;
#include "collision.hpp"
#include "cell_computation.hpp"
#include "lbm_model.hpp"
#include "utils.hpp"
/**
 * Computes the post-collision distribution functions according to the BGK
 * update rule and stores the results again at the same position.
 */
void ComputePostCollisionDistributions(float* current_cell, float tau, const float* const feq) {
    int i;
    float tauinv = 1.0 / tau;
    for (i = 0; i < Q_LBM; i++) {
        current_cell[i] = current_cell[i] - (current_cell[i] - feq[i]) * tauinv;

        /* Probability distribution function can not be less than 0 */
        if (current_cell[i] < 0)
            ERROR("Probability distribution function can not be negative.");
    }
}

void DoCollision(float* collide_field, float tau, const int xmax, const int ymax, const int zmax, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend, vector<int> &bcd) {
    float density, velocity[3], feq[Q_LBM], *currentCell;
    int x, y, z;
    for (z = zstart; z <= zend; z++) {
        for (y = ystart; y <= yend; y++) {
            for (x = xstart; x <= xend; x++) {
               if( bcd[x + y * xmax + z * xmax * ymax] == FLUID) {
                currentCell = &collide_field[Q_LBM * (x + y * xmax + z * xmax * ymax)];

                ComputeDensity(currentCell, &density);
                ComputeVelocity(currentCell, &density, velocity);
                ComputeFeq(&density, velocity, feq);
                ComputePostCollisionDistributions(currentCell, tau, feq);
               }
            }
        }
    }
}
