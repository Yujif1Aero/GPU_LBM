#include "cell_computation.hpp"
#include "lbm_model.hpp"
#include "utils.hpp"
#include <iostream>
using namespace std;

void ComputeDensity(const float* const current_cell, float* density) {
    int i;
    *density = 0.0;
    for (i = 0; i < Q_LBM; i++)
        *density += current_cell[i];
    /* Density should be close to a unit (ρ~1) */
     if ((*density - 1.0) > EPS)
         ERROR("Density dropped below error tolerance.");
}

void ComputeVelocity(const float* const current_cell, const float* const density, float* velocity) {
    int i;
    /* NOTE:Indexes are hard coded to improve program performance */
    velocity[0] = 0.0;
    velocity[1] = 0.0;
    velocity[2] = 0.0;

    for (i = 0; i < Q_LBM; i++) {
        velocity[0] += current_cell[i] * (float)LATTICE_VELOCITIES[i][0];
        velocity[1] += current_cell[i] * (float)LATTICE_VELOCITIES[i][1];
        velocity[2] += current_cell[i] * (float)LATTICE_VELOCITIES[i][2];
    }
    velocity[0] /= *density;
    velocity[1] /= *density;
    velocity[2] /= *density;
    // cout << velocity[0] << endl;
    // cout << velocity[1] << endl;
    // cout << velocity[2] << endl;
}

void ComputeFeq(const float* const density, const float* const velocity, float* feq) {
    int i;
    float dot_prod_cu, dot_prod_cu2, dot_prod_uu; /* summands */
    /* NOTE:Indexes are hard coded to improve program performance */
    for (i = 0; i < Q_LBM; i++) {
        dot_prod_cu = (float)LATTICE_VELOCITIES[i][0] * velocity[0] + (float)LATTICE_VELOCITIES[i][1] * velocity[1] +
                      (float)LATTICE_VELOCITIES[i][2] * velocity[2];
        dot_prod_cu2 = dot_prod_cu * dot_prod_cu;
        dot_prod_uu = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];

        feq[i] = LATTICE_WEIGHTS[i] * (*density) * (1.0 + dot_prod_cu * C_S_POW2_INV + dot_prod_cu2 * C_S_POW4_INV / 2.0 - dot_prod_uu * C_S_POW2_INV / 2.0);
        // cout << feq[i] <<endl;
        /* Probability distribution function can not be less than 0 */
        if (feq[i] < 0)
            ERROR("1Probability distribution function can not be negative.");
    }
}
