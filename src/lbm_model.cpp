#include <cmath>
#include <cstdio>

#include "lbm_model.hpp"
#include "utils.hpp"


void ValidateModel(float* velocity_bc, float RefLength, float tau, float Rho) {
    float mach_number, reynolds_number;
    float nu = C_S_POW2 * (tau - 0.5);
    float nuinv = 1.0 / nu;
    float abs_vel_bc = sqrt((velocity_bc[0]) * (velocity_bc[0]) +
                      (velocity_bc[1]) * (velocity_bc[1]) +
                      (velocity_bc[2]) * (velocity_bc[2]));
    mach_number = abs_vel_bc * C_S_INV;
    reynolds_number = abs_vel_bc * RefLength * nuinv;
    printf("Computed Mach number: %f\n", mach_number);
    printf("Computed Reynolds number: %f\n", reynolds_number);
    printf("Computed relaxation time: %f\n", tau);

    /* Check if characteristic numbers are correct */
    if (mach_number >= 1)
        ERROR("Computed Mach number is too large.");
}
