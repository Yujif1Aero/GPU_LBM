#ifndef _INITIALIZATION_H_
#define _INITIALIZATION_H_
#include "lbm_model.hpp"

/**
 * Reads the parameters for the lid driven cavity scenario from a configuration
 * file
 */
void ReadParameters(int* xlength, float* tau, float& VelInt, float* velocity_wall, int* timesteps, int* timesteps_per_plotting, int argc, char* argv[], int* gpu_enabled);

/**
 * Initializes the particle distribution functions and the flag field
 */
void InitialiseFields(float* collide_field, float* stream_field, int xstep, int ystep, int zstep, int gpu_enabled, float Feq[Q_LBM], float RhoInt, float wall_velocity[D_LBM]);

void InitialiseGrid(int xlength, int& xstep, int& ystep, int& zstep, int& xstart, int& ystart, int& zstart, int& xend, int& yend, int& zend);

#endif
