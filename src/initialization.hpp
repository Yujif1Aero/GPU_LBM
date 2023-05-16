#ifndef _INITIALIZATION_H_
#define _INITIALIZATION_H_
#include "lbm_model.hpp"

/**
 * Reads the parameters for the lid driven cavity scenario from a configuration
 * file
 */
void ReadParameters(float* LengthRef, float* tau, float *nu, float *VelInt, float* RhoInt, float* PresInst, float* velocity_bc, float *pressure_bc, int* timesteps, int* timesteps_per_plotting, char* dirmesh, int argc, char* argv[], int* gpu_enabled);


/**
 * Initializes the particle distribution functions and the flag field
 */
void InitialiseFields(float* collide_field, float* stream_field, int xmax, int ymax, int zmax, int gpu_enabled, float Feq[Q_LBM], float RhoInt, float* VelInt, vector<int>& bcd);


void InitialiseGrid(char *dirmesh, int& jmax, int& kmax, int& lmax, int& xstart, int& ystart, int& zstart, int& xend, int& yend, int& zend, vector<double>& xd, vector<double>& yd, vector<double>& zd, vector<int>& bcd);

#endif
