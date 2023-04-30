#ifndef _MAIN_C_
#define _MAIN_C_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

using namespace std;
#include "boundary.hpp"
#include "collision.hpp"
#include "initialization.hpp"
#include "lbm_model.hpp"
#include "streaming.hpp"
#include "utils.hpp"
#include "visualization.hpp"
#include "stopwatch.h"

int main(int argc, char* argv[]) {
	float* collide_field = NULL, * stream_field = NULL, * swap = NULL, tau, nu,num_cells, mlups_sum;

	/* Boundary condition parameter */
	//float wall_velocity[D_LBM]; //for cavity
	float bc_velocity[D_LBM];
	float bc_pressure;
	int t = 0, timesteps, timesteps_per_plotting, gpu_enabled;

	/*GRID parameter*/
	int xmax, ymax, zmax, xstart, ystart, zstart, xend, yend, zend;
	// int xmax, ymax, zmax, xstart, ystart, zstart, xend, yend, zend; // for cavity
	float LengthRef;
	//int xlength; // for cavity
	vector<double> xd, yd, zd;

	/* INITIAL parameter*/
	//float RhoInt = 1.0, VelInt = 0.0, PresInt = C_S2 * RhoInt, Feq[Q_LBM]; //for cavity
	
	//float RhoInt = 1.0, VelInt[D_LBM], PresInt =  C_S2 * RhoInt, Feq[Q_LBM];
    float RhoInt, VelInt[D_LBM], PresInt, Feq[Q_LBM];

    clock_t mlups_time;
	size_t field_size;
    
	/* process parameters */
	//ReadParameters(&xlength, &tau, VelInt, wall_velocity, &timesteps,
	//	&timesteps_per_plotting, argc, argv, &gpu_enabled); 	// for cavity

	ReadParameters(&LengthRef, &tau, &nu, VelInt, &RhoInt, &PresInt, bc_velocity, &bc_pressure, &timesteps, &timesteps_per_plotting, argc, argv, &gpu_enabled);
	
	/* check if provided parameters are legitimate */
	//ValidateModel(bc_velocity, LengthRef, tau, VelInt, RhoInt); //for cavity
	ValidateModel(bc_velocity, LengthRef, tau, RhoInt); 

	/* 	 grid generation */
	// InitialiseGrid(xlength, xmax, ymax, zmax, xstart, ystart, zstart, xend, yend, zend); // for cavity
	InitialiseGrid(xmax, ymax, zmax , xstart, ystart, zstart, xend, yend, zend, xd, yd, zd);


	/* initializing fields */
	// num_cells = pow(xlength + 2, D_LBM);
	num_cells = xmax * ymax * zmax;
	field_size = Q_LBM * num_cells * sizeof(float);
	collide_field = (float*)malloc(field_size);
	stream_field = (float*)malloc(field_size);
	// flag_field = (int*)malloc(num_cells * sizeof(int));

	/* for cavity */
	// InitialiseFields(collide_field, stream_field, xmax, ymax, zmax, gpu_enabled, Feq, RhoInt, wall_velocity);

	// writeFlagField(flag_field, "../Initial/init", xmax, ymax, zmax, gpu_enabled);
	// TreatBoundary(collide_field, wall_velocity, xstart, ystart,
	//               zstart, xend, yend, zend, xmax, ymax, zmax);

    InitialiseFields(collide_field, stream_field, xmax, ymax, zmax, gpu_enabled,Feq, RhoInt, VelInt);
    WriteFluidVtkOutput(collide_field, "../Initial/initfluidfiled", t, xstart, ystart,
                        zstart, xend, yend, zend, xmax, ymax, zmax);
    //WriteAllVtkOutput(collide_field, "../Initial/initallfiled", t, xmax, ymax, zmax);

	WriteField(collide_field, "../Initial/field", t, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax,
		gpu_enabled);
	WritePhysics(collide_field, "../Initial/physics-filed", t, xstart, ystart,
		zstart, xend, yend, zend, xmax, ymax, zmax, gpu_enabled);
	mlups_sum = 0.f;
	// for (t = 0; t <= timesteps; t++) {
	// 	printf("Time step: #%d\r", t);
	// 	Stopwatch sw;
        //      sw.start();
	// 	/* Compute post collision distributions */
	// 	DoCollision(collide_field, tau, xmax, ymax, zmax, xstart, ystart, zstart, xend, yend, zend);
	// 	/* Copy pdfs from neighbouring cells into collide field */
	// 	DoStreaming(collide_field, stream_field, xmax, ymax, zmax, xstart, ystart, zstart, xend, yend, zend);
	// 	/* Perform the swapping of collide and stream fields */
	// 	swap = collide_field;
	// 	collide_field = stream_field;
	// 	stream_field = swap;
	// 	/* Treat boundaries */
	// 	TreatBoundary(collide_field, wall_velocity, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax);
	//      sw.stop();
        //      mlups_time = sw.getMs();
	// 	/* Print out the MLUPS value */
	// 	mlups_sum += num_cells / (MLUPS_EXPONENT * (float)mlups_time ;
	// 	//if (VERBOSE)
	// 	printf("MLUPS: %f\n", num_cells / (MLUPS_EXPONENT * (float)mlups_time));
	// 	/* Print out vtk output if needed */
	// 	if (!(t % timesteps_per_plotting)) {
	// 		WriteFluidVtkOutput(collide_field, "../img/all-fielda", t, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax);
	// 		WriteFluidVtkOutput(collide_field, "../img/fluid-field", t, xstart + 1, ystart + 1, zstart, xend - 1, yend - 1, zend, xmax, ymax, zmax);
	// 	}
	// }

	// printf("Average MLUPS: %f\n", mlups_sum / (t + 1));

	// if (VERBOSE) {
	// 	WriteField(collide_field, "../Final/disfunc-field", t, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax,
	// 		gpu_enabled);
	// 	WritePhysics(collide_field, "../Final/physics-filed", t, xstart, ystart, zstart, xend, yend, zend, xmax, ymax, zmax, gpu_enabled);

	// }

	// /* Free memory */
	// free(collide_field);
	// free(stream_field);

	// // FreeDeviceFields(&collide_field_d, &stream_field_d, &flag_field_d);

	printf("Simulation complete.\n");
	return 0;

}

#endif
