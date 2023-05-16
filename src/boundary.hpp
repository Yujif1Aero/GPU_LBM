#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_


/**
 * Handles the boundaries in our simulation setup
 */

void TreatBoundary(float *collide_field, float* velocity_bc, int xstart, int ystart, int zstart, int xend, int yend, int zend, int xmax, int ymax, int zmax, vector<int>& bcd, float RhoInt);


#endif
