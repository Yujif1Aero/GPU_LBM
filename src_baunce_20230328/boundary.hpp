#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_


/**
 * Handles the boundaries in our simulation setup
 */
void TreatBoundary(float *collide_field, float *wall_velocity, int xstart, int ystart, int zstart, int xend, int yend, int zend, int xstep, int ystep, int zstep);

#endif
