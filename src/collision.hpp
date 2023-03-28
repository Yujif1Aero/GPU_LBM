#ifndef _COLLISION_H_
#define _COLLISION_H_

/**
 * Carries out the whole local collision process. Computes density and velocity and equilibrium
 * distributions. Carries out BGK update.
 */
void DoCollision(float *collide_field, float tau, const int xstep, const int ystep,const int zstep, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend);

#endif
