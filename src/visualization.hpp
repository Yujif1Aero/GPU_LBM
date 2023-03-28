#ifndef _VISUALIZATION_H_
#define _VISUALIZATION_H_

/**
 * Writes the density and velocity field (derived from the distributions in
 * collideField) to a file determined by 'filename' and timestep 't'. You can
 * re-use parts of the code from visual.c (VTK output for Navier-Stokes solver)
 * and modify it for 3D datasets.
 */
void WriteAllVtkOutput(const float* const collide_field, const char* filename, unsigned int t, const int xstep, const int ystep, const int zstep);

void WriteFluidVtkOutput(const float* const collide_field, const char* filename, unsigned int t, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend, const int xstep, const int ystep, const int zstep);

/**
 * Writes the provided stream or collision field to a file with specified
 * identifiers.
 */
void WriteField(const float* const field, const char* filename, unsigned int t,
                const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend, const int xstep, const int ystep, const int zstep,
                const int rank);

/**
 * Writes the provided flag field to a file with specified identifiers.
 */
// void writeFlagField(const int *const flagField, const char *filename,
//                     const int xstep, const int ystep, const int zstep,
//                     const int rank);

/**
 * Prints out the point by point values of the provided field (4D) to the
 * stdout.
 * @param field
 *          linerized 4D array, with
 * (x,y,z,i)=Q*(x+y*(ncell+2)+z*(ncell+2)*(ncell+2))+i
 * @param ncell
 *          number of inner cells, the ones which are there before adding a
 * boundary layer
 */
void PrintField(float* field, int xstep, int ystep, int zstep);

void WritePhysics(const float* const field, const char* filename, unsigned int t, const int xstart,
                  const int ystart, const int zstart, const int xend, const int yend, const int zend, const int xstep, const int ystep, const int zstep, const int rank);
#endif
