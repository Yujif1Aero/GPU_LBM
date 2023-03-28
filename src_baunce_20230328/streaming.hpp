#ifndef _STREAMING_H_
#define _STREAMING_H_

/**
 * Carries out the streaming step and writes the respective distribution functions from
 * collideField to streamField.
 */
void DoStreaming(float *collide_field, float *stream_field, const int xstep, const int ystep, const int zstep, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend);

#endif
