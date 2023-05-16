#ifndef _VISUALIZATION_H_
#define _VISUALIZATION_H_



void WriteFluidVtkOutput(const float* const collide_field, const char* filename, unsigned int t, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend, const int xmax, const int ymax, const int zmax, vector<double>& xd, vector<double>& yd, vector<double>& zd, vector<int>& bcd);

void WriteAllVtkOutput(const float* const collide_field, const char* filename, unsigned int t,  const int xmax, const int ymax, const int zmax, vector<double>& xd, vector<double>& yd, vector<double>& zd);
/**
 * Writes the provided stream or collision field to a file with specified
 * identifiers.
 */
void WriteField(const float* const field, const char* filename, unsigned int t, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend, const int xmax, const int ymax, const int zmax, const int rank);

/**
 * Writes the provided flag field to a file with specified identifiers.
 */
void writebcd(vector<int> &bcd, const char *filename, const int xmax, const int ymax, const int zmax, const int rank);


void PrintField(float* field, int xmax, int ymax, int zmax);

void WritePhysics(const float* const field, const char* filename, unsigned int t, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend, const int xmax, const int ymax, const int zmax, const int rank);
#endif
