#include "visualization.hpp"

#include <cstdio>
#include <iostream>

#include "cell_computation.hpp"
#include "lbm_model.hpp"
#include "utils.hpp"
using namespace std;

void write_vtkHeader(FILE* fp, int xrange, int yrange, int zrange) {
    if (fp == NULL) {
        char szBuff[80];
        sprintf(szBuff, "Null pointer in write_vtkHeader");
        ERROR(szBuff);
        return;
    }

    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "Header\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "\n");
    fprintf(fp, "DATASET STRUCTURED_GRID\n");
    fprintf(fp, "DIMENSIONS  %i %i %i \n", xrange, yrange, zrange);
    fprintf(fp, "POINTS %i float\n", xrange * yrange * zrange);
    fprintf(fp, "\n");
}

void write_vtkPointCoordinates(FILE* fp, int xstart, int ystart, int zstart, int xend, int yend, int zend, int xrange, int yrange, int zrange) {
    float originX = 0.0, originY = 0.0, originZ = 0.0;
    int x = 0, y = 0, z = 0;

    for (z = zstart; z <= zend; z++)
        for (y = ystart; y <= yend; y++)
            for (x = xstart; x <= xend; x++)
                fprintf(fp, "%f %f %f\n", originX + (x * 1.0 / xrange),
                        originY + (y * 1.0 / yrange),
                        originZ + (z * 1.0 / zrange));
}

void write_vtkPointCoordinates(FILE* fp, int xstep, int ystep, int zstep) {
    float originX = 0.0, originY = 0.0, originZ = 0.0;
    int x = 0, y = 0, z = 0;

    for (z = 0; z < zstep; z++)
        for (y = 0; y < ystep; y++)
            for (x = 0; x < xstep; x++)
                fprintf(fp, "%f %f %f\n", originX + (x * 1.0 / xstep),
                        originY + (y * 1.0 / ystep),
                        originZ + (z * 1.0 / zstep));
}

void WriteAllVtkOutput(const float* const collideField, const char* filename, unsigned int t, int xstep, int ystep, int zstep) {
    int x, y, z;

    float velocity[3], density;

    char szFileName[80];
    FILE* fp = NULL;
    sprintf(szFileName, "%s_%i.vtk", filename, t);
    fp = fopen(szFileName, "w");
    if (fp == NULL) {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", szFileName);
        ERROR(szBuff);
        return;
    }

    write_vtkHeader(fp, xstep, ystep, zstep);
    write_vtkPointCoordinates(fp, xstep, ystep, zstep);

    fprintf(fp, "POINT_DATA %i \n", xstep * ystep * zstep);
    fprintf(fp, "\n");
    fprintf(fp, "VECTORS velocity float\n");
    for (z = 0; z < zstep; z++) {
        for (y = 0; y < ystep; y++) {
            for (x = 0; x < xstep; x++) {
                ComputeDensity(
                    &collideField[Q_LBM * (x + y * ystep + z * zstep * zstep)],
                    &density);
                ComputeVelocity(
                    &collideField[Q_LBM * (x + y * ystep + z * zstep * zstep)],
                    &density, velocity);
                fprintf(fp, "%f %f %f\n", velocity[0], velocity[1],
                        velocity[2]);
            }
        }
    }

    fprintf(fp, "\n");
    fprintf(fp, "SCALARS density float 1 \n");
    fprintf(fp, "LOOKUP_TABLE default \n");
    for (z = 0; z < zstep; z++) {
        for (y = 0; y < ystep; y++) {
            for (x = 0; x < xstep; x++) {
                ComputeDensity(
                    &collideField[Q_LBM * (z * zstep * zstep + y * ystep + x)],
                    &density);
                fprintf(fp, "%f\n", density);
            }
        }
    }
    if (fclose(fp)) {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", szFileName);
        ERROR(szBuff);
    }
}

void WriteFluidVtkOutput(const float* const collideField, const char* filename, unsigned int t, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend, const int xstep, const int ystep, const int zstep) {
    int x, y, z;
    // int len = xlength +2; /* lexicographic order "[ Q * ( z*len*len + y*len +
    // x) + i ]" */

    int xrange = xend - xstart + 1;
    int yrange = yend - ystart + 1;
    int zrange = zend - zstart + 1;

    float velocity[3], density;

    char szFileName[80];
    FILE* fp = NULL;
    sprintf(szFileName, "%s_%i.vtk", filename, t);
    fp = fopen(szFileName, "w");
    if (fp == NULL) {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", szFileName);
        ERROR(szBuff);
        return;
    }

    write_vtkHeader(fp, xrange, yrange, zrange);
    write_vtkPointCoordinates(fp, xstart, ystart, zstart, xend, yend, zend, xrange, yrange, zrange);

    fprintf(fp, "POINT_DATA %i \n", xrange * yrange * zrange);
    fprintf(fp, "\n");
    fprintf(fp, "VECTORS velocity float\n");
    for (z = zstart; z <= zend; z++) {
        for (y = ystart; y <= yend; y++) {
            for (x = xstart; x <= xend; x++) {
                ComputeDensity(
                    &collideField[Q_LBM * (x + y * ystep + z * zstep * zstep)],
                    &density);
                // cout << "(x, y ,z, *density) =" <<  x <<"," <<  y  << "," << z << "," << density << endl;
                ComputeVelocity(
                    &collideField[Q_LBM * (x + y * ystep + z * zstep * zstep)],
                    &density, velocity);
                fprintf(fp, "%f %f %f\n", velocity[0], velocity[1],
                        velocity[2]);
            }
        }
    }

    fprintf(fp, "\n");
    fprintf(fp, "SCALARS density float 1 \n");
    fprintf(fp, "LOOKUP_TABLE default \n");
    for (z = zstart; z <= zend; z++) {
        for (y = ystart; y <= yend; y++) {
            for (x = xstart; x <= xend; x++) {
                ComputeDensity(
                    &collideField[Q_LBM * (x + y * ystep + z * zstep * zstep)],
                    &density);
                fprintf(fp, "%f\n", density);
            }
        }
    }
    if (fclose(fp)) {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", szFileName);
        ERROR(szBuff);
    }
}

void PrintField(float* field, int xstep, int ystep, int zstep) {
    int x, y, z, i;

    for (z = 0; z < zstep; z++) {
        for (y = 0; y < ystep; y++) {
            for (x = 0; x < xstep; x++) {
                printf("(%d,%d,%d): ", x, y, z);
                for (i = 0; i < Q_LBM; i++)
                    printf(
                        "%f ",
                        field[Q_LBM * (x + y * ystep + z * zstep * zstep) + i]);
                printf("\n");
            }
        }
    }
}

void WriteField(const float* const field, const char* filename, unsigned int t,
                const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend, const int xstep, const int ystep, const int zstep, const int rank) {
    int x, y, z, i;
    char szFileName[80];
    FILE* fp = NULL;
    sprintf(szFileName, "%s-rank%i.%i.out", filename, rank, t);
    fp = fopen(szFileName, "w");
    if (fp == NULL) {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", szFileName);
        ERROR(szBuff);
        return;
    }

    for (z = zstart; z <= zend; z++) {
        for (y = ystart; y <= yend; y++) {
            for (x = xstart; x <= xend; x++) {
                fprintf(fp, "(%d,%d,%d): ", x, y, z);
                for (i = 0; i < Q_LBM; i++)
                    fprintf(
                        fp, "%f ",
                        field[Q_LBM * (x + y * ystep + z * zstep * zstep) + i]);
                fprintf(fp, "\n");
            }
        }
    }

    if (fclose(fp)) {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", szFileName);
        ERROR(szBuff);
    }
}


void WritePhysics(const float* const field, const char* filename, unsigned int t, const int xstart,
                  const int ystart, const int zstart, const int xend, const int yend, const int zend, const int xstep, const int ystep, const int zstep, const int rank) {
    int x, y, z;
    float density, velocity[3];

    char szFileName[80];
    FILE* fp = NULL;
    sprintf(szFileName, "%s-rank%i_%i.csv", filename, rank, t);
    fp = fopen(szFileName, "w");

    fopen(szFileName, "w");
    if (fp == NULL) {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", szFileName);
        ERROR(szBuff);
        return;
    }
    fprintf(fp, "x, y, z, vx, vy, vz, rho");
    fprintf(fp, "\n");
    for (z = zstart; z <= zend; z++) {
        for (y = ystart; y <= yend; y++) {
            for (x = xstart; x <= xend; x++) {
                 ComputeDensity(
                    &field[Q_LBM * (x + y * ystep + z * zstep * zstep)],
                    &density);
                ComputeVelocity(
                    &field[Q_LBM * (x + y * ystep + z * zstep * zstep)],
                    &density, velocity);
                fprintf(fp, "%d, %d, %d, %f, %f, %f. %f", x, y, z, velocity[0],
                        velocity[1], velocity[2], density);
                fprintf(fp, "\n");
            }
        }
    }
    if (fclose(fp)) {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", szFileName);
        ERROR(szBuff);
    }
}
