#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std;

#include "cell_computation.hpp"
#include "lbm_model.hpp"
#include "utils.hpp"
#include "visualization.hpp"

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


void write_vtkPointCoordinates(FILE* fp, int idomein, vector<double>& xd, vector<double>& yd, vector<double>& zd, vector<int>& bcd) {


    for (int i = 0; i <= idomein - 1; i++) {
        if (bcd[i] == FLUID)
            fprintf(fp, "%f %f %f\n", xd[i], yd[i], zd[i]);
    }
}

void write_vtkPointCoordinates(FILE* fp, int idomein, vector<double>& xd, vector<double>& yd, vector<double>& zd) {


    for (int i = 0; i < idomein; i++) {
        fprintf(fp, "%f %f %f\n", xd[i], yd[i], zd[i]);
    }
}


void WriteAllVtkOutput(const float* const collideField, const char* filename, unsigned int t, int xmax, int ymax, int zmax, vector<double>& xd, vector<double>& yd, vector<double>& zd) {
    int x, y, z;
    int idomein = xmax * ymax * zmax;
    float velocity[3], density, pressure;

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

    write_vtkHeader(fp, xmax, ymax, zmax);
    write_vtkPointCoordinates(fp, idomein, xd, yd, zd);

    fprintf(fp, "POINT_DATA %i \n", xmax * ymax * zmax);
    fprintf(fp, "\n");
    fprintf(fp, "VECTORS velocity float\n");
    for (z = 0; z < zmax; z++) {
        for (y = 0; y < ymax; y++) {
            for (x = 0; x < xmax; x++) {
                ComputeDensity(
                    &collideField[Q_LBM * (x + y * xmax + z * xmax * ymax)],
                    &density);
                ComputeVelocity(
                    &collideField[Q_LBM * (x + y * xmax + z * xmax * ymax)],
                    &density, velocity);
                fprintf(fp, "%f %f %f\n", velocity[0], velocity[1],
                        velocity[2]);
            }
        }
    }

    fprintf(fp, "\n");
    fprintf(fp, "SCALARS density float 1 \n");
    fprintf(fp, "LOOKUP_TABLE default \n");
    for (z = 0; z < zmax; z++) {
        for (y = 0; y < ymax; y++) {
            for (x = 0; x < xmax; x++) {
                ComputeDensity(
                    &collideField[Q_LBM * (x + y * xmax + z * xmax * ymax)],
                    &density);
                fprintf(fp, "%f\n", density);
            }
        }
    }

    fprintf(fp, "\n");
    fprintf(fp, "SCALARS pressure float 1 \n");
    fprintf(fp, "LOOKUP_TABLE default \n");
    for (z = 0; z < zmax; z++) {
        for (y = 0; y < ymax; y++) {
            for (x = 0; x < xmax; x++) {
                ComputeDensity(
                    &collideField[Q_LBM * (x + y * xmax + z * xmax * ymax)],
                    &density);
                pressure = density * C_S_POW2;
                fprintf(fp, "%f\n", pressure);
            }
        }
    }
    if (fclose(fp)) {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", szFileName);
        ERROR(szBuff);
    }
}



void WriteFluidVtkOutput(const float* const collideField, const char* filename, unsigned int t, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend, const int xmax, const int ymax, const int zmax, vector<double>& xd, vector<double>& yd, vector<double>& zd, vector<int>& bcd) {
    int x, y, z;
    
    int xrange = xend - xstart + 1;
    int yrange = yend - ystart + 1;
    int zrange = zend - zstart + 1;

    int idomein = xmax * ymax * zmax;

    float velocity[3], density, pressure;

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
    write_vtkPointCoordinates(fp, idomein, xd, yd, zd, bcd);

    fprintf(fp, "POINT_DATA %i \n", xrange * yrange * zrange);
    fprintf(fp, "\n");
    fprintf(fp, "VECTORS velocity float\n");
    for (z = zstart; z <= zend; z++) {
        for (y = ystart; y <= yend; y++) {
            for (x = xstart; x <= xend; x++) {
                ComputeDensity(
                    &collideField[Q_LBM * (x + y * xmax + z * xmax * ymax)],
                    &density);
                ComputeVelocity(
                    &collideField[Q_LBM * (x + y * xmax + z * xmax * ymax)],
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
                    &collideField[Q_LBM * (x + y * xmax + z * xmax * ymax)],
                    &density);
                fprintf(fp, "%f\n", density);
            }
        }
    }

    fprintf(fp, "\n");
    fprintf(fp, "SCALARS pressure float 1 \n");
    fprintf(fp, "LOOKUP_TABLE default \n");
    for (z = zstart; z <= zend; z++) {
        for (y = ystart; y <= yend; y++) {
            for (x = xstart; x <= xend; x++) {
                ComputeDensity(
                    &collideField[Q_LBM * (x + y * xmax + z * xmax * ymax)],
                    &density);
                pressure = density * C_S_POW2;
                fprintf(fp, "%f\n", pressure);
            }
        }
    }
    if (fclose(fp)) {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", szFileName);
        ERROR(szBuff);
    }
}

void PrintField(float* field, int xmax, int ymax, int zmax) {
    int x, y, z, i;

    for (z = 0; z < zmax; z++) {
        for (y = 0; y < ymax; y++) {
            for (x = 0; x < xmax; x++) {
                printf("(%d,%d,%d): ", x, y, z);
                for (i = 0; i < Q_LBM; i++)
                    printf(
                        "%f ",
                        field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i]);
                printf("\n");
            }
        }
    }
}

void WriteField(const float* const field, const char* filename, unsigned int t,
                const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend, const int xmax, const int ymax, const int zmax, const int rank) {
    int x, y, z, i;
    char szFileName[80];
    FILE* fp = NULL;
    sprintf(szFileName, "%s-field%i.%i.out", filename, rank, t);
    fp = fopen(szFileName, "w");
    if (fp == NULL) {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", szFileName);
        ERROR(szBuff);
        return;
    }

    for (z = 0; z < zmax; z++) {
        for (y = 0; y < ymax; y++) {
            for (x = 0; x < xmax; x++) {
                fprintf(fp, "(%d,%d,%d): ", x, y, z);
                for (i = 0; i < Q_LBM; i++)
                    fprintf(
                        fp, "%f ",
                        field[Q_LBM * (x + y * xmax + z * xmax * ymax) + i]);
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

void WritePhysics(const float* const field, const char* filename, unsigned int t, const int xstart, const int ystart, const int zstart, const int xend, const int yend, const int zend, const int xmax, const int ymax, const int zmax, const int rank) {
    int x, y, z;
    float density, velocity[3], pressure;

    char szFileName[80];
    FILE* fp = NULL;
    sprintf(szFileName, "%s-pysics%i_%i.csv", filename, rank, t);
    fp = fopen(szFileName, "w");

    fopen(szFileName, "w");
    if (fp == NULL) {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", szFileName);
        ERROR(szBuff);
        return;
    }
    fprintf(fp, "x, y, z, vx, vy, vz, rho, pressure");
    fprintf(fp, "\n");
    for (z = zstart; z <= zend; z++) {
        for (y = ystart; y <= yend; y++) {
            for (x = xstart; x <= xend; x++) {
                ComputeDensity(
                    &field[Q_LBM * (x + y * xmax + z * xmax * ymax)],
                    &density);
                ComputeVelocity(
                    &field[Q_LBM * (x + y * xmax + z * xmax * ymax)],
                    &density, velocity);
                pressure = density * C_S_POW2;
                fprintf(fp, "%d, %d, %d, %f, %f, %f. %f, %f", x, y, z, velocity[0],
                        velocity[1], velocity[2], density, pressure);
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

void writebcd(vector<int>& bcd, const char* filename, const int xmax, const int ymax, const int zmax, const int rank) {
    int x, y, z;
    char szFileName[80];
    FILE* fp = NULL;
    sprintf(szFileName, "%s-bcd%i.out", filename, rank);
    fp = fopen(szFileName, "w");
    if (fp == NULL) {
        char szBuff[80];
        sprintf(szBuff, "Failed to open %s", szFileName);
        ERROR(szBuff);
        return;
    }

    fprintf(fp, "(y,z)\\x  ");
    for (x = 0; x < xmax; x++) {
        fprintf(fp, "%2i ", x);
    }
    fprintf(fp, "\n-------");
    for (x = 0; x < xmax; x++) {
        fprintf(fp, "---");
    }
    fprintf(fp, "\n");

    for (z = 0; z < zmax; z++) {
        for (y = ymax - 1; y >= 0; y--) {
            fprintf(fp, "(%2d,%2d): ", y, z);
            for (x = 0; x < xmax; x++) {
                fprintf(fp, "%2i ", bcd[x + y * xmax + z * xmax * ymax]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }

    if (fclose(fp)) {
        char szBuff[80];
        sprintf(szBuff, "Failed to close %s", szFileName);
        ERROR(szBuff);
    }
}
