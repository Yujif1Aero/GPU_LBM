#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

int main() {
    int jmax, kmax, lmax, idomein;    // number of grid points in xi, eta and zeta direction
    double xl, yl, zl;                // horizontal, vertical, depth computational dommain
    double dnx, dny, dnz, dx, dy, dz; // 1/dnx = dx (grid spacing)
    double unit = 1.0;
    int index;
    int ibinary = 0;

    // parameters set
    // jmax = 201;
    // kmax = 101;
    // lmax = 1;
    xl = 3.0;
    yl = 1.0;
    zl = 0.0;
    dnx = 50.0;
    dny = 100.0;
    dnz = 1.0;
    dx = unit / dnx;
    dy = unit / dny;
    dz = unit / dnz;
    printf("Grid Space: dx = %lf, dy = %lf, dz = %lf", dx, dy, dz);
    jmax = dnx * xl + 1;
    kmax = dny * yl + 1;
    lmax = dnz * zl + 1;

    idomein = jmax * kmax * lmax;
    vector<double> x(idomein, 0.0);
    vector<double> y(idomein, 0.0);
    vector<double> z(idomein, 0.0);
    printf("Number of Grid Poits: jmax = %d, kmax = %d, lmax = %d, total = %d", jmax, kmax, lmax, idomein);

    for (int l = 0, index = 0; l <= lmax - 1; l++) {
        for (int k = 0; k <= kmax - 1; k++) {
            for (int j = 0; j <= jmax - 1; j++, index++) {
                x[index] = dx * double(j);
                y[index] = dy * double(k);
                z[index] = dz * double(l);
            }
        }
    }

    if (ibinary == 0) {
        ofstream file("grid.dat");
        if (!file.is_open()) {
            cerr << "failed to open grid.dat" << endl;
            return -1;
        }
        // ascii
        file << "TITLE = \" mesh data \"" << endl;
        file << "VARIABLES = \"X\", \"Y\", \"Z\"" << endl;
        file << "ZONE T = \"GRID\", I=" << jmax << ", J=" << kmax << ", K=" << lmax << ", F=POINT" << endl;
        for (int l = 0, index = 0; l <= lmax - 1; l++) {
            for (int k = 0; k <= kmax - 1; k++) {
                for (int j = 0; j <= jmax - 1; j++, index++) {
                    file << x[index] << " " << y[index] << " " << z[index] << endl;
                }
            }
        }
    } else {
        // binary
        ofstream file("grid.dat", ios::binary);
        if (!file.is_open()) {
            cerr << "failed to open grid.plt" << endl;
            return -1;
        }
        file << "TITLE = \"mesh data\"" << endl;
        file << "VARIABLES = \"X\", \"Y\", \"Z\"" << endl;
        file << "ZONE T = \"GRID\", I=" << jmax << ", J=" << kmax << ", K=" << lmax << ", F=POINT" << endl;
        // file.write((const char*)&jmax, sizeof(int));
        // file.write((const char*)&kmax, sizeof(int));
        // file.write((const char*)&lmax, sizeof(int));
        file.write((const char*)x.data(), jmax * kmax * lmax * sizeof(double));
        file.write((const char*)y.data(), jmax * kmax * lmax * sizeof(double));
        file.write((const char*)z.data(), jmax * kmax * lmax * sizeof(double));
        file.close();
    }


}
