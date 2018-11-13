#ifndef BHMIE_H
#define BHMIE_H

#include <complex>

using namespace std;

// Function Declarations
int bhmie(double x, complex<double> refrel, int nang, double* Qscat_p, double* Qext_p, double* Qabs_p, double* Qback_p, complex<double>* S1_p, complex<double>* S2_p); // mie calculations

#endif
