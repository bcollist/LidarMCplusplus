#include <cstdlib>
#include <iostream>
#include <complex>

#include "constants.hpp"
#include "bhmie.hpp"

using namespace std;




int main(){
  double x = 10.3825;
  double refRe = 1.05;
  double refIm = 0.0000;
  complex<double> refrel = complex<double>(refRe,refIm);
  int nang = 455;

  double Qscat; double Qext; double Qabs; double Qback;  // Mie scattering efficiencies
  double* Qscat_p = &Qscat; double* Qext_p = &Qext; double* Qabs_p = &Qabs; double* Qback_p = &Qback; // pointers to Mie Scattering efficiencies
  double Qb; double Qc; double Qa; double Qba;
  complex<double> S1[2*nang]; complex<double> S2[2*nang];
  complex<double>* S1_p = S1; complex<double>* S2_p = S2;

  bhmie(x, refrel, nang, Qscat_p, Qext_p, Qabs_p, Qback_p,  S1_p,  S2_p);
}
