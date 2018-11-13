#include <cmath>
#include <iostream>
#include <complex>

#include "constants.hpp"
#include "bhmie.hpp"

using namespace std;

int main(){
  double x;
  complex<double> refrel;
  int nang;

  x = 2.0;
  refrel = complex<double>(1.05,0.0005); // refRel stores the refractive index as a complex double
  nang = 455; // number of angles between 0-90

  double Qscat; double Qext; double Qabs; double Qback;  // Mie scattering efficiencies
  double* Qscat_p = &Qscat; double* Qext_p = &Qext; double* Qabs_p = &Qabs; double* Qback_p = &Qback; // pointers to Mie Scattering efficiencies
  complex<double> S1[2*nang]; complex<double> S2[2*nang];
  complex<double>* S1_p = S1; complex<double>* S2_p = S2;


  bhmie(x, refrel, nang, Qscat_p, Qext_p, Qabs_p, Qback_p, S1_p, S2_p);
  cout << S1[10] << endl;
}
