#ifndef PHOTONTRACKING_H
#define PHOTONTRACKING_H

using namespace std;

// Function Declarations

// Update X Direction Cosine
double updateDirCosX(double theta, double phi, double mux, double muy, double muz);

// Update Y Direction Cosine
double updateDirCosY(double theta, double phi, double mux, double muy, double muz);

// Update Z Direction Cosine
double updateDirCosZ(double theta, double phi, double mux, double muy, double muz);

// calculate the rotation angle
double gammaCalc(double muz1, double muz2, double theta, double phi);

#endif
