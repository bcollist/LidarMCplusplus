#include <cmath>

#include "photon_tracking.hpp"
#include "constants.hpp"


using namespace std;


// Function Declarations
// Update X Direction Cosine
double updateDirCosX(double theta, double phi, double mux, double muy, double muz) {
  double muxPrime;
    if (abs(muz) > 0.999){
      muxPrime = sin(theta) * cos(phi);
    }
        else{
        muxPrime = (1.0/(sqrt(1-muz*muz))) * sin(theta) * (mux * muz * cos(phi)-muy * sin(phi)) + mux * cos(theta);
      }
  return muxPrime;
}

// Update Y Direction Cosine
double updateDirCosY(double theta, double phi, double mux, double muy, double muz) {
  double muyPrime;
    if (abs(muz) > 0.999){
      muyPrime = sin(theta) * sin(phi);
    }
        else{
        muyPrime = (1.0/(sqrt(1-muz*muz))) * sin(theta) * (muy*muz*cos(phi)+mux*sin(phi)) + muy * cos(theta);
      }
  return muyPrime;
}

// Update Z Direction Cosine
double updateDirCosZ(double theta, double phi, double mux, double muy, double muz) {
  double muzPrime;
    if (abs(muz) > 0.999){
      muzPrime = cos(theta) * muz/abs(muz);
    }
        else{
         muzPrime = (-1*sqrt(1-muz*muz)) * sin(theta) * cos(phi) + muz*cos(theta);
       }
  return muzPrime;
}

// Calculate the rotation angle into the new reference frame of the photon
double gammaCalc(double muz1, double muz2, double theta, double phi){

    // Calculates the angle of rotation back into the new photon coordinate space using sperical trig cosine identity
    // cosa = cos(b)cos(c) + sin(b)sin(c)cos(A)
    // set the unknown angle to be A, rearrange, and solve
    // See http://mathworld.wolfram.com/SphericalTrigonometry.html for more details if you need some visual help

    // Variable Initialization
    double gammaCos; double gamma;

    //Function Body

    if (pi < phi < 2*pi){
        gammaCos = (muz1 - muz2*cos(theta)) / sin(theta) * sqrt((1 - muz2 * muz2));
    }
    else{
        gammaCos = (muz1 - muz2*cos(theta)) / (-1 * sin(theta)) * sqrt((1 - muz2 * muz2));

    }
    gamma = acos(gammaCos);

    return gamma;
}
