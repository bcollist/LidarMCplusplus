#ifndef DETECTORPARAMS_H
#define DETECTORPARAMS_H

#include "trig_fxns.hpp"

double detectorRad = 1.5E-1; // detector radius
FOV=deg2Rad(FOV); // convert FOV from degrees to radians

double xd = 0.0; double yd = 0.0; double zd = 0.0; // position of the detector in (m)
double fd; // variable used in detector photon geometry colculations
double anglei; // angle of intersection between photon and detector plane

#endif
