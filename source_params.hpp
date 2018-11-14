#ifndef SOURCEPARAMS_H
#define SOURCEPARAMS_H

double xl; double yl; double zl; // position of the laser source in (m)
double Isrc; double Qsrc; double Usrc; double Vsrc;
double muxsrc; double muysrc; double muzsrc;

double lambda; // lidar wavelength in a vaccuum (um)


// double pulseWidth = 7; // gaussian pulse width FWHM (ns)
// // double pulseWidth_m = pulseWidth * Cwater; // gaussian pulse width FWHM (m)
// double beamRad = 0.5E-2; // beam radius (m)
// double beamDiv = 2E-3; // beam divergence (radians)

#endif
