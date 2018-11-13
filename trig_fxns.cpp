#include "trig_fxns.hpp"

using namespace std;

// Convert Radians to Degrees
double rad2Deg(double xrad) {
    double xdeg = xrad*(180/pi);
  return xdeg;
}


// Convert Degrees to Radians
double deg2Rad(double xdeg) {
    double xrad = xdeg *(pi/180);
  return xrad;
}
