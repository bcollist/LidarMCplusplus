//
//  test.cpp
//  LidarMCplusplus
//
//  Created by Brian Collister on 12/23/17.
//  Copyright © 2017 Brian Collister. All rights reserved.
//
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <vector>
#include <armadillo>
#include "spline.h" // https://github.com/ttk592/spline/

using namespace std;
using namespace arma;

double intersectionAngle(double x1,double y1,double z1,double x2,double y2,double z2);

int main (){
    
    //double hmm = intersectionAngle(1.123,2.3,.13,5.3,3.2,-1.1);
    
return 0;
}




double intersectionAngle(double x1,double y1,double z1,double x2,double y2,double z2){
    // intersectionAngle(c1,c2) - Calculates the intersection angle between a photon trajectory and the plane made by the lidar detector
    // In order to determine if a photon has entered the detector within the FOV of the detector, this function calculates the
    // angle between the unit vector normal to the detector plane and the propegation vector
    colvec c1; colvec c2; colvec u; colvec v;
    double angle;
    
    c1   << x1 << endr  // an array containing the first x,y,z points of the propegation vector
    << y1 << endr  // vector declarations are made using c++94 style
    << z1 << endr;
    
    c2   << x2 << endr // an array containing the first x,y,z points of the propegation vector
    << y2 << endr  // vector declarations are made using c++94 style
    << z2 << endr; // an array containing the second x,y,z points of the propegation vector
    
    u   << 0 << endr
    << 0 << endr
    << 1 << endr;   // create a unit vector normal to the plane of the detector at the origin
    
    v = c2 - c1;   // convert the propegation vector to a unit vector
    
    angle = atan2(norm(cross(u,v)),dot(u,v)); // calculate the angle between the propegation and the vector normal to the detector plane
    
    return angle;
}
