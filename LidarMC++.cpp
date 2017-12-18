#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <random>
#include <chrono>
#include <armadillo>
#include "spline.h" // https://github.com/ttk592/spline/

using namespace std;
using namespace arma;

/////////////////////// Function Prototypes ////////////////////////////

// Trig Stuuuff
double rad2Deg(double); // convert radians to degrees
double deg2Rad(double); // convert degrees to radians

// Direction Cosine Functions
double updateDirCosX(double theta, double phi, double mux, double muy, double muz);
double updateDirCosY(double theta, double phi, double mux, double muy, double muz);
double updateDirCosZ(double theta, double phi, double mux, double muy, double muz);

double randArray(int); // create an array of random numbers from


// main function
int main (){

//////////////////////////// define constants //////////////////////////////////

double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;


////// Define Lidar Parameters//////

// Detector Size and FOV
double detectorRad = 1.5e-1; // number of photons to trace
double detectorDiam = detectorRad * 2; // detector diameter (m)
double detectorArea = pi * detectorRad * detectorRad;
double scatLimit = 4; // number of scattering events to trace
double FOV = deg2Rad(20); // half-angle FOV; enter in degrees -> converts to rad

// detector position
double xd = 0.04; double yd = 0; double zd = 0; // position of the detector in (m)
double fd; // variable used in detector photon geometry colculations
double anglei; // angle of intersection between photon and detector plane

// Define water column IOPs
double a = 0.08; //absorption coefficient (m^-^1)
double b = 0.01; //scattering coefficient (m^-^1)
double c = a+b; //bema attenuation coefficient (m^-^1)
double omega = b/c; // single scattering albedo

//
double thetaArray[55] = {0.1,0.12589,0.15849,0.19953,0.25119,0.31623,0.39811,0.50119,0.63096,0.79433,1.0,1.2589,
              1.5849,1.9953,2.5119,3.1623,3.9811,5.0119,6.3096,7.9433,10,15,20,25,30,35,40,45,50,
              55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,
              165,170,175,180])

double pTheta[55] = {0.043292475, 0.051470904, 0.061194794, 0.07278701,  0.086643078,
                   0.103058065, 0.122316295, 0.144714728, 0.170373413, 0.199273394,
                   0.23138308,  0.266595059, 0.304726101, 0.345302331, 0.387470778,
                   0.431035572, 0.475832438, 0.5216971,   0.568433689, 0.615808429,
                   0.663612814, 0.749541922, 0.806027674, 0.847981298, 0.877235105,
                   0.898464649, 0.915776837, 0.929298035, 0.940418273, 0.949074367,
                   0.956340431, 0.962153282, 0.967271119, 0.971314842, 0.975042649,
                   0.978012258, 0.9808555,   0.983003728, 0.985215139, 0.986857901,
                   0.988690213, 0.990017059, 0.991533455, 0.992607569, 0.993934416,
                   0.99481898 , 0.995893094, 0.996524926, 0.997472673, 0.997978139,
                   0.998736337, 0.999052252, 0.999620901, 0.999747267, 1.0};

// Mont Carlo parameters
//nPhotons = 1000 // number of photons to trace // < 1 second
//nPhotons = 10000 // number of photons to trace // ~ 1.5 seconds
//Photons = 100000 // number of photons to trace // ~10 seconds
//nPhotons = 1000000 // number of photons to trace // ~1.5 min
nPhotons = 1000000 // number of photons to trace // ~10 min

// define the number of dpeht bins to aggregate photons into
//double depthBin [(int)((maxDepth-minDepth)/dDepth)+1] = {0};

// create depth bin variable with depth bin values
//for ( int i = 0; i<((sizeof(depthBin)/sizeof(depthBin[0])+1)); i++){
//    depthBin [i] = {(double) (i*0.25)+0.25};
//}

//int nBin = sizeof(depthBin)/sizeof(depthBin[0]); // size of bin array

//double Signal[nBin] = {0};


/////////////////// VSF Probability /////////////////////

double psi[55] = {0.1, 0.12589, 0.15849, 0.19953, 0.25119, 0.31623,
                  0.39811, 0.50119, 0.63096, 0.79433, 1.0, 1.2589, 1.5849,
                   1.9953, 2.5119, 3.1623, 3.9811, 5.0119, 6.3096, 7.9433,
                  10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90,
                  95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155,
                  160, 165, 170, 175, 180};
double pPsi[55] = {0.043292475, 0.051470904, 0.061194794, 0.07278701,  0.086643078,
                   0.103058065, 0.122316295, 0.144714728, 0.170373413, 0.199273394,
                   0.23138308,  0.266595059, 0.304726101, 0.345302331, 0.387470778,
                   0.431035572, 0.475832438, 0.5216971,   0.568433689, 0.615808429,
                   0.663612814, 0.749541922, 0.806027674, 0.847981298, 0.877235105,
                   0.898464649, 0.915776837, 0.929298035, 0.940418273, 0.949074367,
                   0.956340431, 0.962153282, 0.967271119, 0.971314842, 0.975042649,
                   0.978012258, 0.9808555,   0.983003728, 0.985215139, 0.986857901,
                   0.988690213, 0.990017059, 0.991533455, 0.992607569, 0.993934416,
                   0.99481898 , 0.995893094, 0.996524926, 0.997472673, 0.997978139,
                   0.998736337, 0.999052252, 0.999620901, 0.999747267, 1.0};

//////////////////////////////////////////////////////////////////////////////////////



// Predefined Working Variables


// Main Code

  for (int i = 0; i = nPhotons ; i++){      // loop through each individual photon

    // Photon Position and Direction Initialization
    double x1 = 0; double y1 = 0; double z1 = 0; // initialize photon position 1
    double x2 = 0; double y2 = 0; double z2 = 0; // initialize calculation positions for photons

    double mux1 = 0; double muy1 = 0, double muz1 = 1; // initialize new direction cosine variables
    double mux2 = 0; double muy2 = 0, double muz2 = 1; // initialize new direction cosine calculation variables


    // Photon Status variable
    int status = 1; // status variable 1 = alive 0 = DEAD
    double rTotal = 0; // total pathlength variable
    double weight = 1; // current weight of photon (omega^nscat)
    int nScat = 0; // number of scattering events so far
      while (status == 1 && nScat < 10) {   // while the photon is still alive.....
        // Move Photon
        double r = -1 * log(rand())/c; // generate a random propegation distance
        x2 = x1 + mux1 * r; // update the photon's new x position
        y2 = x1 + muy1 * r; // update the photon's new y position
        z2 = x1 + muz1 * r; // update the photon's new z position

        // Update Pathlength Storage Variable
        rTotal = rTotal + r;

        // Did the photon cross the plane of the detector?
        if (z2 < zd){
          fd = (zd - z1) / (z2 - z1); // calculate the multiplicative factor for the distance along the photon trajectory to the detector
          xT = x1 + fd * (x2 - x1); // calculate x-location that photon hits plane
          yT = x1 + fd * (y2 - y1); // calculate y-location that photon hits plane
          hitRad = sqrt((xT-xd) * (xT-xd) + (yT-yd) * (yT-yd)) // distance from detector center

          if (hitRad > detectorRad){
            status = 0;
            }
            else{
              anglei = pi - intersectionAngle(x1,y1,z1,x2,y2,z2); // calculate the angle betweenthe
            }
          }
        }
      }
}





return 0;
}

////////// Function Definitions ////////////

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

// Update X Direction Cosine
double updateDirCosX(double theta, double phi, double mux, double muy, double muz) {
  double muxPrime;
    if (abs(muz) > 0.99){
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
    if (abs(muz > 0.99)){
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
    if (abs(muz > 0.99)){
      muzPrime = cos(theta) * muz/abs(muz);
    }
        else{
         muzPrime = (-1*sqrt(1-muz*muz)) * sin(theta) * cos(phi) + muz*cos(theta);
       }
  return muzPrime;
}

double intersectionAngle(double x1,double y1,double z1,double x2,double y2,double z2){
// intersectionAngle(c1,c2) - Calculates the intersection angle between a photon trajectory and the plane made by the lidar detector
// In order to determine if a photon has entered the detector within the FOV of the detector, this function calculates the
// angle between the unit vector normal to the detector plane and the propegation vector
    vec c1; vec c2;


    c1 = ([x1,y1,z1]);   // an array containing the first x,y,z points of the propegation vector
    c2 = np.array([x2,y2,z2]);   // an array containing the second x,y,z points of the propegation vector
    u = np.array([0,0,1]);   // create a unit vector normal to the plane of the detector at the origin
    v = c2-c1;   // convert the propegation vector to a unit vector
    angle = atan2(LA.norm(np.cross(u,v)),np.dot(u,v)); // calculate the angle between the propegation and the vector normal to the detector plane
    return angle;
}


// double gammaCalc(double muz1, double muz2, double theta, double phi){
//   double gammaCos;
//   if (pi < phi && phi < 2 * pi){
//     gammaCos = (muz2 * cos(theta) - muz1) / sqrt((1 - cos(theta) * cos(theta)) * (1 - muz2 * muz2));
//     below = sqrt((1 - cos(theta) * cos(theta)) * (1 - muz2 * muz2));
//   }
//       else{
//         gammaCos = (muz2 * cos(theta) - muz1) / (-1 * sqrt((1 - cos(theta) * cos(theta)) * (1 - muz2 * muz2)));
//         gamma = acos(gammaCos);
//         //below = (-1 * math.sqrt((1-math.cos(theta)**2) * (1 - muz2**2)))
//       }
//     return gamma;
// }

// double gammaCalc(muz1, muz2, theta, phi){
//     double gammaCos
//     double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;
//
//     if pi < phi < 2 * pi:
//         gammaCos = (muz2*cos(theta) - muz1) / sqrt((1 - cos(theta)**2) * (1 - muz2**2))
//         //below = sqrt((1 - cos(theta)**2) * (1 - muz2**2))
//     else:
//         gammaCos = (muz2*math.cos(theta) - muz1) / (-1 * math.sqrt((1-math.cos(theta)**2) * (1 - muz2**2)))
//         //below = (-1 * math.sqrt((1-math.cos(theta)**2) * (1 - muz2**2)))
//     //print('above' + str((muz2*math.cos(theta) - muz1)))
//     //print('below'+ str(below))
//     //print(gammaCos)
//   return gamma
// }

// double intersectionAngle(x1,y1,z1,x2,y2,z2){
// # intersectionAngle(c1,c2) - Calculates the intersection angle between a photon trajectory and the plane made by the lidar detector
//     # In order to determine if a photon has entered the detector within the FOV of the detector, this function calculates the
//     # angle between the unit vector normal to the detector plane and the propegation vector
//     c1 = np.array([x1,y1,z1])   # an array containing the first x,y,z points of the propegation vector
//     c2 = np.array([x2,y2,z2])   # an array containing the second x,y,z points of the propegation vector
//     u = np.array([0,0,1])   # create a unit vector normal to the plane of the detector at the origin
//     v = c2-c1   # convert the propegation vector to a unit vector
//     angle = math.atan2(LA.norm(np.cross(u,v)),np.dot(u,v))  # calculate the angle between the propegation and the vector normal to the detector plane
//     return angle
// }


// // Generate a random number array of size (x) with values from 0-1
// double randArray(double r[], int arrayRow, int arrayCol){
//
//   // construct a random generator engine from a time-based seed:
//   unsigned seed = chrono::system_clock::now().time_since_epoch().count(); // this seeds the random number generato differently for each itteration
//
//   default_random_engine generator (seed);  // this uses the default random generator in the <random> header file
//
//   uniform_real_distribution<double> distribution(0,1); // this distributes the random numbers to a uniform distribution
//
//   double r[arrayRow][arrayCol] = {{0}};
//
//   for (int i = 0; arrayRow-1; ++i){
//         for(int j = 0; arrayCol; ++j){
//                 double r[i][j] = {distribution(generator)};  // now generate the array of size (x) from the generator and distributor
//         }
//
//   }
//
//
//   return r;
// }
//
