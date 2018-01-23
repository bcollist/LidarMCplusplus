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
#include <complex>

using namespace std;
//using namespace arma;

// Global Variables
double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;
complex <double> imaginary = sqrt(complex<double>(-1,0));

/////////////////////// Function Prototypes ////////////////////////////

// Trig Stuuuff
double rad2Deg(double); // convert radians to degrees
double deg2Rad(double); // convert degrees to radians

// Direction Cosine Functions
double updateDirCosX(double theta, double phi, double mux, double muy, double muz); // update X direction cosine
double updateDirCosY(double theta, double phi, double mux, double muy, double muz); // update Y direction cosine
double updateDirCosZ(double theta, double phi, double mux, double muy, double muz); // update Z direction cosine

// Detector FOV Calculations
double intersectionAngle(double x1,double y1,double z1,double x2,double y2,double z2); // angle of intersection between detector and photon trajectory

double intersectionPointX(double x1,double y1,double z1,double x2,double y2,double z2); // x location of intersection with detector plane
double intersectionPointY(double x1,double y1,double z1,double x2,double y2,double z2); // y location of intersection with detector plane

// Mueller Matrix Math
arma::mat updateStokes(arma::mat stokes, arma::mat mueller, double phi, double gamma); // update photon stokes vector
double gammaCalc(double muz1, double muz2, double theta, double phi); // calculate the rotation angle

// Mie Calculations
int bhmie(double x,double refrel,int nang, double* Qscat_p, double* Qext_p, double* Qback_p, complex<double>* S1_p, complex<double>* S2_p, double* diffPoint);
double trapz(double x[], double y[], int size);

// main function
int main (){

    //////////////////////////// define constants //////////////////////////////////

    ////// Define Lidar Parameters//////

    // Detector Size and FOV //
    double detectorRad = 1.5e-1; // number of photons to trace
    //double scatLimit = 4; // number of scattering events to trace
    double FOV = deg2Rad(20); // half-angle FOV; enter in degrees -> converts to rad

    // detector position
    double xd = 0.04; double yd = 0; double zd = 0; // position of the detector in (m)
    double fd; // variable used in detector photon geometry colculations
    double anglei; // angle of intersection between photon and detector plane

    // Define water column IOPs //
    double a = 0.4; //absorption coefficient (m^-^1)
    double b = 0.5; //scattering coefficient (m^-^1)
    double c = a + b; //bema attenuation coefficient (m^-^1)
    double omega = b/c; // single scattering albedo

    // Define Mie Parameters //

    // Refractive Index
    double refMed = 1.33;
    double refPart = 1.45;//(1.3, 0.008); // relative refractive index
    double refRel = refPart/refMed;

    // Wavelength
    double lambda = 0.532; // lidar wavelength in a vaccuum (um)
    double lambdaMed = lambda/refMed; // Lidar Wavelength in Medium (um)
    double kMed = 2*pi/lambdaMed; // Lidar wavenumber in Medium;

    // Angles
    int nang = 91; // number of angles between 0-90
    int nang1 = nang*2-1; // number of angles 0-180

    // Particle Size Distribution Parameters

    // Set PSD parameters
    double Dmin = 0.1; // minimum particle diameter (um);
    double Dmax = 150.0; // maximum particle diameter (um);
    int diamBin = 100; // # of diameter bins;
    double fac = pow((Dmax/Dmin),(1.0/(diamBin-1.0))); // exponential factor necessary for defining logarithmically spaced diameter bins

    // Initialize Particle Size Arrays
    double D[diamBin]; double radius[diamBin]; double sizeParam[diamBin];
    double diffNumDistribution[diamBin];
    double* diffPoint = &diffNumDistribution[0];

    // Set PSD array values/Users/Brian/Documents/C++/LidarMCplusplus/LidarMC++.cpp

    // Diameter Array
    for (int i=0; i<diamBin; i++){ // generate an array of particle diameters
      D[i]=Dmin*pow(fac,(i)); // define the diameter bins
    }
    // Radius Array
    for (int i=0; i<diamBin; i++){ // generate an array of particle diameters
      radius[i]=D[i]/2; // define the diameter bins
    }
    // Size Parameter Array
    for (int i=0; i<diamBin; i++){ // generate an array of size parameters
      sizeParam[i] = 2*pi*radius[i]*refMed / lambda; // mie theory size parameter
    }

    // Define Mie Output Variables and Pointers
    double Qscat; double Qext; double Qback;  // Mie scattering efficiencies
    double* Qscat_p = &Qscat; double* Qext_p = &Qext; double* Qback_p = &Qback; // pointers to Mie Scattering efficiencies

    complex<double> S1[2*nang]; complex<double> S2[2*nang];
    complex<double>* S1_p = S1; complex<double>* S2_p = S2;

    // Mueller Matrix elements
    double s11[2*nang-1][diamBin];
    double s12[2*nang-1][diamBin];
    double s33[2*nang-1][diamBin];
    double s34[2*nang-1][diamBin];

    double integrandS11[2*nang-1][diamBin];
    double integrandS12[2*nang-1][diamBin];
    double integrandS33[2*nang-1][diamBin];
    double integrandS34[2*nang-1][diamBin];

    double integrandArray11[diamBin];
    double integrandArray12[diamBin];
    double integrandArray33[diamBin];
    double integrandArray34[diamBin];

    double s11bar[2*nang-1];
    double s12bar[2*nang-1];
    double s33bar[2*nang-1];
    double s34bar[2*nang-1];

    // Define Distribution //
    double k = 5E18; // differential number concentration at particle size D0

    double jungeSlope = 4.0; // slope of the junge distribution

    for (int i = 0; i<diamBin; i++){
      diffNumDistribution[i] = k*pow((D[i]/D[0]),(-1*jungeSlope)); // # of particles m^-3 um^-1
    }
/////////////////// Bulk Mie Calculations /////////////////////
    // Mie Calculations for Each Size Parameter in the distribution
    //j+1 is used to convert from fortran indexing to c++indexing
    for (int i = 0; i<diamBin; i++){
      bhmie(sizeParam[i],refRel,nang,Qscat_p, Qext_p, Qback_p, S1_p, S2_p, diffPoint);

      for (int j = 0; j<2*nang-1; j++){
        s11[j][i] = 0.5 * (pow(abs(S2[j+1]),2) + pow(abs(S1[j+1]),2));
        s12[j][i] = 0.5 * (pow(abs(S2[j+1]),2) - pow(abs(S1[j+1]),2));
        s33[j][i] = real(S1[j+1]*conj(S2[j+1]));
        s34[j][i] = imag(S2[j+1]*conj(S1[j+1]));
      }
    }

    // Define integrand to calculate bulk mueller atrix properties
    for (int i=0; i<diamBin; i++){
      for (int j=0; j<2*nang-1; j++){
        integrandS11[j][i] = diffNumDistribution[i] * s11[j][i]; // good
        integrandS12[j][i] = diffNumDistribution[i] * s12[j][i]; // good
        integrandS33[j][i] = diffNumDistribution[i] * s33[j][i]; // good
        integrandS34[j][i] = diffNumDistribution[i] * s34[j][i]; // good
    }
    }

    for (int i=0; i<2*nang-1; i++){
      for (int j=0; j<diamBin; j++){
        integrandArray11[j] = integrandS11[i][j];
        integrandArray12[j] = integrandS12[i][j];
        integrandArray33[j] = integrandS33[i][j];
        integrandArray34[j] = integrandS34[i][j];
      }
    s11bar[i] = (1.0/(kMed*kMed))*trapz(sizeParam,integrandArray11,diamBin);
    s12bar[i] = (1.0/(kMed*kMed))*trapz(sizeParam,integrandArray12,diamBin);
    s33bar[i] = (1.0/(kMed*kMed))*trapz(sizeParam,integrandArray33,diamBin);
    s34bar[i] = (1.0/(kMed*kMed))*trapz(sizeParam,integrandArray34,diamBin);
    }

    // Photon Tracking and Position Variables //
    double xT; double yT; // variable used to describe the x and y location where a photon intersects the detector plane
    double hitRad; // radial distance away from detector center that photon crosses detector plane
    double r;   // photon pathlength variable
    double theta; // off-axis scattering angle
    double phi; // scattering angle around the azimuth; if vector is coming at you, this is a counterclockwise rotation
    double gamma; // rotation angle back into photon frame of reference
    double detectAngle; // rotation angle back into the detector reference frame

    // Polarization
    arma::mat stokes; arma::mat mueller;
    arma::mat stokesDetect; arma::mat rotationDetector; //stokes vector of detected photon & rotation matrix to translate to detector reference frame


    // Signal Variables
    double dBin; // a varible describing the width of each signal bin
    double max; // maximum distance traveled by a photon
    double bd; // temporary variable use to hold the value of distance while binning the signal
    double coPol; // temporary variable used to hold the co-polarized value
    double crossPol; // temporary variable used to hold the cross-polarized value

    vector<double> binEdges; // upper edges of signal bins
    vector<double> signalWeight; // a variable used to hold the weight of each photon reaching the detector
    vector<double> distance; // distance associated with each signal bin
    vector<double> signalCOweight; // distance associated with each signal bin
    vector<double> signalCROSSweight; // distance associated with each signal bin

    // VSF Probability
    static const double thetaArray[]={0.1,0.12589,0.15849,0.19953,0.25119,0.31623,0.39811,0.50119,0.63096,0.79433,1.0,1.2589,
        1.5849,1.9953,2.5119,3.1623,3.9811,5.0119,6.3096,7.9433,10,15,20,25,30,35,40,45,50,
        55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,
        165,170,175,180};
    static const double pThetaArray[] = {0.043292475, 0.051470904, 0.061194794, 0.07278701,  0.086643078,
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

    vector<double> thetaVec (thetaArray, thetaArray + sizeof(thetaArray) / sizeof(thetaArray[0])); // vector containing angles between 0 and 180 that are used to define VSF
    vector<double> pThetaVec (pThetaArray, pThetaArray + sizeof(pThetaArray) / sizeof(pThetaArray[0])); // vector containing the cumulative probability of scattering at angles defined above

    // Initialize the spline interpolation function
    tk::spline spl; // define the spline interpolation function
    spl.set_points(pThetaVec,thetaVec); // set the x and y values the the spline interpolation will operate on

    // Mont Carlo parameters
    //nPhotons = 1000 // number of photons to trace
    //nPhotons = 10000 // number of photons to trace
    //Photons = 100000 // number of photons to trace
    //nPhotons = 1000000 // number of photons to trace
    int nPhotons = 1000000; // number of photons to trace

    // Predefined Working Variables
    mt19937::result_type seed = chrono::high_resolution_clock::now().time_since_epoch().count(); // seed the random number generator
    auto real_rand = std::bind(std::uniform_real_distribution<double>(0,1),
                             mt19937(seed));


    // Main Code
    for (int i = 0; i < nPhotons; ++i){      // loop through each individual photon

        // Photon Position and Direction Initialization
        double x1 = 0.0; double y1 = 0.0; double z1 = 0.0; // initialize photon position 1
        double x2 = 0.0; double y2 = 0.0; double z2 = 0.0; // initialize calculation positions for photons

        double mux1 = 0.0; double muy1 = 0.0; double muz1 = 1.0; // initialize new direction cosine variables
        double mux2 = 0.0; double muy2 = 0.0; double muz2 = 0.0; // initialize new direction cosine calculation variables


        // Photon Status variable
        int status = 1; // status variable 1 = alive 0 = DEAD
        double rTotal = 0; // total pathlength variable
        double nScat = 0; // number of scattering events so far
        double weight = 1; // weight of a photon
        double threshold = 0.1; // 1/10 photons will survive the roulette sequence

        // Polarization
        stokes << 1 << arma::endr     // initialize vertically polarized photon
               << -1 << arma::endr
               << 0 << arma::endr
               << 0 << arma::endr;

        while (status == 1 && nScat < 10) {   // while the photon is still alive.....



            // Move Photon
            r = -1 * log(((double) rand() / (RAND_MAX)))/c; // generate a random propegation distance
            x2 = x1 + mux1 * r; // update the photon's new x position
            y2 = y1 + muy1 * r; // update the photon's new y position
            z2 = z1 + muz1 * r; // update the photon's new z position
            // Update Pathlength Storage Variable
            rTotal = rTotal + r;

            // Did the photon cross the plane of the detector?
            if (z2 < zd){ // if the photons position is above the plane of the detector.... then it crossed the plane
                status = 0; //kill the photon

                fd = (zd - z1) / (z2 - z1); // calculate the multiplicative factor for the distance along the photon trajector to the detector
                xT = x1 + fd * (x2 - x1); // calculate x-location that photon hits plane
                yT = y1 + fd * (y2 - y1); // calculate y-location that photon hits plane
                hitRad = sqrt((xT-xd) * (xT-xd) + (yT-yd) * (yT-yd)); // distance from detector center


                // Did the photon hit the detector?
                if (hitRad < detectorRad){      // yes, if the photon hits within the radius of the detector
                    anglei = intersectionAngle(x1,y1,z1,x2,y2,z2); // calculate the angle between the detector plane and photon trajectory

                    // Did the photon hit the detector within the FOV?
                    if(anglei <= FOV){      // yes, if the intersection angle is less than the 1/2 angle FOV


                        // Create unpolarized signal
                        rTotal = rTotal - (r-(fd *r )); // calculate the distance;
                        signalWeight.push_back(weight); // append photon weight (omega^n) to signal weight vector
                        distance.push_back(rTotal); // append the total distance travelled by the photon to the distance vector
                        // Create polarized signal

                        // rotate into the detector reference frame
                        detectAngle = atan2(mux1,muy1);
                        rotationDetector  << 1 << 0 << 0 << 0 << arma::endr
                                          << 0 << (cos(2*detectAngle)) << (sin(2*detectAngle)) << 0 << arma::endr
                                          << 0 << (-1*sin(2*detectAngle)) << (cos(2*detectAngle)) << 0 << arma::endr
                                          << 0 << 0 << 0 << 1 << arma::endr;

                        stokesDetect = rotationDetector * stokes; // rotate stokes vector into the referecne frame of detector

                        coPol = (0.5 * stokesDetect[0] - 0.5 * stokesDetect[1]) * weight;
                        crossPol = (0.5 * stokesDetect[0] + 0.5 * stokesDetect[1]) * weight;

                        signalCOweight.push_back(coPol);
                        signalCROSSweight.push_back(crossPol);
                    }
                }
            }
                else{
                    theta = deg2Rad(spl(((double) rand() / (RAND_MAX)))); // generate a random scattering angle from the VSF
                    phi = 2 * pi * ((double) rand() / (RAND_MAX));  // generate a

                    mux2 = updateDirCosX(theta, phi, mux1, muy1, muz1); // update the photon X direction cosine
                    muy2 = updateDirCosY(theta, phi, mux1, muy1, muz1); // update the photon Y direction cosine
                    muz2 = updateDirCosZ(theta, phi, mux1, muy1, muz1); // update the photon Z direction cosine

                    // Update Polarization Variables
                    gamma = gammaCalc(muz1, muz2, theta, phi); // update reference frame rotation angle

                    mueller << 1 << (-sin(theta)*sin(theta)) / (1 + cos(theta) * cos(theta)) << 0 << 0 << arma::endr
                            << (-sin(theta)*sin(theta)) / (1 + cos(theta) * cos(theta)) << 1 << 0 << 0 << arma::endr
                            << 0 << 0 << 2*cos(theta) / (1+cos(theta)*cos(theta)) << 0 << arma::endr
                            << 0 << 0 << 0 << 2*cos(theta) / (1+cos(theta)*cos(theta)) << arma::endr;



//                    mueller << 1 << -1 << 0 << 0 << endr
//                            << -1 << 1 << 0 << 0 << endr
//                            << 0 << 0 << 0 << 0 << endr
//                            << 0 << 0 << 0 << 0 << endr;

                    stokes = updateStokes(stokes, mueller, phi, gamma);
                    // reset position variables
                    x1 = x2;
                    y1 = y2;
                    z1 = z2;

                    mux1 = mux2;
                    muy1 = muy2;
                    muz1 = muz2;

                    //cout << mux1*mux1 + muy1*muy1 + muz1*muz1 << endl;

                    nScat = nScat+1; // update the number of scattering events

                    // Update Photon Weight
                    nScat = nScat+1; // update the number of scattering events
                    weight = weight * omega; // update weight variable

                    // Photon Termination Roulette - allows for conservation of energy with unbiased photon termination //

                    if (weight < 0.01){ // unbiased roulette termination
                        if (rand() < threshold){
                            weight = weight * threshold;
                        }
                        else {
                            status = 0;
                        }
                    }

                }

        }
    }

    dBin = 0.25; //bin widths (m)
    max = *max_element(distance.begin(), distance.end()); //find the maximum value in the distance vector
    for (int i=0; i<ceil(max/dBin); i++){
        binEdges.push_back(i*dBin+dBin); //create a variable containing the upper bin edge of each distance bin
    }

    vector <double> signal(binEdges.size(),0.0); // initialize the final signal vector based off of the size of the distance bins vector
    vector <double> signalCO(binEdges.size(),0.0); // initialize the final signal vector based off of the size of the distance bins vector
    vector <double> signalCROSS(binEdges.size(),0.0); // initialize the final signal vector based off of the size of the distance bins vector

    for (int i=0; i<distance.size(); i++){  //loop through each element of the distance bin......
        bd = (ceil(distance[i]/dBin)*dBin); //.....find the value of the distance bin that the photon belongs to.....
        signal.at(int(bd/0.25)-1) = signal[(int(bd/0.25)-1)] + signalWeight[i]; //...add the value of the photon weight to the signal variable at the correct index for its distance bin
        signalCO.at(int(bd/0.25)-1) = signalCO[(int(bd/0.25)-1)] + signalCOweight[i]; //...add the value of the photon weight to the signal variable at the correct index for its distance bin
        signalCROSS.at(int(bd/0.25)-1) = signalCROSS[(int(bd/0.25)-1)] + signalCROSSweight[i]; //...add the value of the photon weight to the signal variable at the correct index for its distance bin
    }

    ofstream myfile;
    // Write File Header
    myfile.open ("/Users/Brian/Documents/C++/LidarMCplusplus/LidarMC.csv");
    myfile << "LidarMCplusplus.cpp output file:\n";
    myfile << "Detector Parameters:\n";
    myfile << "Radius(m) = ";
    myfile << detectorRad;
    myfile << "\n";
    myfile << "FOV(rad) = ";
    myfile << FOV;
    myfile << "\n";
    myfile << "Medium Parameters:";
    myfile << "\n";
    myfile << "a(m^-1) = ";
    myfile << a;
    myfile << "\n";
    myfile << "b(m^-1) = ";
    myfile << b;
    myfile << "\n";
    myfile << "c(m^-1) = ";
    myfile << c;
    myfile << "\n";
    myfile << "Run Parameters:\n";
    myfile << "# of photons = \n";
    myfile << nPhotons;
    myfile << "\n";
    myfile << "distance,signal,co,cross\n";
    // Write Signal
    for (int j=0; j<(signal.size()+1); j++){
        myfile << binEdges[j]/2;
        myfile << ",";
        myfile << signal[j];
        myfile << ",";
        myfile << signalCO[j];
        myfile << ",";
        myfile << signalCROSS[j];
        myfile << "\n";
    }

    myfile.close();
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
    if (abs(muz) > 0.99){
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
    if (abs(muz) > 0.99){
      muzPrime = cos(theta) * muz/abs(muz);
    }
        else{
         muzPrime = (-1*sqrt(1-muz*muz)) * sin(theta) * cos(phi) + muz*cos(theta);
       }
  return muzPrime;
}

// Determine the angle of intersection between the photon trajectory and the plane of the detector
double intersectionAngle(double x1,double y1,double z1,double x2,double y2,double z2){
// intersectionAngle(c1,c2) - Calculates the intersection angle between a photon trajectory and a vector normal to the plane of the detector
// Validated against calculations performed by hand
// In order to determine if a photon has entered the detector within the FOV of the detector, this function calculates the
// angle between the unit vector normal to the detector plane and the propegation vector
    arma::colvec c1; arma::colvec c2; arma::colvec u; arma::colvec v;
  double angle;

  c1   << x1 << arma::endr  // an array containing the first x,y,z points of the propegation vector
       << y1 << arma::endr  // vector declarations are made using c++94 style
       << z1 << arma::endr;

  c2   << x2 << arma::endr // an array containing the first x,y,z points of the propegation vector
       << y2 << arma::endr  // vector declarations are made using c++94 style
       << z2 << arma::endr; // an array containing the second x,y,z points of the propegation vector

  u   <<  0 << arma::endr
      <<  0  << arma::endr
      << -1  << arma::endr;   // create a unit vector normal to the plane of the detector at the origin

  v = c2 - c1;   // convert the propegation vector to a unit vector

  angle = atan2(norm(cross(u,v)),dot(u,v)); // calculate the angle between the propegation vector and the vector normal to the detector plane

  return angle;
}

// Determine the x location on the detector plane that the photon intersects
double intersectionPointX(double x1, double y1, double z1, double x2, double y2, double z2){
// Find the point of intersection between a line and the plane occupied by the detector (ie. the plane created by the x and y axis
// Parametric equation for a line r(t) = <x1,y1,z1> + t<x2-x1, y2-y1, z2-z1>
// or <x1 + t*(x2-x1), y1 + t(y2-y1), z1 + t(z2-z1)>
// Step (1) - plug z coordinate into plane equation and solve for t

    //Variable Initialization
    double t; double xi;

    //Function Body
    t = -1 * z1 / (z2-z1);
    xi = x1 + t*(x2-x1);

    return xi;
}

//Determine the Y location on the detector plane that the photon intersects
double intersectionPointY(double x1, double y1, double z1, double x2, double y2, double z2){
// Find the point of intersection between a line and the plane occupied by the detector (ie. the plane created by the x and y axis
// Parametric equation for a line r(t) = <x1,y1,z1> + t<x2-x1, y2-y1, z2-z1>
// or <x1 + t*(x2-x1), y1 + t(y2-y1), z1 + t(z2-z1)>
// Step (1) - plug z coordinate into plane equation and solve for t

    //Variable Initialization
    double t; double yi;

    //Function Body
    t = -1 * z1 / (z2-z1);
    yi = y1 + t*(y2-y1);

    return yi;
}

// Update the Stokes vector
arma::mat updateStokes(arma::mat stokes, arma::mat mueller, double phi, double gamma){
    // Variable Initialization
    arma::mat rotationIn; arma::mat rotationOut; arma::mat stokesPrime;

    //Function Body
    rotationIn  << 1 << 0 << 0 << 0 << arma::endr // Rotates the photon counterclockwise (if photon is coming at you) by angle phi into the scattering plane
                << 0 << (cos(-2*phi)) << (sin(-2*phi)) << 0 << arma::endr
                << 0 << (-sin(-2*phi)) << (cos(-2*phi)) << 0 << arma::endr
                << 0 << 0 << 0 << 1 << arma::endr;

    rotationOut  << 1 << 0 << 0 << 0 << arma::endr // Rotates the photon counterclockwise (if photon is coming at you) by angle gamma into the new scattering plane
                 << 0 << (cos(-2*gamma)) << (sin(-2*gamma)) << 0 << arma::endr
                 << 0 << (-sin(-2*gamma)) << (cos(-2*gamma)) << 0 << arma::endr
                 << 0 << 0 << 0 << 1 << arma::endr;

    stokesPrime = rotationOut * mueller * rotationIn * stokes;


    stokesPrime[1] = stokesPrime[1] / stokesPrime[0];
    stokesPrime[2] = stokesPrime[2] / stokesPrime[0];
    stokesPrime[3] = stokesPrime[3] / stokesPrime[0];
    stokesPrime[0] = 1.0;


    return stokesPrime;
}

// Calculate the rotation angle into the new reference frame of the photon
double gammaCalc(double muz1, double muz2, double theta, double phi){

    // Calculates the angle of rotation back into the new photon coordinate space using sperical trig cosine identity
    //cosa = cos(b)cos(c) + sin(b)sin(c)cos(A)
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

// Mie Calculations translated from Bohren and Huffman 1998

int bhmie(double x, double refrel, int nang, double* Qscat_p, double* Qext_p, double* Qback_p, complex<double>* S1_p, complex<double>* S2_p, double * diffPoint){

    // Variable Definitions
    complex<double> y; double dx; double nstop; double ymod; int nmx; double dang; double theta;
    int nn;
    double RN; double DN; double FN;
    complex<double> AN; complex<double> BN;
    double PSI; double PSI0; double PSI1;
    double CHI; double CHI0; double CHI1;
    double APSI; double APSI0; double APSI1;
    complex <double> XI; complex <double> XI0; complex <double> XI1;
    double P; double T;

    // Assign variables that will be exported to the value of their pointers
    double Qscat_f = *Qscat_p; double Qext_f = *Qext_p; double Qback_f = *Qback_p;

    // Array Definitions
    double PI[nang+1]; double PI0[nang+1]; double PI1[nang+1];
    complex<double> S1[2*nang+1]; complex<double> S2[2*nang+1];
    double AMU[nang+1]; double TAU[nang+1];

    dx = x;
    y = x * refrel;

    nstop = ceil(x + 4 * pow(x,0.3333) +2);
    ymod = abs(y);
    nmx = max(nstop,ceil(ymod))+15;
    complex<double> D[nmx];
    dang = pi/2/(nang-1);


    for (int i = 1; i<=nang; i++){
        theta = (double)(i-1) * dang;
        AMU[i] =  cos(theta);
    }

    //Logarithmic derivative D(j) calculated by downward recurence
    D[nmx]= complex <double>(0,0);
    nn = nmx-1;

    for (int n = 1; n<=nn; n++){
        RN = nmx-n+1;
        D[nmx-n]=(RN/y)-(1.0/(D[nmx-n+1]+RN/y));
    }

    for (int j = 1; j <= nang; j++){
        PI0[j] = 0.0;
        PI1[j] = 1.0;
    }

    nn = 2*nang-1;

    for (int j = 1; j <= nang; j++){
        S1[j] = 0.0;
        S2[j] = 0.0;
    }

    //Riccati-Bessel functins with real argument x calculated by upward recurence
    PSI0=cos(dx);
    PSI1=sin(dx);
    CHI0=-sin(x);
    CHI1=cos(x);
    APSI0=PSI0;
    APSI1=PSI1;
    XI0=APSI0-CHI0*imaginary;
    XI1=APSI1-CHI1*imaginary;
    Qscat_f=0.0;

    int n=1;
    while(n-1-nstop<0){
        DN=n;
        RN=n;
        FN=(2*RN+1)/(RN*(RN+1));
        PSI=(2*DN-1)*PSI1/dx-PSI0;
        APSI=PSI;
        CHI=(2*RN-1)*CHI1/x-CHI0;
        XI=APSI-CHI*imaginary;
        AN=((D[n]/refrel+RN/x)*APSI-APSI1)/((D[n]/refrel+RN/x)*XI-XI1);
        BN=((refrel*D[n]+RN/x)*APSI-APSI1)/((D[n]*refrel+RN/x)*XI-XI1);
        Qscat_f=Qscat_f+(2.0*RN+1.0)*(abs(AN)*abs(AN)+abs(BN)*abs(BN));


        for (int j=1; j<=nang; j++){
            int jj=2*nang-j;
            PI[j]=PI1[j];
            TAU[j]=RN*AMU[j]*PI[j]-(RN+1)*PI0[j];
            P=pow(-1.0,(n-1));
            S1[j]=S1[j]+FN*(AN*PI[j]+BN*TAU[j]);
            T=pow((-1),n);
            S2[j]=S2[j]+FN*(AN*TAU[j]+BN*PI[j]);
            if(j != jj){
                S1[jj]=S1[jj]+FN*(AN*PI[j]*P+BN*TAU[j]*T);
                S2[jj]=S2[jj]+FN*(BN*PI[j]*P+AN*TAU[j]*T);

            }
        }

        PSI0=PSI1;
        PSI1=PSI;
        APSI1=PSI1;
        CHI0=CHI1;
        CHI1=CHI;
        XI1=APSI1-CHI1*imaginary;
        n=n+1;
        RN=n;
        for (int j=1; j<=nang; j++){
            PI1[j]=((2*RN-1)/(RN-1))*AMU[j]*PI[j]-RN*PI0[j]/(RN-1);
            PI0[j]=PI[j];
        }
    }

    Qscat_f=(2.0/(x*x))*Qscat_f;
    Qext_f=(4.0/(x*x))*real(S1[1]);
    //Note: Qback is not Qbb, but the radar back scattering.
    Qback_f=(4.0/(x*x))*(abs(S1[2*nang-1])*abs(S1[2*nang-1]));
    //cout << "Qext_f = " << Qext_f << endl;


    // Move scattering efficiencies out of the function
    *Qscat_p = Qscat_f; // update the value of Qscat in main using a pointer
    *Qext_p  = Qext_f; // update the value of Qext in main using a pointer
    *Qback_p = Qback_f; // update the value of Qback in main using a pointer

    // Move S1 and S2 out of the function
    for (int i=1; i<=nang*2-1; i++){
      S1_p[i] = S1[i];
    }

    for (int i=1; i<=nang*2-1; i++){
      S2_p[i] = S2[i];
    }
    return 0;
}

double trapz(double x[], double y[], int size){
  double s;
  double sTemp = 0.0;
  for (int i = 0; i<size-1; i++){
    sTemp = sTemp+(x[i+1]-x[i])*(y[i+1]+y[i]);
  }
  s = 0.5*sTemp;
  return s;
}




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
