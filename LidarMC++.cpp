#include <cstdlib> // contains the standard c++ libraries
#include <cmath> // contains the c++ math libraries
#include <iostream> // contains input/output functions (ie. cout <<)
#include <fstream> // contains function that allow file manipulation
#include <random> // contains random number generation functions
#include <chrono> // contains functions that deal with time
#include <vector> // contains vector notations
#include <armadillo> // Armadillo linear algebra
#include <complex> // allows for complex number notation
#include <string> // string stuff

#include "constants.hpp" // defines constants (pi and i)
#include "source_params.hpp" // defines source parameters
#include "bhmie.hpp" // Mie Code translated from Bohren and Huffman
#include "photon_tracking.hpp" // contains photon tracking functions
#include "IOPs.hpp"
// #include "detector_params.hpp" // defines source parameters

#include "spline.hpp" // https://github.com/ttk592/spline/
#include "erfinv.hpp" // https://gist.github.com/lakshayg/d80172fe5ae3c5d2c2aedb53c250320e


using namespace std;

/////////////////////// Function Prototypes ////////////////////////////

// Trig Stuuuff
double rad2Deg(double); // convert radians to degrees
double deg2Rad(double); // convert degrees to radians

// Detector FOV Calculations
double intersectionAngle(double x1,double y1,double z1,double x2,double y2,double z2); // angle of intersection between detector and photon trajectory
double intersectionPointX(double x1,double y1,double z1,double x2,double y2,double z2); // x location of intersection with detector plane
double intersectionPointY(double x1,double y1,double z1,double x2,double y2,double z2); // y location of intersection with detector plane

// Mueller Matrix Math
arma::mat updateStokes(arma::mat stokes, arma::mat mueller, double phi, double gamma); // update photon stokes vector

// Mie Calculations
double trapz(double x[], double y[], int size); // trapezoidal integration function

// main function
int main (){

    cout << "start" << endl;
    auto start=chrono::system_clock::now(); // start a timer

      ///////////////////////////////////
     //// Input File Parameters ////////
    ///////////////////////////////////

    string dummyLine; // use this variable to dump variable description lines
    string temp; // temporary string variable

    // Run Parameters
    string fileID;
    int runType;
    int nPhotons;

    // detection parameters
    double FOV;

    // mie parameters
    double nRe; // real refractive index
    double nIm; // imaginary refractive index

    // IOPs
    double a; // absorption coefficient (m-1)
    double b; // scattering coefficient (m-1)
    double c; // attenuation coefficient (m-1)
    double omega; // single scattering albedo

    // Particle Size Distribution
    double Dmin;
    double Dmax;
    double k;
    double jungeSlope;

    string photonFile("photon"); // first part of photon tracing filename
    string signalFile("signal"); // first part of signal trcking filename
    string muellerFile("mueller"); // first part of signal trcking filename

    ifstream lidarMCinputCSV("lidarMCinput.csv"); // open up a file stream
    if (lidarMCinputCSV.is_open()){
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,fileID,'\n'); // load fileID to name output files
      cout << fileID << endl;

      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load runType into temporary variable
      runType = stoi(temp); // store runType as an integer

      // Source Position
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load # of photons to be run into temporary variable
      xl = stod(temp); // store x-position of source as an double
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load # of photons to be run into temporary variable
      yl = stod(temp); // store y-position of source as an double
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load # of photons to be run into temporary variable
      zl = stod(temp); // store z-position  of source as an double

      // Source Polarization
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load # of photons to be run into temporary variable
      Isrc = stod(temp); // store number of photons as an integer
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load # of photons to be run into temporary variable
      Qsrc = stod(temp); // store number of photons as an integer
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load # of photons to be run into temporary variable
      Usrc = stod(temp); // store number of photons as an integer
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load # of photons to be run into temporary variable
      Vsrc = stod(temp); // store number of photons as an integer

      // Source Initial Direction
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load # of photons to be run into temporary variable
      muxsrc = stod(temp); // store number of photons as an integer
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load # of photons to be run into temporary variable
      muysrc = stod(temp); // store number of photons as an integer
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load # of photons to be run into temporary variable
      muzsrc = stod(temp); // store number of photons as an integer

      // Source Wavelength
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load # of photons to be run into temporary variable
      lambda = stod(temp); // store number of photons as an integer

      // Number of photons to trace
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load # of photons to be run into temporary variable
      nPhotons = stoi(temp); // store number of photons as an integer

      // Detector Field of View
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load FOV into temporary variable
      FOV = stod(temp); // store FOV as a double

      // Particle real refractive index
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load the real refractive index into temporary variable
      nRe = stod(temp); // store the real refractive index as a double
      // Particle imaginary refractive index
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load imaginary refractive indes into temporary variable
      nIm = stod(temp); // store the imaginary refractive index as a double

      if (runType == 1){ // skip IOPs; calculated from Mie Theory
      getline(lidarMCinputCSV,dummyLine,'\n'); // throw away variable description
      getline(lidarMCinputCSV,dummyLine,'\n'); // load file  variable
      }

      else if (runType == 2){ // load IOPs; ignore Mie Calculations of IOPs
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load file  variable
      a = stod(temp);
      getline(lidarMCinputCSV,dummyLine,','); // load file  variable
      getline(lidarMCinputCSV,temp,'\n'); // load file  variable
      b = stod(temp);
      }

      // Particle minimum diameter (um)
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load file  variable
      Dmin = stod(temp);
      // Particle maximum diameter (um)
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load file  variable
      Dmax = stod(temp);
      // Differential number concentration of particles at Dmin
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load file  variable
      k = stod(temp);
      // Junge slope of particle size distribution
      getline(lidarMCinputCSV,dummyLine,','); // throw away variable description
      getline(lidarMCinputCSV,temp,'\n'); // load file  variable
      jungeSlope = stod(temp);
    }
      ///////////////////////////////////
     ///// Define Lidar Parameters /////
    ///////////////////////////////////


    // Detector Parameters //
    double detectorRad = 1.5E-1; // detector radius(m)
    FOV=deg2Rad(FOV); // convert FOV from degrees to radians

    double xd = 0.0; double yd = 0.0; double zd = 0.0; // position of the detector in (m)
    double fd; // variable used in detector photon geometry colculations
    double anglei; // angle of intersection between photon and detector plane

    // Beam Parameters //
    // double Cwater = 0.225; // speed of light in water (m ns^-1)
    // double xl = 0.0; double yl = 0.0; double zl = 0.0; // position of the laser source in (m)
    // double pulseWidth = 7; // gaussian pulse width FWHM (ns)
    // double pulseWidth_m = pulseWidth * Cwater; // gaussian pulse width FWHM (m)
    // double beamRad = 0.5E-2; // beam radius (m)
    // double beamDiv = 2E-3; // beam divergence (radians)


      //////////////////////////////////////////
     ///// Define Water Column Parameters//////
    //////////////////////////////////////////

    // Mie Parameters //
    double refMed = 1.33; // Refractive Index of Water
    complex<double> refRel = complex<double>(nRe,nIm); // refRel stores the refractive index as a complex double
    double lambdaMed = lambda/refMed; // Lidar Wavelength in Medium (um)
    double kMed = 2*pi/(lambdaMed*1e-6); // Lidar wavenumber in Medium; convert lambda to (m)

    // Angles
    const int nang = 455; // number of angles between 0-90
    const int nangTot = (2*nang-1);
    double dang = pi/2/(nang-1);
    double angles[nangTot];
    int degreei; // index for scattering angle theta

    for (int i=0; i<nangTot; i++){
        angles[i] = (double)i*dang; // create an array of angles for mie calculations
    }

      /////////////////////////////////////////////////
     ///// Particle Size Distribution Parameters//////
    /////////////////////////////////////////////////

    // PSD bin definitions
    int diamBin = 100; // # of diameter bins;
    double fac = pow((Dmax/Dmin),(1.0/(diamBin-1.0))); // exponential factor necessary for defining logarithmically spaced diameter bins

    // Initialize Particle Size Arrays
    double radius[diamBin]; double D[diamBin]; double Dm[diamBin]; // initialize particle radius and diameters in um and m
    double sizeParam[diamBin]; // initialie mie size parameter
    double diffNumDistribution[diamBin]; //initialize array containing a differential number distribution of particles diameters

    // Diameter Array
    for (int i=0; i<diamBin; i++){ // generate an array of particle diameters
      D[i]=Dmin*pow(fac,(i)); // define the diameter bins
      Dm[i] = D[i] * 1E-6; // define a diameter bin in units of meters
    }
    // Radius Array
    for (int i=0; i<diamBin; i++){ // generate an array of particle diameters
      radius[i]=D[i]/2; // define the diameter bins
    }
    // Size Parameter Array
    for (int i=0; i<diamBin; i++){ // generate an array of size parameters
      sizeParam[i] = 2*pi*radius[i]*refMed / lambda; // mie theory size parameter
    }

    // Define Jungian Distribution //
    for (int i = 0;i<diamBin; i++){
        diffNumDistribution[i] = k*pow((D[i]/D[0]),(-1*jungeSlope)); // # of particles m^-3 um^-1
    }
    ///////////////////////////////////////////////////
    ////// Mie Output Variables and Pointers /////////
    /////////////////////////////////////////////////


    //IOP Stufffff
    double aInt[diamBin];
    double bInt[diamBin];
    double cInt[diamBin];

    // Mueller Matrix elements
    double s11[nangTot][diamBin];
    double s12[nangTot][diamBin];
    double s33[nangTot][diamBin];
    double s34[nangTot][diamBin];

    double integrandS11[nangTot][diamBin];
    double integrandS12[nangTot][diamBin];
    double integrandS33[nangTot][diamBin];
    double integrandS34[nangTot][diamBin];

    double integrandArray11[diamBin];
    double integrandArray12[diamBin];
    double integrandArray33[diamBin];
    double integrandArray34[diamBin];

    double s11bar[nangTot];
    double s12bar[nangTot];
    double s33bar[nangTot];
    double s34bar[nangTot];

    double compFunction[nangTot];
    double compFunctionI;
    double compFunctionC[nangTot];

      /////////////////////////////////////////////////
     ///// Variables loaded from other files /////////
    /////////////////////////////////////////////////

    vector<string> storageS11;
    vector<string> storageS12;
    vector<string> storageS33;

    string strS11;
    string strS12;
    string strS33;

    // Load Mueller Matrix
    ifstream seawaterMuellerCSV("seawaterVSFZHH.csv"); // csv file containing the mueller matrix of seawater
    int i = 1;
    if(seawaterMuellerCSV.is_open()){
      while (seawaterMuellerCSV.good()){
          getline(seawaterMuellerCSV,strS11,',');
          getline(seawaterMuellerCSV,strS12,',');
          getline(seawaterMuellerCSV,strS33,'\n');

          storageS11.push_back(strS11);
          storageS12.push_back(strS12);
          storageS33.push_back(strS33);

          i++;
      }
    }
    else{
      cout << "Error Opening" << endl;
    }
    double seawaterS11[storageS11.size()];
    double seawaterS12[storageS11.size()];
    double seawaterS33[storageS11.size()];

    for (int i = 0; i<storageS11.size()-1; i++){
      seawaterS11[i] = stod(storageS11[i]);
      seawaterS12[i] = stod(storageS12[i]);
      seawaterS33[i] = stod(storageS33[i]);
    }

      /////////////////////////////////////////////////
     ////////////// Mie Calculations /////////////////
    /////////////////////////////////////////////////

    // Mie Calculations for Each Size Parameter in the distribution
    //j+1 is used to convert from fortran indexing to c++indexing
    double Qscat; double Qext; double Qabs; double Qback;  // Mie scattering efficiencies
    double* Qscat_p = &Qscat; double* Qext_p = &Qext; double* Qabs_p = &Qabs; double* Qback_p = &Qback; // pointers to Mie Scattering efficiencies
    double Qb[diamBin]; double Qc[diamBin]; double Qa[diamBin]; double Qba[diamBin];
    complex<double> S1[2*nang]; complex<double> S2[2*nang];
    complex<double>* S1_p = S1; complex<double>* S2_p = S2;

    for (int i = 0; i<diamBin; i++){
      bhmie(sizeParam[i], refRel, nang, Qscat_p, Qext_p, Qabs_p, Qback_p, S1_p, S2_p);

      for (int j = 0; j<nangTot; j++){
        s11[j][i] = 0.5 * (pow(abs(S2[j+1]),2) + pow(abs(S1[j+1]),2));
        s12[j][i] = 0.5 * (pow(abs(S2[j+1]),2) - pow(abs(S1[j+1]),2));
        s33[j][i] = real(S1[j+1]*conj(S2[j+1]));
        s34[j][i] = imag(S2[j+1]*conj(S1[j+1]));
      }
      Qb[i]=Qscat_p[0]; // scattering efficiency
      Qa[i]=Qabs_p[0]; // absorption efficiency
      Qc[i]=Qext_p[0]; // extinction efficiency
      Qba[i]=Qback_p[0];
    }

    // Define integrand to calculate bulk mueller matrix properties
    for (int i=0; i<diamBin; i++){
      for (int j=0; j<nangTot; j++){
        integrandS11[j][i] = diffNumDistribution[i] * s11[j][i]; // good
        integrandS12[j][i] = diffNumDistribution[i] * s12[j][i]; // good
        integrandS33[j][i] = diffNumDistribution[i] * s33[j][i]; // good
        integrandS34[j][i] = diffNumDistribution[i] * s34[j][i]; // good
      }
    }

    for (int i=0; i<nangTot; i++){
      for (int j=0; j<diamBin; j++){
        integrandArray11[j] = integrandS11[i][j];
        integrandArray12[j] = integrandS12[i][j];
        integrandArray33[j] = integrandS33[i][j];
        integrandArray34[j] = integrandS34[i][j];
      }
    // If you integrate over particle diameter, you get the VSF
    // or you can integrate over the size parameters
    s11bar[i] = (1.0/(kMed*kMed)) * trapz(Dm,integrandArray11,diamBin);
    s11bar[i] += seawaterS11[i]; // add the contribution from seawater
    s12bar[i] = (1.0/(kMed*kMed)) * trapz(Dm,integrandArray12,diamBin);
    s12bar[i] += seawaterS12[i]; // add the contribution from seawater
    s33bar[i] = (1.0/(kMed*kMed)) * trapz(Dm,integrandArray33,diamBin);
    s33bar[i] += seawaterS33[i]; // add the contribution from seawater
    s34bar[i] = (1.0/(kMed*kMed)) * trapz(Dm,integrandArray34,diamBin);
    // add the elements of the mueller matrix of seawater here
    compFunction[i] = (s11bar[i] + abs(s12bar[i])) * sin(angles[i]);
    }

// Rejection Method From Jallion and Saint-James 2003
    compFunctionI = trapz(angles,compFunction,nangTot); // find integral of compFunction

    compFunctionC[0] = 0.0; //initialize first element of cumulative comp function

    for (int i = 1; i<nangTot; i++){ // find the cumulative integral of comp function
      compFunctionC[i] = compFunctionC[i-1]+(angles[i]-angles[i-1])*0.5*(compFunction[i]+compFunction[i-1]);
    }

    for (int i = 1; i<nangTot; i++){ // normalize the cumulative comp function to the integral
      compFunctionC[i] = compFunctionC[i]/compFunctionI;
    }

    // Calculate IOPs
    for (int i = 0; i<diamBin; i++){
      bInt[i] = Qb[i] * pi/4 *(D[i]*1E-6)*(D[i]*1E-6) * diffNumDistribution[i];
      aInt[i] = Qa[i] * pi/4 *(D[i]*1E-6)*(D[i]*1E-6) * diffNumDistribution[i];
      cInt[i] = Qc[i] * pi/4 *(D[i]*1E-6)*(D[i]*1E-6) * diffNumDistribution[i];
    }
    
    for (int i = 0; i<diamBin; i++){
      Dm[i] = D[i] * 1E-6; // diameter converted to meters
    }

    if (runType == 1){
      // IOPs - calculated from Mie Theory
      a = trapz(Dm,aInt,diamBin); // absorption coefficient(m^-1)
      double aw;
      aw = getPopeandFry(lambda);
      a += aw; // add the absorption coefficient of water (Pope and Fry)
      b = trapz(Dm,bInt,diamBin); // scattering coefficient (m^-1)
      //c = trapz(Dm,cInt,diamBin); // beam attenuation coefficient (m^-1)

    }
    c = a+b;
    omega = b/c; // single scattering albedo

    vector<double> compFunctionVec (compFunctionC, compFunctionC+sizeof(compFunctionC) / sizeof(compFunctionC[0]));
    vector<double> s11barVec (s11bar, s11bar+sizeof(s11bar) / sizeof(s11bar[0]));
    vector<double> s12barVec (s12bar, s12bar+sizeof(s12bar) / sizeof(s12bar[0]));
    vector<double> s33barVec (s33bar, s33bar+sizeof(s33bar) / sizeof(s33bar[0]));
    vector<double> s34barVec (s34bar, s34bar+sizeof(s34bar) / sizeof(s34bar[0]));
    vector<double> anglesVec (angles, angles+sizeof(angles) / sizeof(angles[0]));


    tk::spline splComp; //define the spline for the comparison function
    splComp.set_points(compFunctionVec, anglesVec); // input probability output angle

    tk::spline splS11;
    splS11.set_points(anglesVec, s11barVec); // input:angle output:S11

    tk::spline splS12;
    splS12.set_points(anglesVec, s12barVec); // input:angle output:S12

    tk::spline splS33;
    splS33.set_points(anglesVec, s12barVec); // input:angle output:S12

    tk::spline splS34;
    splS34.set_points(anglesVec, s12barVec); // input:angle output:S12


    // Photon Tracking and Position Variables //
    double xT; double yT; // variable used to describe the x and y location where a photon intersects the detector plane
    double hitRad; // radial distance away from detector center that photon crosses detector plane
    double r;   // photon pathlength variable
    double theta; // off-axis scattering angle
    double phi; // scattering angle around the azimuth; if vector is coming at you, this is a counterclockwise rotation
    double gamma; // rotation angle back into photon frame of reference
    double detectAngle; // rotation angle back into the detector reference frame

    // Polarization //
    arma::mat stokes; arma::mat mueller;
    arma::mat stokesDetect; arma::mat rotationDetector; //stokes vector of detected photon & rotation matrix to translate to detector reference frame


    // Signal Variables //
    double dBin; // a varible describing the width of each signal bin
    int nBin; // number of bins in signal output
    double max; // maximum distance traveled by a photon
    double bd; // temporary variable use to hold the value of distance while binning the signal
    double coPol; // temporary variable used to hold the co-polarized value
    double crossPol; // temporary variable used to hold the cross-polarized value

    vector<double> distance; // distance associated with each signal bin
    vector<double> signalWeight; // a variable used to hold the weight of each photon reaching the detector
    vector<double> signalCOweight; // distance associated with each signal bin
    vector<double> signalCROSSweight; // distance associated with each signal bin
    vector<double> scatHist; // a variable used to hold the number of scattering events of each photon
    vector<double> angleDet; // angle of incidence on the detector
    vector<double> binEdges; // edges of signal bins

    double I; double Irand;

//    mt19937::result_type seed = chrono::high_resolution_clock::now().time_since_epoch().count(); // seed the random number generator
//    auto real_rand = std::bind(std::uniform_real_distribution<double>(0,1),
//                             mt19937(seed));

///////////////////////////////////////////////////
/////////////// Photon Tracking  /////////////////
/////////////////////////////////////////////////

    double threshold = 0.01; // photon weight threshold to enter roulette
    double rouletteWeight = 0.1; // fraction of photons surviving roulette

    // Main Code
    for (int i = 0; i < nPhotons; ++i){      // loop through each individual photon

        // Initial Photon Position and Direction
        double x1 = 0.0; double y1 = 0.0; double z1 = 0.0; // initialize photon position 1
        double x2; double y2; double z2 ; // initialize calculation positions for photons

        double mux1 = muxsrc;
        double muy1 = muysrc;
        double muz1 = muzsrc; // initialize new direction cosine variables
        double mux2; double muy2; double muz2; // initialize new direction cosine calculation variables


        // Photon Status variable
        int status = 1; // status variable 1 = alive 0 = DEAD
        double rTotal = 0; // total pathlength variable
        int nScat = 0; // number of scattering events so far
        double weight = 1; // weight of a photon

        // Polarization
        stokes <<  Isrc  << arma::endr     // initialize vertically polarized photon
               <<  Qsrc  << arma::endr
               <<  Usrc  << arma::endr
               <<  Vsrc  << arma::endr;

        while (status == 1) {   // while the photon is still alive.....

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

                        // Add Pathlength to Distance Vector
                        rTotal -= (r-(fd *r )); // calculate the distance;
                        distance.push_back(rTotal); // append the total distance travelled by the photon to the distance vector

                        // Create unpolarized signal //
                        signalWeight.push_back(weight); // append photon weight (omega^n) to signal weight vector

                        // Create polarized signal //

                        // rotate into the detector reference frame
                        detectAngle = atan2(mux1,muy1);
                        rotationDetector  << 1 << 0 << 0 << 0 << arma::endr
                                          << 0 << (cos(2*detectAngle)) << (sin(2*detectAngle)) << 0 << arma::endr
                                          << 0 << (-1*sin(2*detectAngle)) << (cos(2*detectAngle)) << 0 << arma::endr
                                          << 0 << 0 << 0 << 1 << arma::endr;

                        stokesDetect = rotationDetector * stokes; // rotate stokes vector into the referecne frame of detector

                        // partition signalinto co- and cross- polarized using linear polarization matrices
                        // ***Warning*** currently does not account for circularly polarized light.
                        coPol = (0.5 * stokesDetect[0] - 0.5 * stokesDetect[1]) * weight;
                        crossPol = (0.5 * stokesDetect[0] + 0.5 * stokesDetect[1]) * weight;

                        signalCOweight.push_back(coPol);
                        signalCROSSweight.push_back(crossPol);

                        // Create nScat histogram
                        scatHist.push_back(nScat);

                        // Create detection angle vector
                        angleDet.push_back(anglei);
                    }
                }
            }
                else{

                  // Rejection Method - Jallion et al 2011
                    do{
                        theta = splComp((double) rand() / (RAND_MAX)); //select a random polar scatering angle based off of the comparison function
                        phi = ((double) rand() / (RAND_MAX))*2.0*pi; //select a random azimuthal scattering randomly distributed from 0-2pi
                        I = splS11(theta)+splS12(theta)*(stokes[1]*cos(2*phi)+stokes[2]*sin(2*phi))/stokes[0]; //calculate phase function using theta and phi
                        Irand = ((double) rand() / (RAND_MAX)) * splS11(theta)+abs(splS12(theta)); //rendomly draw from the comparison function
                    }while(Irand>=I); //if random draw from comparison function is >= phase function, do it again

                    mux2 = updateDirCosX(theta, phi, mux1, muy1, muz1); // update the photon X direction cosine
                    muy2 = updateDirCosY(theta, phi, mux1, muy1, muz1); // update the photon Y direction cosine
                    muz2 = updateDirCosZ(theta, phi, mux1, muy1, muz1); // update the photon Z direction cosine

                    // Update Polarization Variables //

                    // update reference frame rotation angle
                    gamma = gammaCalc(muz1, muz2, theta, phi);

                    // update mueller matrix
                    mueller << splS11(theta) << splS12(theta) << 0 << 0 << arma::endr
                            << splS12(theta) << splS11(theta) << 0 << 0 << arma::endr
                            << 0 << 0 << splS33(theta) << splS34(theta) << arma::endr
                            << 0 << 0 <<-1*splS34(theta)<< splS33(theta) << arma::endr;
                    // Update Stokes vector
                    stokes = updateStokes(stokes, mueller, phi, gamma);

                    // reset position variables
                    x1 = x2;
                    y1 = y2;
                    z1 = z2;

                    mux1 = mux2;
                    muy1 = muy2;
                    muz1 = muz2;

                    // Update Photon Weight
                    nScat += 1; // update the number of scattering events
                    weight *= omega; // update weight variable

                     //Photon Termination Roulette - allows for conservation of energy with unbiased photon termination //
                     if (weight < threshold){ // unbiased roulette termination
                         if (rand() < rouletteWeight){ // if random # is less than threshold, survives roulette
                           weight /= rouletteWeight; //
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
    nBin = ceil(max/dBin); //number of bins to store signal in
    for (int i=0; i<nBin; i++){
        binEdges.push_back(i*dBin+dBin); //create a variable containing the upper bin edge of each distance bin
    }

    vector <double> signal(binEdges.size(),0.0); // initialize the final signal vector based off of the size of the distance bins vector
    vector <double> signalCO(binEdges.size(),0.0); // initialize the final signal vector based off of the size of the distance bins vector
    vector <double> signalCROSS(binEdges.size(),0.0); // initialize the final signal vector based off of the size of the distance bins vector

    for (int i=0; i<distance.size(); i++){  //loop through each element of the distance bin......
        bd = (ceil(distance[i]/dBin)*dBin); //.....find the value of the distance bin that the photon belongs to.....
        signal.at(int(bd/dBin)-1) += signalWeight[i]; //...add the value of the photon weight to the signal variable at the correct index for its distance bin
        signalCO.at(int(bd/dBin)-1) += signalCOweight[i]; //...add the value of the photon weight to the signal variable at the correct index for its distance bin
        signalCROSS.at(int(bd/dBin)-1) += signalCROSSweight[i]; //...add the value of the photon weight to the signal variable at the correct index for its distance bin
    }

    ofstream muellerCSV;
    // Write File Header
    muellerCSV.open (muellerFile+fileID+".csv");
    muellerCSV << "LidarMCplusplus.cpp input file:\n";
    muellerCSV << "Radius(m),";
    muellerCSV << detectorRad;
    muellerCSV << "\n";
    muellerCSV << "FOV(rad),";
    muellerCSV << FOV;
    muellerCSV << "\n";
    muellerCSV << "a(m^-1),";
    muellerCSV << a;
    muellerCSV << "\n";
    muellerCSV << "b(m^-1),";
    muellerCSV << b;
    muellerCSV << "\n";
    muellerCSV << "c(m^-1),";
    muellerCSV << c;
    muellerCSV << "\n";
    muellerCSV << "Dmin,";
    muellerCSV << Dmin;
    muellerCSV << "\n";
    muellerCSV << "Dmax,";
    muellerCSV << Dmax;
    muellerCSV << "\n";
    muellerCSV << "k,";
    muellerCSV << k;
    muellerCSV << "\n";
    muellerCSV << "Junge,";
    muellerCSV << jungeSlope;
    muellerCSV << "\n";
    muellerCSV << "bulk ref index,";
    muellerCSV << refRel;
    muellerCSV << "\n";
    muellerCSV << "#photons,";
    muellerCSV << nPhotons;
    muellerCSV << "\n";
    muellerCSV << "distance,signal,co,cross\n";
    // Write Signal
    for (int j=0; j<(nangTot); j++){
        muellerCSV << angles[j];
        muellerCSV << ",";
        muellerCSV << s11bar[j];
        muellerCSV << ",";
        muellerCSV << s12bar[j];
        muellerCSV << ",";
        muellerCSV << s33bar[j];
        muellerCSV << ",";
        muellerCSV << s34bar[j];
        muellerCSV << "\n";

    }
    muellerCSV.close();

    ofstream signalCSV;
    // Write File Header
    signalCSV.open (signalFile+fileID+".csv");
    signalCSV << "LidarMCplusplus.cpp output file:\n";
    signalCSV << "Radius(m),";
    signalCSV << detectorRad;
    signalCSV << "\n";
    signalCSV << "FOV(rad),";
    signalCSV << FOV;
    signalCSV << "\n";
    signalCSV << "a(m^-1),";
    signalCSV << a;
    signalCSV << "\n";
    signalCSV << "b(m^-1),";
    signalCSV << b;
    signalCSV << "\n";
    signalCSV << "c(m^-1),";
    signalCSV << c;
    signalCSV << "\n";
    signalCSV << "Dmin,";
    signalCSV << Dmin;
    signalCSV << "\n";
    signalCSV << "Dmax,";
    signalCSV << Dmax;
    signalCSV << "\n";
    signalCSV << "k,";
    signalCSV << k;
    signalCSV << "\n";
    signalCSV << "Junge,";
    signalCSV << jungeSlope;
    signalCSV << "\n";
    signalCSV << "bulk ref index,";
    signalCSV << refRel;
    signalCSV << "\n";
    signalCSV << "#photons,";
    signalCSV << nPhotons;
    signalCSV << "\n";
    signalCSV << "distance,signal,co,cross\n";
    // Write Signal
    for (int j=0; j<(signal.size()); j++){
        signalCSV << binEdges[j];
        signalCSV << ",";
        signalCSV << signal[j];
        signalCSV << ",";
        signalCSV << signalCO[j];
        signalCSV << ",";
        signalCSV << signalCROSS[j];
        signalCSV << "\n";

    }
    signalCSV.close();




    ofstream photonCSV;
    // Write File Header
    photonCSV.open (photonFile+fileID+".csv");
    photonCSV << "LidarMCplusplus.cpp output file:\n";
    photonCSV << "Radius(m),";
    photonCSV << detectorRad;
    photonCSV << "\n";
    photonCSV << "FOV(rad),";
    photonCSV << FOV;
    photonCSV << "\n";
    photonCSV << "a(m^-1),";
    photonCSV << a;
    photonCSV << "\n";
    photonCSV << "b(m^-1),";
    photonCSV << b;
    photonCSV << "\n";
    photonCSV << "c(m^-1),";
    photonCSV << c;
    photonCSV << "\n";
    photonCSV << "Dmin,";
    photonCSV << Dmin;
    photonCSV << "\n";
    photonCSV << "Dmax,";
    photonCSV << Dmax;
    photonCSV << "\n";
    photonCSV << "k,";
    photonCSV << k;
    photonCSV << "\n";
    photonCSV << "Junge,";
    photonCSV << jungeSlope;
    photonCSV << "\n";
    photonCSV << "bulk ref index,";
    photonCSV << refRel;
    photonCSV << "\n";
    photonCSV << "#photons,";
    photonCSV << nPhotons;
    photonCSV << "\n";
    photonCSV << "distance,signal,co,cross,nScat,angleDet\n";
    // Write Signal
    for (int j=0; j<(distance.size()); j++){
        photonCSV << distance[j];
        photonCSV << ",";
        photonCSV << signalWeight[j];
        photonCSV << ",";
        photonCSV << signalCOweight[j];
        photonCSV << ",";
        photonCSV << signalCROSSweight[j];
        photonCSV << ",";
        photonCSV << scatHist[j];
        photonCSV << ",";
        photonCSV << angleDet[j];
        photonCSV << "\n";

    }
    photonCSV.close();


    auto end = chrono::system_clock::now();     // End Time

    std::chrono::duration<double> elapsed_seconds = end-start;
    cout<< "elapsed time: " << elapsed_seconds.count() << "s\n";
    return 0;
}



///////////////////////////////////////////////////
//////////// Function Definitions ////////////////
/////////////////////////////////////////////////

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

  c2   << x2 << arma::endr  // an array containing the first x,y,z points of the propegation vector
       << y2 << arma::endr  // vector declarations are made using c++94 style
       << z2 << arma::endr; // an array containing the second x,y,z points of the propegation vector

  u   <<  0  << arma::endr
      <<  0  << arma::endr
      << -1  << arma::endr;   // create a unit vector normal to the plane of the detector at the origin

  v = c2 - c1;   // convert the propagation vector to a unit vector

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
                << 0 << (cos(2*phi)) << (sin(2*phi)) << 0 << arma::endr
                << 0 << (-sin(2*phi)) << (cos(2*phi)) << 0 << arma::endr
                << 0 << 0 << 0 << 1 << arma::endr;

    rotationOut  << 1 << 0 << 0 << 0 << arma::endr // Rotates the photon counterclockwise (if photon is coming at you) by angle gamma into the new scattering plane
                 << 0 << (cos(2*gamma)) << (sin(2*gamma)) << 0 << arma::endr
                 << 0 << (-sin(2*gamma)) << (cos(2*gamma)) << 0 << arma::endr
                 << 0 << 0 << 0 << 1 << arma::endr;

    stokesPrime = rotationOut * mueller * rotationIn * stokes;


    stokesPrime[1] = stokesPrime[1] / stokesPrime[0];
    stokesPrime[2] = stokesPrime[2] / stokesPrime[0];
    stokesPrime[3] = stokesPrime[3] / stokesPrime[0];
    stokesPrime[0] = 1.0;


    return stokesPrime;
}

// Rejection method From Jallion et al. 2003


// Trapezoidal Integration Function
double trapz(double x[], double y[], int size){
  double s;
  double sTemp = 0.0;
  for (int i = 0; i<size-1; i++){
    sTemp = sTemp+(x[i+1]-x[i])*(y[i+1]+y[i]);
  }
  s = 0.5*sTemp;
  return s;
}
