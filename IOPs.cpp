#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "spline.hpp"
#include "IOPs.hpp"

using namespace std;

// Input a wavelength 380nm-724nm, and this function returns the value of the absorption coefficient
// for pure water according to Pope and Fry (1997)
double getPopeandFry(double lambda){

  if((0.380<lambda) && (lambda<0.724)){
    lambda = lambda*pow(10,3);
    cout << "getPopeandFry - wavelength input in um and converted to nm for pope and fry" << endl;
  }

  // String storage variables
  vector<string> lambdaPF_storage;
  vector<string> aPF_storage;

  // String storage variables
  string lambdaPFstr;
  string aPFstr;

  ifstream popeandfryCSV("popeandfry.csv"); // csv file containing the mueller matrix of seawater

  int i = 1;
  if(popeandfryCSV.is_open()){
    while (popeandfryCSV.good()){
        getline(popeandfryCSV,lambdaPFstr,',');
        getline(popeandfryCSV,aPFstr,'\n');

        lambdaPF_storage.push_back(lambdaPFstr);
        aPF_storage.push_back(aPFstr);

        i++;
    }
  }
  else{
    cout << "Error Opening - is the pope and fry data in your directory?" << endl;
  }

  vector<double> lambdaPF;
  vector<double> aPF;

  for (int i = 0; i<lambdaPF_storage.size()-1; i++){
    lambdaPF.push_back(stod(lambdaPF_storage[i]));
    aPF.push_back(stod(aPF_storage[i]));
  }


  tk::spline splPF; //define the spline for the comparison function
  splPF.set_points(lambdaPF, aPF);

  double a_water = splPF(lambda);
  a_water = a_water*pow(10,2);

  return a_water;
}
