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

vector<double> getS11water(double lambda){
  int index;
  string temp;
  int counter;
  string s11waterstr;
  vector<double> s11water;


  ifstream s11waterCSV("seawaterVSFZHHM11.csv"); // csv file containing the mueller matrix of seawater
  if((0.380<lambda) && (lambda<0.724)){
    lambda = lambda*pow(10,3);
    cout << "getPopeandFry - wavelength input in um and converted to nm for pope and fry" << endl;
  }
  // get the index for the desired wavelength
  index = getWavelengthIndexZHH(lambda);

  // Extract desired phase function from csv matrix
  if(s11waterCSV.is_open()){
    while (s11waterCSV.good()){
      counter = 0;
      do{
        getline(s11waterCSV,temp,',');
        counter+=1;
      }while(counter<index);
      getline(s11waterCSV,s11waterstr,',');
      getline(s11waterCSV,temp,'\n');
      s11water.push_back(stod(s11waterstr));
    }
  }
  else{
    cout << "Error Opening - is the pope and fry data in your directory?" << endl;
  }
  return s11water;
}

vector<double> getS12water(double lambda){
  int index;
  string temp;
  int counter;
  string s12waterstr;
  vector<double> s12water;


  ifstream s12waterCSV("seawaterVSFZHHM12.csv"); // csv file containing the mueller matrix of seawater
  if((0.380<lambda) && (lambda<0.724)){
    lambda = lambda*pow(10,3);
    cout << "getPopeandFry - wavelength input in um and converted to nm for pope and fry" << endl;
  }
  // get the index for the desired wavelength
  index = getWavelengthIndexZHH(lambda);

  // Extract desired phase function from csv matrix
  if(s12waterCSV.is_open()){
    while (s12waterCSV.good()){
      counter = 0;
      do{
        getline(s12waterCSV,temp,',');
        counter+=1;
      }while(counter<index);
      getline(s12waterCSV,s12waterstr,',');
      getline(s12waterCSV,temp,'\n');
      s12water.push_back(stod(s12waterstr));
    }
  }
  else{
    cout << "Error Opening - is the pope and fry data in your directory?" << endl;
  }
  return s12water;
}

vector<double> getS33water(double lambda){
  int index;
  string temp;
  int counter;
  string s33waterstr;
  vector<double> s33water;


  ifstream s33waterCSV("seawaterVSFZHHM33.csv"); // csv file containing the mueller matrix of seawater
  if((0.380<lambda) && (lambda<0.724)){
    lambda = lambda*pow(10,3);
    cout << "getPopeandFry - wavelength input in um and converted to nm for pope and fry" << endl;
  }
  // get the index for the desired wavelength
  index = getWavelengthIndexZHH(lambda);

  // Extract desired phase function from csv matrix
  if(s33waterCSV.is_open()){
    while (s33waterCSV.good()){
      counter = 0;
      do{
        getline(s33waterCSV,temp,',');
        counter+=1;
      }while(counter<index);
      getline(s33waterCSV,s33waterstr,',');
      getline(s33waterCSV,temp,'\n');
      s33water.push_back(stod(s33waterstr));
    }
  }
  else{
    cout << "Error Opening - is the pope and fry data in your directory?" << endl;
  }
  return s33water;
}


int getWavelengthIndexZHH(double lambda){

  int nwavelength = 301;
  double wavelength_start = 400; // first wavelength in array
  double wavelength_int = 1; // 1nm interval between wavelengths
  int index = 0;

  index += lambda - wavelength_start*wavelength_int;
  return index;
}
