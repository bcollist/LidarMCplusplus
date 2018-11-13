#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

int main(){

  int i = 1;
  vector<string> storageS11;
  vector<string> storageS12;
  vector<string> storageS33;

  string strS11;
  string strS12;
  string strS33;

  double dbl;

  ifstream seawaterCSV("seawaterVSFZHH.csv");
  if(seawaterCSV.is_open()){
    while (seawaterCSV.good()){
        getline(seawaterCSV,strS11,',');
        getline(seawaterCSV,strS12,',');
        getline(seawaterCSV,strS33,'\n');

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

    cout << seawaterS33[i] << endl;
    
  }
}
