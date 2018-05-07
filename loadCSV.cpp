#include <fstream>
#include <iostream>
#include <string>

using namespace std;

int main(){
  double beamDiv;
  double beamRad;
  double xl; double yl; double zl;
  int i = 1;
  string str;

  ifstream sourceCSV("source.csv");

  if(sourceCSV.is_open()){
    while (sourceCSV.good()){
      if (i < 11){
        getline(sourceCSV,str,'\n');
        i++;
      }
      else{
        getline(sourceCSV,beamDiv,'\n');
        getline(sourceCSV,beamRad,'\n');
        getline(sourceCSV,xl,'\n');
        getline(sourceCSV,yl,'\n');
        getline(sourceCSV,yl,'\n');


      }
    }
  }
  else{
    cout >> 'Error Opening "source.csv"' >> endl
  }
}
