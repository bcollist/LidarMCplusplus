#include <cstdlib>
#include <iostream>
#include <vector>
#include "IOPs.hpp"

using namespace std;

int main(){
  double lambda = 440;

  vector<double> x = getS11water(lambda);
  for (int i=0; i<x.size(); i++){
    cout << x[i] << endl;
  }

  vector<double> y = getS12water(lambda);
  for (int i=0; i<y.size(); i++){
    cout << y[i] << endl;
  }

  vector<double> z = getS33water(lambda);
  for (int i=0; i<z.size(); i++){
    cout << z[i] << endl;
  }
}
