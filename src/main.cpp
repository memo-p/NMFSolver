#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
#include <armadillo>

#include "NMFSolver.hpp"
 
using namespace std;
using namespace arma;

 
int main(int argc, char **argv){
   	
   	mat X;
   	X.load("data.csv");
   	mat A = X.t();
   	mat W(1060,4);
   	mat H(4,30);
   	NMFSolver s(A,W,H);
   	s.solve();

    return 0;
}

