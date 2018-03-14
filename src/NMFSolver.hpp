#pragma once

#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
#include <armadillo>
#include <string>
#include <ctime>

#include "ChronoP.hpp"

using namespace arma;
using namespace std;

class NMFSolver
{
public:
	NMFSolver(mat &A, mat &W, mat &H, 
    int gradient_method = 0,
    int init_method = 1,
    double sparsity_coefficient = 0.001,
    int time_out_in_second = 3,
    int number_of_iteration_step = 200,
    double convergence_stop = 1e-6,
    bool verbose = false
    );

	// Utils functions
	void normalizeA();
	void normalizeW();
	void normalizeH();
	void normalizeHAndUpdateW();

	// main update
	void gradientUpdate();
	void gradientUpdateL2();
	void gradientUpdateL1();
	void gradientUpdateL21();

	// Data processing
	void reconstruction();
	double lossFunction();
	double lossFunctionL2();
	double lossFunctionL1();
	double lossFunctionL21();

	void init_from_random();
	void init_from_data();

	// Run the solver on the current instance
	void solve();

	/* data */
	mat &A;			// Matrix to factorize
	mat &W;			// W in A = W * H
	mat &H;			// H in A = W * H
	mat R;			// Reconstructed version (W*H)
	int N;			// number of sample points
	int Q;			// number of wavelength
	int K;			// number of phase in W
	int M;			// number of shifted version

    time_t init_time; // initial time
	int time_out_in_second;
	int number_of_iteration_step;
	int init_method;
	int gradient_method;
	bool verbose;
	
    double final_loss;
	const static double epsilon; // 0 threshold
	double convergence_stop;

	double sparsity_coefficient;
	double alpha; 	// gradient step
	double beta; 	// gradient step modifier
    double last_cost;

	// Data for printing for W differences(To remove)
	mat W_init; 	// initialized in solve().

};