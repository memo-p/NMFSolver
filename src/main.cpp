/*
 * Copyright (C) 2018 Guillaume Perez
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; If not, see <http://www.gnu.org/licenses/>.
*/
 
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
   	mat W(A.n_rows,4);
   	mat H(4,A.n_cols);
   	NMFSolver s(A,W,H,
   		0,		// KL loss
    	1, 		// Random values
    	0.001,  // Sparsity coefficient
    	2,		// Time out
    	200,	// Number of iteration
    	1e-6,	// Convergence rate
    	false,	// W fiw
   		false,	// H fix
    	false);	// Verbose
   	s.solve();
   	W.save("W.csv",csv_ascii);
   	H.save("H.csv",csv_ascii);
   	printf("Final cost    : %e   \n", s.final_loss);
   	printf("Solving time  : %d ms\n", s.chrono.ellapsed_m_second());
   	printf("#iterations   : %d   \n", s.iter_solving);

    return 0;
}

