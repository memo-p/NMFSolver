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

#include "NMFSolver.hpp"

const double NMFSolver::epsilon = 1e-10;

NMFSolver::NMFSolver(mat &A, mat &W, mat &H, 
    int gradient_method,
    int init_method,
    double sparsity_coefficient,
    int time_out_in_second,
    int number_of_iteration_step,
    double convergence_stop ,
    bool verbose
    ):
    A(A),  // A = Q * N matrix
    W(W),
    H(H),
    sparsity_coefficient(sparsity_coefficient),
    K(W.n_cols),
    Q(A.n_rows),
    N(A.n_cols),
    last_cost(0.),
    alpha(0.0001),
    beta(0.1),
    init_method(init_method),
    verbose(verbose),
    convergence_stop(convergence_stop),
    gradient_method(gradient_method),
    time_out_in_second(time_out_in_second),
    number_of_iteration_step(number_of_iteration_step)
    {
        srand(time(0));

    }


void NMFSolver::solve(){
    init_time = time(NULL);
    ChronoP chrono;
    
    // normalizeA();

    switch(init_method){
        case 0:
            init_from_data();
            break;
        case 1:
            init_from_random();
            break;
        default:
            init_from_random();
            break;
    }
    
    
    reconstruction();

    if(gradient_method == 0){ last_cost = lossFunction(); }
    else if(gradient_method == 1){last_cost = lossFunctionL1();}
    else if(gradient_method == 2){last_cost = lossFunctionL2();}
    else if(gradient_method == 3){last_cost = lossFunctionL21();}

    double current_cost = last_cost;
    chrono.Start();
    if (verbose){   printf("%f \t cost\n",current_cost); }
    

    for (int iter = 0; iter < number_of_iteration_step; ++iter) {

        if(gradient_method == 0) {       gradientUpdate();}
        else if(gradient_method ==1){    gradientUpdateL1();}
        else if(gradient_method ==2){    gradientUpdateL2();}
        else if(gradient_method ==3){    gradientUpdateL21();}
    
        reconstruction();
        if(gradient_method == 0){ current_cost = lossFunction(); }
        else if(gradient_method == 1){current_cost = lossFunctionL1();}
        else if(gradient_method == 2){current_cost = lossFunctionL2();}
        else if(gradient_method == 3){current_cost = lossFunctionL21();}
        if (verbose){   printf("%f \t cost\n",current_cost); }
        chrono.Stop();
        if( current_cost > last_cost
            || abs(last_cost - current_cost) < convergence_stop
            || chrono.ellapsed_second() > time_out_in_second){
            break;
        }
        
        last_cost = current_cost;
    }
    final_loss = last_cost;
    normalizeHAndUpdateW();
}



void NMFSolver::gradientUpdate(){
	normalizeW();
    reconstruction();
    // update H
    mat O = ones<mat>(Q, N); // an all-1 matrix
    mat Aux1 = A / (R + epsilon);
    H = H % ((W.t() * Aux1)/(W.t() * O + sparsity_coefficient));

	// update W
	reconstruction();
    if (sparsity_coefficient <= epsilon){
        Aux1 = A / (R + epsilon); 
        W = W % ((Aux1 * H.t())/(O * H.t()));
    }else{
        Aux1 = A / (R + epsilon); 
        mat Wxp = Aux1 * H.t();     // (A/R) * H
        mat Wyp = O * H.t();        // 11^t * H
        rowvec h(K);
        // numerator
        h = ones<mat>(1,Q) * (Wyp % W);     //h[k] = sum of colmun of (11^t * H) % W
        mat T1 = ones<mat>(Q, K);
        mat W_x = Wxp + (T1.each_row() % h) % W;
        // denominator
        h = ones<mat>(1,Q) * (Wxp % W);     //h[k] = sum of colmun of ((A/R) * H) % W
        mat T2 = ones<mat>(Q, K);
        mat W_y =  Wyp + (T2.each_row() % h) % W;
        mat GradW = W_x / (W_y + epsilon);
        W = W % GradW;
    }
    
}


void NMFSolver::gradientUpdateL2(){
    printf("Not Yet Implemented\n");
}



void NMFSolver::gradientUpdateL1(){
    printf("Not Yet Implemented\n");
}


void NMFSolver::gradientUpdateL21(){
    printf("Not Yet Implemented\n");
}


double NMFSolver::lossFunctionL2(){
    double l2 = norm(A-R, "fro");
    return l2*l2;
}

double NMFSolver::lossFunctionL1(){
    mat AmR = A-R;
    double l1 = 0;
    for (int q = 0; q < Q; q++)
        for (int n = 0; n < N; n++)
            l1 += abs(AmR(q,n));
    return l1;
}

double NMFSolver::lossFunctionL21(){
    mat AmR = A-R;
    AmR.transform([](double val) { return (val < 0)? val*val : val; } );
    return norm(AmR,1);
}

double NMFSolver::lossFunction(){
    double ckl = 0.0; // KL divergence(cost)
    for (int q = 0; q < Q; q++) {
        for (int m = 0; m < N; m++) {
            double a  = A(q,m);
            double r = R(q,m);
            ckl += a * log((a + epsilon)/(r + epsilon)) - a + r; 
        }
    }
    return ckl;
}






void NMFSolver::init_from_random(){
    for (int k = 0; k < K; ++k){
        for (int q = 0; q < Q; ++q){
            W(q,k) = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) + epsilon;
        }
    }
    for (int k = 0; k < K; ++k){
        for (int m = 0; m < N; ++m){
            H(k,m) = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) + epsilon;
        }
    }
    
}

void NMFSolver::init_from_data(){
    std::srand ( unsigned ( std::time(0) ) );
    std::vector<int> sampleIDs;
    for (int i=0; i<N; ++i){ sampleIDs.push_back(i); }
    std::random_shuffle ( sampleIDs.begin(), sampleIDs.end() ); 
    for (int k = 0; k < K; ++k){
        for (int q = 0; q < Q; ++q){
            W(q,k) = A(q, sampleIDs[k]);
        }
    }
    for (int k = 0; k < K; ++k){
        for (int m = 0; m < N; ++m){
            H(k,m) = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) + epsilon;
        }
    }
}





void NMFSolver::reconstruction(){
    R = W * H;
}


void NMFSolver::normalizeW(){
    for (int j = 0; j < K; j++) {
        double W2norm = norm(W.cols(j, j), "fro"); // calculate frobenius norm of the j-th column
        for (int i = 0; i < Q; i++) {
            if (W(i,j) < epsilon) continue;
            W(i,j) /= W2norm;
        }
    }
}

void NMFSolver::normalizeA(){
    A = normalise(A,1);
}

void NMFSolver::normalizeH(){
    H = normalise(H);
}

void NMFSolver::normalizeHAndUpdateW(){ 
    for (int k = 0; k < K; k++) {
        double maxH = max(H.row(k)); // get max of the line of H
        if(maxH <= epsilon){continue;}

        for (int i = 0; i < N; i++) { // normalize the line
            H(k,i) = H(k,i) / maxH;
        }
        
        for (int i = 0; i < Q; i++) { // Normalize W
            W(i,k) = W(i,k) * maxH;
        }
    }
}




