//
//  ChronoP.hpp
//  MDD
//
//  Created by guillaume perez on 04/11/2015.
//  Copyright © 2015 MemoCop. All rights reserved.
//

#pragma once

#include <iostream>
#include <chrono>
#include <ctime>

class ChronoP {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    
public:
    void Start(){  	start = std::chrono::system_clock::now(); }
	void Stop(){ 	end = std::chrono::system_clock::now();   }
    void Restart(){ start = std::chrono::system_clock::now(); }

    int64_t ellapsed_second(){
    return std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
	}
    int64_t ellapsed_m_second(){
    return std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
	}
    int64_t ellapsed_u_second(){
    return std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	}
    
    int64_t compute_granularity(){
	    std::chrono::time_point<std::chrono::system_clock> tmp;
	    std::chrono::time_point<std::chrono::system_clock> st= std::chrono::system_clock::now();
	    while ((tmp =std::chrono::system_clock::now()) == st);
	    return std::chrono::duration_cast<std::chrono::nanoseconds>(tmp-st).count();
	}
};

