//
//  main.cpp
//  TestData
//
//  Created by Aristomenis Donos on 22/11/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#define _USE_MATH_DEFINES
//#define EIGEN_NO_DEBUG

#include <iostream>
#include <vector>
#include <limits.h>
#include <istream>
//#include "omp.h"
#include "../ReadWrite.h"
#include "../function.h"
#include "../MyTypes.h"

template<typename type>
class ccomplex : public std::complex<type>{
public:
    ccomplex(type p): std::complex<type>(p){};
    ~ccomplex();
    type re(){return (*this).real();}
    type im(){return (*this).imag();}
};

int main(int argc, const char * argv[])
{
	
    std::cout << std::numeric_limits<long>::max() << std::endl;
    std::cout << std::numeric_limits<long long>::max() << std::endl;
    std::cout << std::numeric_limits<long>::min() << std::endl;
    std::cout << std::numeric_limits<long long>::min() << std::endl;
    
    MKL_Complex16 a(9.,2.);
    std::complex<double> c(0,1);
    
    a.real=c.real();
    
    return 0;
}

