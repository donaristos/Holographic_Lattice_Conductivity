//
//  function.cpp
//  PDEs
//
//  Created by Aristomenis Donos on 02/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "function.h"

//using namespace mpfr;

//Overload complex

std::complex<double> operator*(int i, std::complex<double> z){
	return (double)i*z;
}

std::complex<double> operator+(int i, std::complex<double> z){
	return (double)i+z;
}

std::complex<double> operator/(int i, std::complex<double> z){
	return (double)i/z;
}

std::complex<double> operator/(std::complex<double> z,int i){
	return z/(double)i;
}

std::complex<double> operator-(int i, std::complex<double> z){
	return (double)i-z;
}

std::complex<float128> operator*(double i, std::complex<float128> z){
	return (float128)i*z;
}

std::complex<float128> operator*(std::complex<float128> z,double i){
	return (float128)i*z;
}

std::complex<float128> operator+(double i, std::complex<float128> z){
	return (float128)i+z;
}

std::complex<float128> operator+(std::complex<float128> z,double i){
	return (float128)i+z;
}

std::complex<float128> operator/(double i, std::complex<float128> z){
	return (float128)i/z;
}

std::complex<float128> operator/(std::complex<float128> z,double i){
	return z/(float128)i;
}

std::complex<float128> operator-(double i, std::complex<float128> z){
	return (float128)i-z;
}

std::complex<float128> operator-(std::complex<float128> z,double i){
	return z-(float128)i;
}

std::complex<float128> pow(std::complex<float128> __x, int __n){

    if(__n >= 0){
    std::complex<float128> __y(__x);
    __y = (__n%2)?__x:std::complex<float128>(1);
    
        while (__n >>= 1)
        {
            __x = __x * __x;
            if (__n % 2)
                __y = __y * __x;
        }
        
        return __y;}
    else return 1/pow(__x, -__n);
}

std::complex<mpreal> operator*(double i, std::complex<mpreal> z){
	return (mpreal)i*z;
}

std::complex<mpreal> operator*(std::complex<mpreal> z,double i){
	return (mpreal)i*z;
}

std::complex<mpreal> operator+(double i, std::complex<mpreal> z){
	return (mpreal)i+z;
}

std::complex<mpreal> operator+(std::complex<mpreal> z,double i){
	return (mpreal)i+z;
}

std::complex<mpreal> operator/(double i, std::complex<mpreal> z){
	return (mpreal)i/z;
}

std::complex<mpreal> operator/(std::complex<mpreal> z,double i){
	return z/(mpreal)i;
}

std::complex<mpreal> operator-(double i, std::complex<mpreal> z){
	return (mpreal)i-z;
}

std::complex<mpreal> operator-(std::complex<mpreal> z,double i){
	return z-(mpreal)i;
}

std::complex<mpreal> pow(std::complex<mpreal> __x,  int __n){

        
    if(__n >= 0){
        std::complex<mpreal> __y(__x);
        __y = (__n%2)?__x:std::complex<mpreal>(1);
        
        while (__n >>= 1)
        {
            __x = __x * __x;
            if (__n % 2)
                __y = __y * __x;
        }
        
        return __y;}
    else return 1/pow(__x, -__n);
    }


std::complex<long double> operator*(double i, std::complex<long double> z){
	return (long double)i*z;
}

std::complex<long double> operator*(std::complex<long double> z,double i){
	return (long double)i*z;
}

std::complex<long double> operator+(double i, std::complex<long double> z){
	return (long double)i+z;
}

std::complex<long double> operator+(std::complex<long double> z,double i){
	return (long double)i+z;
}

std::complex<long double> operator/(double i, std::complex<long double> z){
	return (long double)i/z;
}

std::complex<long double> operator/(std::complex<long double> z,double i){
	return z/(long double)i;
}

std::complex<long double> operator-(double i, std::complex<long double> z){
	return (long double)i-z;
}

std::complex<long double> operator-(std::complex<long double> z,double i){
	return z-(long double)i;
}

std::complex<long double> operator*(int i, std::complex<long double> z){
	return (long double)i*z;
}

std::complex<long double> operator+(int i, std::complex<long double> z){
	return (long double)i+z;
}

std::complex<long double> operator/(int i, std::complex<long double> z){
	return (long double)i/z;
}

std::complex<long double> operator/(std::complex<long double> z,int i){
	return z/(long double)i;
}

std::complex<long double> operator-(int i, std::complex<long double> z){
	return (long double)i-z;
}

// Overload abs function

template<> double modulus<double, double> (double z){
	return fabs(z);
}

template<> double modulus<double, std::complex<double> > (std::complex<double> z){
	return abs(z);
}


template<> float128 modulus<float128, float128> (float128 z){
	return fabs(z);
}

template<> float128 modulus<float128, std::complex<float128> > (std::complex<float128> z){
	return abs(z);
}

template<> mpfr::mpreal modulus<mpfr::mpreal, mpfr::mpreal> (mpfr::mpreal z){
	return fabs(z);
}

template<> mpfr::mpreal modulus<mpfr::mpreal, std::complex<mpfr::mpreal> > (std::complex<mpfr::mpreal> z){
	return abs(z);
}

template<> long double modulus<long double, long double> (long double z){
	return fabsl(z);
}

template<> long double modulus<long double, std::complex<long double> > (std::complex<long double> z){
	return abs(z);
}

//Overload long double functions

long double sqrt(long double x){
	return sqrtl(x);
}

long double exp(long double x){
	return expl(x);
}

long double fabs(long double x){
	return fabsl(x);
}

long double pow(long double x, int i){
	return powl(x,i);
}

long double pow(long double x, long double y){
	return powl(x,y);
}

long double pow(long double x, double y){
	return powl(x, (long double)y);
}