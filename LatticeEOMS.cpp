//
//  EOMS.cpp
//  PDEs
//
//  Created by Aristomenis Donos on 05/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "LatticeEOMS.h"

//using namespace mpfr;

template<>
void equations<double, double>(double* elem, Cfunction<double, double> *F1, grid<double>& grid1,  double* params, My_Int i, My_Int j){
	Cfunction<double,double> &Qtt=F1[0];
	Cfunction<double,double> &Qrr=F1[1];
	Cfunction<double,double> &Q11=F1[2];
	Cfunction<double,double> &Qr1=F1[3];
	Cfunction<double,double> &Q22=F1[4];
	Cfunction<double,double> &a0=F1[5];
	Cfunction<double,double> &h=F1[6];
	grid<double>& rgrid=grid1;
	double mu=params[0];
	double L=params[1]/mu;
	double &B1=params[6];
	double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);

#include "SourceFiles/EOM.c"
}

template<>
void bconditions<double, double>(double* elem, Cfunction<double, double> *F1,double* params,  My_Int i, My_Int j){
	Cfunction<double, double> &Qtt=F1[0];
	Cfunction<double, double> &Qrr=F1[1];
	Cfunction<double, double> &Q11=F1[2];
	Cfunction<double, double> &Qr1=F1[3];
	Cfunction<double, double> &Q22=F1[4];
	Cfunction<double, double> &a0=F1[5];
	Cfunction<double, double> &h=F1[6];
	double mu=params[0];
	double L=params[1]/mu;
	double &B1=params[6];
	double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);

#include "SourceFiles/BC.c"
}

template<>
void ksi2<double, double>(double* elem, Cfunction<double, double> *F1, grid<double>& grid1, double* params, My_Int i, My_Int j){
	Cfunction<double, double> &Qtt=F1[0];
	Cfunction<double, double> &Qrr=F1[1];
	Cfunction<double, double> &Q11=F1[2];
	Cfunction<double, double> &Qr1=F1[3];
	Cfunction<double, double> &Q22=F1[4];
	Cfunction<double, double> &a0=F1[5];
	Cfunction<double, double> &h=F1[6];
	grid<double>& rgrid=grid1;
	double mu=params[0];
	double L=params[1]/mu;
	double &B1=params[6];
	double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/Ksi2.c"
}


template<>
void equations<float128, float128>(float128* elem, Cfunction<float128, float128> *F1, grid<float128>& grid1,  float128* params, My_Int i, My_Int j){
	Cfunction<float128,float128> &Qtt=F1[0];
	Cfunction<float128,float128> &Qrr=F1[1];
	Cfunction<float128,float128> &Q11=F1[2];
	Cfunction<float128,float128> &Qr1=F1[3];
	Cfunction<float128,float128> &Q22=F1[4];
	Cfunction<float128,float128> &a0=F1[5];
	Cfunction<float128,float128> &h=F1[6];
	grid<float128>& rgrid=grid1;
	float128 mu=params[0];
	float128 L=params[1]/mu;
	float128 &B1=params[6];
	float128 &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/EOM.c"
}

template<>
void bconditions<float128, float128>(float128* elem, Cfunction<float128, float128> *F1,float128* params,  My_Int i, My_Int j){
	Cfunction<float128, float128> &Qtt=F1[0];
	Cfunction<float128, float128> &Qrr=F1[1];
	Cfunction<float128, float128> &Q11=F1[2];
	Cfunction<float128, float128> &Qr1=F1[3];
	Cfunction<float128, float128> &Q22=F1[4];
	Cfunction<float128, float128> &a0=F1[5];
	Cfunction<float128, float128> &h=F1[6];
	float128 mu=params[0];
	float128 L=params[1]/mu;
	float128 &B1=params[6];
	float128 &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/BC.c"
}

template<>
void ksi2<float128, float128>(float128* elem, Cfunction<float128, float128> *F1, grid<float128>& grid1, float128* params, My_Int i, My_Int j){
	Cfunction<float128, float128> &Qtt=F1[0];
	Cfunction<float128, float128> &Qrr=F1[1];
	Cfunction<float128, float128> &Q11=F1[2];
	Cfunction<float128, float128> &Qr1=F1[3];
	Cfunction<float128, float128> &Q22=F1[4];
	Cfunction<float128, float128> &a0=F1[5];
	Cfunction<float128, float128> &h=F1[6];
	grid<float128>& rgrid=grid1;
	float128 mu=params[0];
	float128 L=params[1]/mu;
	float128 &B1=params[6];
	float128 &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/Ksi2.c"
}

template<>
void equations<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1,  mpreal* params, My_Int i, My_Int j){
	Cfunction<mpreal,mpreal> &Qtt=F1[0];
	Cfunction<mpreal,mpreal> &Qrr=F1[1];
	Cfunction<mpreal,mpreal> &Q11=F1[2];
	Cfunction<mpreal,mpreal> &Qr1=F1[3];
	Cfunction<mpreal,mpreal> &Q22=F1[4];
	Cfunction<mpreal,mpreal> &a0=F1[5];
	Cfunction<mpreal,mpreal> &h=F1[6];
	grid<mpreal>& rgrid=grid1;
	mpreal mu=params[0];
	mpreal L=params[1]/mu;
	mpreal &B1=params[6];
	mpreal &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/EOM.c"
}

template<>
void bconditions<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1,mpreal* params,  My_Int i, My_Int j){
	Cfunction<mpreal, mpreal> &Qtt=F1[0];
	Cfunction<mpreal, mpreal> &Qrr=F1[1];
	Cfunction<mpreal, mpreal> &Q11=F1[2];
	Cfunction<mpreal, mpreal> &Qr1=F1[3];
	Cfunction<mpreal, mpreal> &Q22=F1[4];
	Cfunction<mpreal, mpreal> &a0=F1[5];
	Cfunction<mpreal, mpreal> &h=F1[6];
	mpreal mu=params[0];
	mpreal L=params[1]/mu;
	mpreal &B1=params[6];
	mpreal &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/BC.c"
}
 
template<>
void ksi2<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, mpreal* params, My_Int i, My_Int j){
	Cfunction<mpreal, mpreal> &Qtt=F1[0];
	Cfunction<mpreal, mpreal> &Qrr=F1[1];
	Cfunction<mpreal, mpreal> &Q11=F1[2];
	Cfunction<mpreal, mpreal> &Qr1=F1[3];
	Cfunction<mpreal, mpreal> &Q22=F1[4];
	Cfunction<mpreal, mpreal> &a0=F1[5];
	Cfunction<mpreal, mpreal> &h=F1[6];
	grid<mpreal>& rgrid=grid1;
	mpreal mu=params[0];
	mpreal L=params[1]/mu;
	mpreal &B1=params[6];
	mpreal &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/Ksi2.c"
}

template<>
void equations<long double, long double>(long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1,  long double* params, My_Int i, My_Int j){
	Cfunction<long double,long double> &Qtt=F1[0];
	Cfunction<long double,long double> &Qrr=F1[1];
	Cfunction<long double,long double> &Q11=F1[2];
	Cfunction<long double,long double> &Qr1=F1[3];
	Cfunction<long double,long double> &Q22=F1[4];
	Cfunction<long double,long double> &a0=F1[5];
	Cfunction<long double,long double> &h=F1[6];
	grid<long double>& rgrid=grid1;
	long double mu=params[0];
	long double L=params[1]/mu;
	long double &B1=params[6];
	long double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/EOM.c"
}

template<>
void bconditions<long double, long double>(long double* elem, Cfunction<long double, long double> *F1,long double* params,  My_Int i, My_Int j){
	Cfunction<long double, long double> &Qtt=F1[0];
	Cfunction<long double, long double> &Qrr=F1[1];
	Cfunction<long double, long double> &Q11=F1[2];
	Cfunction<long double, long double> &Qr1=F1[3];
	Cfunction<long double, long double> &Q22=F1[4];
	Cfunction<long double, long double> &a0=F1[5];
	Cfunction<long double, long double> &h=F1[6];
	long double mu=params[0];
	long double L=params[1]/mu;
	long double &B1=params[6];
	long double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/BC.c"
}

template<>
void ksi2<long double, long double>(long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1, long double* params, My_Int i, My_Int j){
	Cfunction<long double, long double> &Qtt=F1[0];
	Cfunction<long double, long double> &Qrr=F1[1];
	Cfunction<long double, long double> &Q11=F1[2];
	Cfunction<long double, long double> &Qr1=F1[3];
	Cfunction<long double, long double> &Q22=F1[4];
	Cfunction<long double, long double> &a0=F1[5];
	Cfunction<long double, long double> &h=F1[6];
	grid<long double>& rgrid=grid1;
	long double mu=params[0];
	long double L=params[1]/mu;
	long double &B1=params[6];
	long double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/Ksi2.c"
}