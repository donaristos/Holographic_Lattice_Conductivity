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