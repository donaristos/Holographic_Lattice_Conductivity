//
//  InitLinearOp.cpp
//  TimeDependence
//
//  Created by Aristomenis Donos on 17/06/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "LatticeLinearOp.h"

//using namespce mpfr;


template<>
void BulkLinearOp<long double, long double>(long double **elem, derivativeCol<long double>& D,Cfunction<long double, long double> *F1, grid<long double>& grid1,  long double* params, My_Int i, My_Int j, My_Int i1, My_Int j1){
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
	
#include "SourceFiles/LO.c"
}

template<>
void BoundaryLinearOp<long double, long double>(long double **elem, derivativeCol<long double>& D,Cfunction<long double, long double> *F1, long double* params,  My_Int i, My_Int j, My_Int i1, My_Int j1){
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
	
#include "SourceFiles/BCLO.c"
}