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
void BulkLinearOp<float128, float128>(float128 **elem, derivativeCol<float128>& D,Cfunction<float128, float128> *F1, grid<float128>& grid1,  float128* params, My_Int i, My_Int j, My_Int i1, My_Int j1){
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
	
#include "SourceFiles/LO.c"
}

template<>
void BoundaryLinearOp<float128, float128>(float128 **elem, derivativeCol<float128>& D,Cfunction<float128, float128> *F1, float128* params,  My_Int i, My_Int j, My_Int i1, My_Int j1){
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
	
#include "SourceFiles/BCLO.c"
}