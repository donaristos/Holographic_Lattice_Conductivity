//
//  TimeEvoLinearOp.cpp
//  TimeDependence
//
//  Created by Aristomenis Donos on 17/06/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "PertsEOMSLinearOp.h"

template<>
void BulkLinearOp<double, Complex>(Complex **elem, derivativeCol<double>& D,Cfunction<double, std::complex<double> > *F1, grid<double>& grid1,  double* params, My_Int i, My_Int j, My_Int i1, My_Int j1){
	Cfunction<double, std::complex<double> > &Ax=F1[0];
	Cfunction<double, std::complex<double> > &h11=F1[1];
	Cfunction<double, std::complex<double> > &Dh=F1[2];
	Cfunction<double, std::complex<double> > &htt=F1[3];
	Cfunction<double, std::complex<double> > &h22=F1[4];
	Cfunction<double, std::complex<double> > &ht1=F1[5];
	Cfunction<double, std::complex<double> > &A0=F1[6];
	Cfunction<double, std::complex<double> > &Qtt=F1[7];
	Cfunction<double, std::complex<double> > &Qrr=F1[8];
	Cfunction<double, std::complex<double> > &Q11=F1[9];
	Cfunction<double, std::complex<double> > &Qr1=F1[10];
	Cfunction<double, std::complex<double> > &Q22=F1[11];
	Cfunction<double, std::complex<double> > &a0=F1[12];
	Cfunction<double, std::complex<double> > &h=F1[13];
	typedef std::complex<double> complex;
	grid<double>& rgrid=grid1;
	double mu=params[0];
	double L=params[1]/mu;
	double w=params[7];
	double muj=1.;
	double &B1=params[6];
	double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);

#include "SourceFiles/PLO.c"

}

template<>
void BoundaryLinearOp<double, Complex>(Complex **elem, derivativeCol<double>& D,Cfunction<double, std::complex<double> > *F1, double* params,  My_Int i, My_Int j, My_Int i1, My_Int j1){
	Cfunction<double, std::complex<double> > &Ax=F1[0];
	Cfunction<double, std::complex<double> > &h11=F1[1];
	Cfunction<double, std::complex<double> > &Dh=F1[2];
	Cfunction<double, std::complex<double> > &htt=F1[3];
	Cfunction<double, std::complex<double> > &h22=F1[4];
	Cfunction<double, std::complex<double> > &ht1=F1[5];
	Cfunction<double, std::complex<double> > &A0=F1[6];
	Cfunction<double, std::complex<double> > &Qtt=F1[7];
	Cfunction<double, std::complex<double> > &Qrr=F1[8];
	Cfunction<double, std::complex<double> > &Q11=F1[9];
	Cfunction<double, std::complex<double> > &Qr1=F1[10];
	Cfunction<double, std::complex<double> > &Q22=F1[11];
	Cfunction<double, std::complex<double> > &a0=F1[12];
	Cfunction<double, std::complex<double> > &h=F1[13];
	typedef std::complex<double> complex;
	double mu=params[0];
	double L=params[1]/mu;
	double w=params[7];
	double muj=1.;
	double &B1=params[6];
	double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5), d12(6);

#include "SourceFiles/PBCLO.c"
#include "SourceFiles/PAsBCLO.c"

}

template<>
void BulkLinearOp<float128, ComplexQ>(ComplexQ **elem, derivativeCol<float128>& D,Cfunction<float128, std::complex<float128> > *F1, grid<float128>& grid1,  float128* params, My_Int i, My_Int j, My_Int i1, My_Int j1){
	Cfunction<float128, std::complex<float128> > &Ax=F1[0];
	Cfunction<float128, std::complex<float128> > &h11=F1[1];
	Cfunction<float128, std::complex<float128> > &Dh=F1[2];
	Cfunction<float128, std::complex<float128> > &htt=F1[3];
	Cfunction<float128, std::complex<float128> > &h22=F1[4];
	Cfunction<float128, std::complex<float128> > &ht1=F1[5];
	Cfunction<float128, std::complex<float128> > &A0=F1[6];
	Cfunction<float128, std::complex<float128> > &Qtt=F1[7];
	Cfunction<float128, std::complex<float128> > &Qrr=F1[8];
	Cfunction<float128, std::complex<float128> > &Q11=F1[9];
	Cfunction<float128, std::complex<float128> > &Qr1=F1[10];
	Cfunction<float128, std::complex<float128> > &Q22=F1[11];
	Cfunction<float128, std::complex<float128> > &a0=F1[12];
	Cfunction<float128, std::complex<float128> > &h=F1[13];
	typedef std::complex<float128> complex;
	grid<float128>& rgrid=grid1;
	float128 mu=params[0];
	float128 L=params[1]/mu;
	float128 w=params[7];
	float128 muj=1.q;
	float128 &B1=params[6];
	float128 &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/PLO.c"
	
}

template<>
void BoundaryLinearOp<float128, ComplexQ>(ComplexQ **elem, derivativeCol<float128>& D,Cfunction<float128, std::complex<float128> > *F1, float128* params,  My_Int i, My_Int j, My_Int i1, My_Int j1){
	Cfunction<float128, std::complex<float128> > &Ax=F1[0];
	Cfunction<float128, std::complex<float128> > &h11=F1[1];
	Cfunction<float128, std::complex<float128> > &Dh=F1[2];
	Cfunction<float128, std::complex<float128> > &htt=F1[3];
	Cfunction<float128, std::complex<float128> > &h22=F1[4];
	Cfunction<float128, std::complex<float128> > &ht1=F1[5];
	Cfunction<float128, std::complex<float128> > &A0=F1[6];
	Cfunction<float128, std::complex<float128> > &Qtt=F1[7];
	Cfunction<float128, std::complex<float128> > &Qrr=F1[8];
	Cfunction<float128, std::complex<float128> > &Q11=F1[9];
	Cfunction<float128, std::complex<float128> > &Qr1=F1[10];
	Cfunction<float128, std::complex<float128> > &Q22=F1[11];
	Cfunction<float128, std::complex<float128> > &a0=F1[12];
	Cfunction<float128, std::complex<float128> > &h=F1[13];
	typedef std::complex<float128> complex;
	float128 mu=params[0];
	float128 L=params[1]/mu;
	float128 w=params[7];
	float128 muj=1.q;
	float128 &B1=params[6];
	float128 &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5), d12(6);
	
#include "SourceFiles/PBCLO.c"
#include "SourceFiles/PAsBCLO.c"
	
}

template<>
void BulkLinearOp<mpreal, ComplexMP>(ComplexMP **elem, derivativeCol<mpreal>& D,Cfunction<mpreal, std::complex<mpreal> > *F1, grid<mpreal>& grid1,  mpreal* params, My_Int i, My_Int j, My_Int i1, My_Int j1){
	Cfunction<mpreal, std::complex<mpreal> > &Ax=F1[0];
	Cfunction<mpreal, std::complex<mpreal> > &h11=F1[1];
	Cfunction<mpreal, std::complex<mpreal> > &Dh=F1[2];
	Cfunction<mpreal, std::complex<mpreal> > &htt=F1[3];
	Cfunction<mpreal, std::complex<mpreal> > &h22=F1[4];
	Cfunction<mpreal, std::complex<mpreal> > &ht1=F1[5];
	Cfunction<mpreal, std::complex<mpreal> > &A0=F1[6];
	Cfunction<mpreal, std::complex<mpreal> > &Qtt=F1[7];
	Cfunction<mpreal, std::complex<mpreal> > &Qrr=F1[8];
	Cfunction<mpreal, std::complex<mpreal> > &Q11=F1[9];
	Cfunction<mpreal, std::complex<mpreal> > &Qr1=F1[10];
	Cfunction<mpreal, std::complex<mpreal> > &Q22=F1[11];
	Cfunction<mpreal, std::complex<mpreal> > &a0=F1[12];
	Cfunction<mpreal, std::complex<mpreal> > &h=F1[13];
	typedef std::complex<mpreal> complex;
	grid<mpreal>& rgrid=grid1;
	mpreal mu=params[0];
	mpreal L=params[1]/mu;
	mpreal w=params[7];
	mpreal muj=1;
	mpreal &B1=params[6];
	mpreal &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/PLO.c"
	
}

template<>
void BoundaryLinearOp<mpreal, ComplexMP>(ComplexMP **elem, derivativeCol<mpreal>& D,Cfunction<mpreal, std::complex<mpreal> > *F1, mpreal* params,  My_Int i, My_Int j, My_Int i1, My_Int j1){
	Cfunction<mpreal, std::complex<mpreal> > &Ax=F1[0];
	Cfunction<mpreal, std::complex<mpreal> > &h11=F1[1];
	Cfunction<mpreal, std::complex<mpreal> > &Dh=F1[2];
	Cfunction<mpreal, std::complex<mpreal> > &htt=F1[3];
	Cfunction<mpreal, std::complex<mpreal> > &h22=F1[4];
	Cfunction<mpreal, std::complex<mpreal> > &ht1=F1[5];
	Cfunction<mpreal, std::complex<mpreal> > &A0=F1[6];
	Cfunction<mpreal, std::complex<mpreal> > &Qtt=F1[7];
	Cfunction<mpreal, std::complex<mpreal> > &Qrr=F1[8];
	Cfunction<mpreal, std::complex<mpreal> > &Q11=F1[9];
	Cfunction<mpreal, std::complex<mpreal> > &Qr1=F1[10];
	Cfunction<mpreal, std::complex<mpreal> > &Q22=F1[11];
	Cfunction<mpreal, std::complex<mpreal> > &a0=F1[12];
	Cfunction<mpreal, std::complex<mpreal> > &h=F1[13];
	typedef std::complex<mpreal> complex;
	mpreal mu=params[0];
	mpreal L=params[1]/mu;
	mpreal w=params[7];
	mpreal muj=1;
	mpreal &B1=params[6];
	mpreal &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5), d12(6);
	
#include "SourceFiles/PBCLO.c"
#include "SourceFiles/PAsBCLO.c"
	
}

template<>
void BulkLinearOp<long double, ComplexLD>(ComplexLD **elem, derivativeCol<long double>& D,Cfunction<long double, std::complex<long double> > *F1, grid<long double>& grid1,  long double* params, My_Int i, My_Int j, My_Int i1, My_Int j1){
	Cfunction<long double, std::complex<long double> > &Ax=F1[0];
	Cfunction<long double, std::complex<long double> > &h11=F1[1];
	Cfunction<long double, std::complex<long double> > &Dh=F1[2];
	Cfunction<long double, std::complex<long double> > &htt=F1[3];
	Cfunction<long double, std::complex<long double> > &h22=F1[4];
	Cfunction<long double, std::complex<long double> > &ht1=F1[5];
	Cfunction<long double, std::complex<long double> > &A0=F1[6];
	Cfunction<long double, std::complex<long double> > &Qtt=F1[7];
	Cfunction<long double, std::complex<long double> > &Qrr=F1[8];
	Cfunction<long double, std::complex<long double> > &Q11=F1[9];
	Cfunction<long double, std::complex<long double> > &Qr1=F1[10];
	Cfunction<long double, std::complex<long double> > &Q22=F1[11];
	Cfunction<long double, std::complex<long double> > &a0=F1[12];
	Cfunction<long double, std::complex<long double> > &h=F1[13];
	typedef std::complex<long double> complex;
	grid<long double>& rgrid=grid1;
	long double mu=params[0];
	long double L=params[1]/mu;
	long double w=params[7];
	long double muj=1;
	long double &B1=params[6];
	long double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5);
	
#include "SourceFiles/PLO.c"
	
}

template<>
void BoundaryLinearOp<long double, ComplexLD>(ComplexLD **elem, derivativeCol<long double>& D,Cfunction<long double, std::complex<long double> > *F1, long double* params,  My_Int i, My_Int j, My_Int i1, My_Int j1){
	Cfunction<long double, std::complex<long double> > &Ax=F1[0];
	Cfunction<long double, std::complex<long double> > &h11=F1[1];
	Cfunction<long double, std::complex<long double> > &Dh=F1[2];
	Cfunction<long double, std::complex<long double> > &htt=F1[3];
	Cfunction<long double, std::complex<long double> > &h22=F1[4];
	Cfunction<long double, std::complex<long double> > &ht1=F1[5];
	Cfunction<long double, std::complex<long double> > &A0=F1[6];
	Cfunction<long double, std::complex<long double> > &Qtt=F1[7];
	Cfunction<long double, std::complex<long double> > &Qrr=F1[8];
	Cfunction<long double, std::complex<long double> > &Q11=F1[9];
	Cfunction<long double, std::complex<long double> > &Qr1=F1[10];
	Cfunction<long double, std::complex<long double> > &Q22=F1[11];
	Cfunction<long double, std::complex<long double> > &a0=F1[12];
	Cfunction<long double, std::complex<long double> > &h=F1[13];
	typedef std::complex<long double> complex;
	long double mu=params[0];
	long double L=params[1]/mu;
	long double w=params[7];
	long double muj=1;
	long double &B1=params[6];
	long double &A1=params[5];
	My_Int d10(1),d20(2),d11(3),d01(4),d02(5), d12(6);
	
#include "SourceFiles/PBCLO.c"
#include "SourceFiles/PAsBCLO.c"
	
}