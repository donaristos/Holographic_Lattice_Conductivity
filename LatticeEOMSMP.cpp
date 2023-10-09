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