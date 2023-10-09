//
//  main.cpp
//  Lattice
//
//  Created by Aristomenis Donos on 05/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include "omp.h"
#define _USE_MATH_DEFINES
//#define EIGEN_NO_DEBUG

#include "../function.h"
#include "../NewtonMethod.h"
#include "../ReadWrite.h"
#include <ctime>

#include "../LatticeEOMS.h"
#include "../LatticeLinearOp.h"

#include "../MyTypes.h"

#ifdef F128
typedef float128 MyType;
#elif MP
typedef mpreal MyType;
#elif LONGD
typedef long double MyType;
#else
typedef double MyType;
#endif

MyType **cheby_d1;
MyType **cheby_d2;
MyType **fourier_d1;
MyType **fourier_d2;
MyType **idr=0;
MyType **idx=0;

clock_t start1,start2;
double omp_start, omp_end;
MyType diff;


int main(int argc, const char * argv[])
{
	mpfr::mpreal::set_default_prec(76);
	
	My_Int* Rdimensions;
	My_Int* ROrders;
	My_Int* Xdimensions;
	
	My_Int rlength=GetLength<My_Int>(1, "Rdimensions.dat");
	
	Rdimensions= new My_Int[rlength];
	ROrders= new My_Int[rlength];
	Xdimensions= new My_Int[1];
	
	ReadIntParameters(Rdimensions, rlength, "Rdimensions.dat");
	ReadIntParameters(ROrders, rlength, "ROrders.dat");
	ReadIntParameters(Xdimensions, 1, "Xdimensions.dat");
	
	My_Int tnp(0);
	My_Int tnx=Xdimensions[0];
	
	for(My_Int i=0;i<rlength;++i){
		tnp+=Rdimensions[i];
	}
	
	std::cout << "Np= " << tnp << "\t, Nx= " << tnx << std::endl;
	
	My_Int XOrders[]={Xdimensions[0]};
	
	MyType* trgrid;
	MyType* txgrid;
	
	trgrid= new MyType[tnp];
	txgrid= new MyType[tnx];
	
//	setxgrid(txgrid, tnx);
	
	ReadDoubleParameters(trgrid, tnp, "rgrid.dat");
	ReadDoubleParameters(txgrid, tnx, "xgrid.dat");
	
	//Define the Chebyshev differentiation matrix and its square
	
	cheby_d1= new MyType*[tnp];
	cheby_d2= new MyType*[tnp];
	for (My_Int i1=0; i1<tnp;i1++){cheby_d1[i1]=new MyType[tnp];cheby_d2[i1]=new MyType[tnp];}
	
	//Define the Fourier differenations matrix and its square
	
	fourier_d1= new MyType*[tnx];
	fourier_d2= new MyType*[tnx];
	for (My_Int i1=0; i1<tnx ;i1++){fourier_d1[i1]=new MyType[tnx];fourier_d2[i1]=new MyType[tnx];}
	
	ReadDiffFromBin<MyType>(cheby_d1,  tnp, "cheby_d1.dat");
	ReadDiffFromBin<MyType>(cheby_d2,  tnp, "cheby_d2.dat");
	ReadDiffFromBin<MyType>(fourier_d1, tnx, "fourier_d1.dat");
	ReadDiffFromBin<MyType>(fourier_d2, tnx, "fourier_d2.dat");
	
	grid<MyType> rgrid(Rdimensions, trgrid, ROrders,rlength);
	grid<MyType> xgrid(Xdimensions, txgrid, XOrders,1);
	
	delete [] trgrid;
    delete [] txgrid;
	
	derivative<MyType> d10(cheby_d1,idx,&rgrid,&xgrid);
	derivative<MyType> d20(cheby_d2,idx,&rgrid,&xgrid);
	derivative<MyType> d01(idr,fourier_d1,&rgrid,&xgrid);
	derivative<MyType> d02(idr,fourier_d2,&rgrid,&xgrid);
	derivative<MyType> d11(cheby_d1,fourier_d1,&rgrid,&xgrid);
	
	derivative<MyType> d00(idr,idx,&rgrid,&xgrid);
	
	//Start playing
	
	MyType ps[7];
	
	ReadDoubleParameters(ps, 7, "ps.dat");
	
	MyType **test;
	MyType **test2;
	MyType **test3;
	
	test= new MyType*[tnp];
	test2= new MyType*[tnp];
	test3= new MyType*[tnp];
	
    for (My_Int i=0; i<tnp; i++) {
        test[i]=new MyType[tnx];
		test2[i]=new MyType[tnx];
		test3[i]=new MyType[tnx];
    }
	
	
	MyType lambda(ps[2]),Phi0(ps[3]),Phi1(ps[4]);
	
	MyType PI=acos(MyType(-1));
	
    for (My_Int i=0; i<tnp; i++) {
        for (My_Int j=0; j<tnx; j++) {
			test[i][j]=1.;
			test2[i][j]=1.+lambda*cos(2*PI*xgrid[j]);
			test3[i][j]=Phi0+Phi1*cos(2*PI*xgrid[j]);
        }}
	
	ReadDoubleParameters(test2[0], tnx, "chempot.dat");
	
	for(My_Int ii=1; ii<tnp; ii++){
		for(My_Int jj=0; jj<tnx ; jj++){
			test2[ii][jj]=test2[0][jj];
		}
	}
	
	Cfunction<MyType, MyType> Q[7];
	
	
    function<MyType, MyType> testf(test,&rgrid,&xgrid);
	function<MyType, MyType> testf2(test2,&rgrid,&xgrid);
	function<MyType, MyType> testf3(tnp,tnx);
	function<MyType, MyType> testf4(test3,&rgrid,&xgrid);
	
	derivativeCol<MyType> D(5);
	
	D[0]=d10; D[1]=d20; D[2]=d11; D[3]=d01; D[4]=d02;
	
	My_Int index[]={0,1,2,3,4,5,6};
	
	Q[0]=Cfunction<MyType, MyType>(testf3,D);
	Q[1]=Cfunction<MyType, MyType>(testf3,D);
	Q[2]=Cfunction<MyType, MyType>(testf3,D);
	Q[3]=Cfunction<MyType, MyType>(testf3,D);
	Q[4]=Cfunction<MyType, MyType>(testf3,D);
	Q[5]=Cfunction<MyType, MyType>(testf2,D);
	Q[6]=Cfunction<MyType, MyType>(testf4,D);
	
	ReadFromBin(Q, index, 7, tnp, tnx, "data");
	
	
	for (My_Int ii=0; ii<7; ++ii) {
		Q[ii].update(D);
	}
	
	
	{MyType w1=0, w2=0;
		w2=MaxEOM<MyType,MyType>(equations, 7, Q, rgrid,xgrid, ps);
		w1=MaxBCond<MyType,MyType>(bconditions, 7, Q, ps, xgrid, rgrid.TotalLength()-1);
		w2=w1>w2?w1:w2;
		std::cout <<"MaxError= " << w2 << std::endl;}
	
	{MyType w1=0, w2=0, w3=0, tw3=0;
		w2=MaxEOM<MyType,MyType>(equations, 7, Q, rgrid,xgrid, ps);
		w1=MaxBCond<MyType,MyType>(bconditions, 7, Q, ps, xgrid, rgrid.TotalLength()-1);
		
		w3=MaxEOM<MyType,MyType>(ksi2, 1, Q, rgrid, xgrid, ps);
		do
		{	tw3=w3;
			start1 = clock();
            omp_start = omp_get_wtime();
			UpdateFunctions<MyType,MyType>(equations, bconditions, BulkLinearOp, BoundaryLinearOp, Q, D, rgrid, xgrid, 7, 7, 7, ps );
            omp_end = omp_get_wtime();
			diff = ( std::clock() - start1) / (double)CLOCKS_PER_SEC;
			std::cout<<"CPU time: "<< diff <<'\n';
            std::cout<<"Wall time: "<< omp_end-omp_start <<'\n';
			
			w2=MaxEOM<MyType,MyType>(equations, 7, Q, rgrid,xgrid, ps);
			w1=MaxBCond<MyType,MyType>(bconditions, 7, Q, ps, xgrid, rgrid.TotalLength()-1);
			w2=w1>w2?w1:w2;
			w3=MaxEOM<MyType,MyType>(ksi2, 1, Q, rgrid, xgrid, ps);
			std::cout << "MaxError= " << w2 << std::endl;
			std::cout << "ksi= " << w3 << std::endl; }while((w2>1.E-3) || fabs(1-(tw3/w3))>1. );}

	WriteToBin(Q, index, 7, tnp, tnx, "data");
	
//	WriteStamp(t,"timestamp");
	
//	std::cout.precision(15);
	
//	for(My_Int i=0;i<7;++i){
//		std::cout << Q[i] << std::endl;
//	}
	
    return 0;
}

