//
//  NewtonMethod.h
//  PDEs
//
//  Created by Aristomenis Donos on 04/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#ifndef __PDEs__NewtonMethod__
#define __PDEs__NewtonMethod__

#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <complex>

#include "function.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <umfpack.h>

#include "mkl_pardiso.h"
#include "mkl_spblas.h"

typedef std::complex<double> Complex;



template <typename Matrix, typename Vec, typename dType,typename Type>
void ConstructLinearOp(Matrix &TriplVec, Vec &V
					   ,void (*eom) (Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, dType* params, My_Int i, My_Int j)
					   ,void (*b_cond) (Type* elem, Cfunction<dType, Type> *F1, dType* params, My_Int i, My_Int j)
					   ,void (*l_eom)(Type **elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, grid<dType>& grid1, dType* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
					   ,void (*l_b_cond)(Type **elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, dType* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<dType, Type>* f, derivativeCol<dType>& D, grid<dType>& Rgrid, grid<dType>& Xgrid ,My_Int N_functions,My_Int N_eoms, My_Int N_bconds,dType* params, const dType tolerance){
    My_Int N_unkowns=N_eoms;
	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
	My_Int nvars = (Rgrid.MaxOrder()*N_eoms+(N_bconds-N_eoms))*Xgrid.MaxOrder();
	My_Int nbdof = (Np-1)*Nx*N_eoms+(N_bconds-N_eoms)*Nx;
	My_Int res_estimate= (My_Int)nvars*(My_Int)nbdof;


#pragma omp parallel
{
	Type **telem1;
	Type **telem2;
	Type *telem3;
	Type *telem4;
	
	telem1 = new Type*[N_eoms];
	telem2 = new Type*[N_bconds];
	telem3 = new Type[N_eoms];
	telem4= new Type[N_bconds];
	
	for (My_Int ii=0; ii<N_eoms; ++ii) {
		telem1[ii]= new Type[N_unkowns];
	}
	
	for (My_Int ii=0; ii<N_bconds; ++ii) {
		telem2[ii]= new Type[N_unkowns];
	}

	Matrix TempTriplVec;
	TempTriplVec.reserve(res_estimate);

// i1 and i2 are labels of the points which we are varying
	
	#pragma omp for
	for (My_Int i1=1; i1<Np ; i1++) {
		
		My_Int lower= ((i1-Rgrid.OrderAtPoint(i1)) >= Rgrid.LeftPatchBoundary(i1))?(i1-Rgrid.OrderAtPoint(i1)):(Rgrid.LeftPatchBoundary(i1)-1);
		My_Int upper= ((i1+Rgrid.OrderAtPoint(i1)) <=  Rgrid.RightPatchBoundary(i1))?(i1+Rgrid.OrderAtPoint(i1)):(2+Rgrid.RightPatchBoundary(i1));
		
		
		lower= (lower>=1)?lower:1;
		upper= (upper<=Np-1)?upper:Np-1;

		
		for (My_Int i2=0; i2<Nx; i2++) {
			
			
			for (My_Int i=lower; i<upper; i++) {
				for (My_Int j=0; j<Nx; j++) {
					
					
					l_eom(telem1, D,f, Rgrid, params,i,j,i1, i2);
					
					for (My_Int l3=0; l3<N_unkowns; l3++) {
						for(My_Int l=0; l<N_eoms; l++){
							
													if(modulus<dType, Type>(telem1[l][l3])>tolerance)
							{
							//#pragma omp critical
							TempTriplVec.push_back(Eigen::Triplet<Type, My_Int>(j+(i-1)*Nx+l*Nx*(Np-2),i2+(i1-1)*Nx+l3*(Np-1)*Nx,telem1[l][l3]));}
							
							
						}}
				}}
			
			
			if(i1<=Rgrid.OrderAtPoint(0) || i1>Np-Rgrid.OrderAtPoint(Np-1)-1){
			for (My_Int j=0; j<Nx; j++) {
				
				
				l_b_cond(telem2, D, f, params, Np-1, j, i1, i2);
				
				for (My_Int l3=0; l3<N_unkowns; l3++) {
					for(My_Int l=0; l<N_bconds; l++){
												if(modulus<dType, Type>(telem2[l][l3])>tolerance)
						{
						//#pragma omp critical
						TempTriplVec.push_back(Eigen::Triplet<Type, My_Int>(j+l*Nx+N_unkowns*(Np-2)*Nx,i2+(i1-1)*Nx+l3*(Np-1)*Nx,telem2[l][l3]));}

					}}}
            }}}

	#pragma omp for
		for (My_Int i2=0; i2<Nx; i2++) {
			
			for (My_Int i=1; i<((Rgrid.OrderAtPoint(0)<=Rgrid.RightPatchBoundary(0))?Rgrid.OrderAtPoint(0):Rgrid.RightPatchBoundary(0)); i++) {
				for (My_Int j=0; j<Nx; j++) {
					
					
					l_eom(telem1, D,f, Rgrid, params,i,j,0, i2);
					
					for (My_Int l3=0; l3<N_bconds-N_eoms; l3++) {
						for(My_Int l=0; l<N_eoms; l++){
						
													if(modulus<dType, Type>(telem1[l][l3])>tolerance)
							{
							//#pragma omp critical
							TempTriplVec.push_back(Eigen::Triplet<Type, My_Int>(j+(i-1)*Nx+l*Nx*(Np-2),N_eoms*(Np-1)*Nx+l3*Nx+i2,telem1[l][l3]));}
							
							
						}}
				}}
			
			
			for (My_Int j=0; j<Nx; j++) {
				
				
				l_b_cond(telem2, D, f, params, Np-1, j, 0, i2);
				
				for (My_Int l3=0; l3<N_bconds-N_eoms; l3++) {
					for(My_Int l=0; l<N_bconds; l++){
						
												if(modulus<dType, Type>(telem2[l][l3])>tolerance)
						{
						//#pragma omp critical
						TempTriplVec.push_back(Eigen::Triplet<Type, My_Int>(j+l*Nx+N_unkowns*(Np-2)*Nx,N_eoms*(Np-1)*Nx+l3*Nx+i2,telem2[l][l3]));}

					}}}
		}
	
    
	#pragma omp for
	for (My_Int i=1; i<Np-1; i++) {
		for (My_Int j=0; j<Nx; j++) {
			eom(telem3,f, Rgrid, params, i,  j);
			for(My_Int l=0; l<N_eoms; l++){
				V[j+(i-1)*Nx+l*Nx*(Np-2)]=telem3[l];
			}}}
	
	#pragma omp for
	for (My_Int j=0; j<Nx; j++) {
		b_cond(telem4, f, params, Np-1,  j);
		for(My_Int l=0; l<N_bconds; l++){
			V[j+l*Nx+N_unkowns*(Np-2)*Nx]=telem4[l];
		}}
    

	for (My_Int ii=0; ii<N_eoms; ++ii) {
		delete[] telem1[ii];
	}
	
	for (My_Int ii=0; ii<N_bconds; ++ii) {
		delete[] telem2[ii];
	}
	delete[] telem1;
	delete[] telem2;
	delete[] telem3;
	delete[] telem4;
    
	#pragma omp critical
	TriplVec.insert(TriplVec.end(),TempTriplVec.begin(),TempTriplVec.end());
	}
}


template<typename dType, typename Type>
void UpdateFunctions( void (*eom) (Type* elem, Cfunction<dType,Type> *F1, grid<dType>& grid, dType* params, My_Int i, My_Int j)
					 ,void (*b_cond) (Type* elem, Cfunction<dType,Type> *F1,dType* params, My_Int i, My_Int j)
					 ,void (*l_eom)(Type **elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, grid<dType>& grid, dType* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
					 ,void (*l_b_cond)(Type **elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, dType* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<dType, Type>* f, derivativeCol<dType>& D, grid<dType>& Rgrid, grid<dType>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds,dType* params){ }

template<>
void UpdateFunctions<double, double>( void (*eom) (double* elem, Cfunction<double, double> *F1, grid<double>& grid, double* params, My_Int i, My_Int j)
					 ,void (*b_cond) (double* elem, Cfunction<double, double> *F1,double* params, My_Int i, My_Int j)
					 ,void (*l_eom)(double **elem, derivativeCol<double>& Der,Cfunction<double, double> *F1, grid<double>& grid, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
							 ,void (*l_b_cond)(double **elem,  derivativeCol<double>& Der,Cfunction<double, double> *F1, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<double, double>* f, derivativeCol<double>& D, grid<double>& Rgrid,grid<double>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds, double* params);

template<>
void UpdateFunctions<double, Complex>( void (*eom) (Complex* elem, Cfunction<double, Complex> *F1 , grid<double>& grid, double* params, My_Int i, My_Int j)
					 ,void (*b_cond) (Complex* elem, Cfunction<double, Complex> *F1,double* params, My_Int i, My_Int j)
					 ,void (*l_eom)(Complex **elem,  derivativeCol<double>& Der,Cfunction<double, Complex> *F1, grid<double>& grid, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
							  ,void (*l_b_cond)(Complex **elem,  derivativeCol<double>& Der,Cfunction<double, Complex> *F1, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<double, Complex>* f, derivativeCol<double>& D, grid<double>& Rgrid, grid<double>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds,double* params);


template<>
void UpdateFunctions<float128, float128>( void (*eom) (float128* elem, Cfunction<float128, float128> *F1, grid<float128>& grid, float128* params, My_Int i, My_Int j)
										 ,void (*b_cond) (float128* elem, Cfunction<float128, float128> *F1,float128* params, My_Int i, My_Int j)
										 ,void (*l_eom)(float128 **elem, derivativeCol<float128>& Der,Cfunction<float128, float128> *F1, grid<float128>& grid, float128* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
										 ,void (*l_b_cond)(float128 **elem,  derivativeCol<float128>& Der,Cfunction<float128, float128> *F1, float128* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<float128, float128>* f, derivativeCol<float128>& D, grid<float128>& Rgrid,grid<float128>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds, float128* params);

template<>
void UpdateFunctions<float128, ComplexQ>( void (*eom) (ComplexQ* elem, Cfunction<float128, ComplexQ> *F1 , grid<float128>& grid, float128* params, My_Int i, My_Int j)
										 ,void (*b_cond) (ComplexQ* elem, Cfunction<float128, ComplexQ> *F1,float128* params, My_Int i, My_Int j)
										 ,void (*l_eom)(ComplexQ **elem,  derivativeCol<float128>& Der,Cfunction<float128, ComplexQ> *F1, grid<float128>& grid, float128* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
										 ,void (*l_b_cond)(ComplexQ **elem,  derivativeCol<float128>& Der,Cfunction<float128, ComplexQ> *F1, float128* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<float128, ComplexQ>* f, derivativeCol<float128>& D, grid<float128>& Rgrid, grid<float128>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds,float128* params);


template<>
void UpdateFunctions<mpreal, mpreal>( void (*eom) (mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid, mpreal* params, My_Int i, My_Int j)
									 ,void (*b_cond) (mpreal* elem, Cfunction<mpreal, mpreal> *F1,mpreal* params, My_Int i, My_Int j)
									 ,void (*l_eom)(mpreal **elem, derivativeCol<mpreal>& Der,Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid, mpreal* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
									 ,void (*l_b_cond)(mpreal **elem,  derivativeCol<mpreal>& Der,Cfunction<mpreal, mpreal> *F1, mpreal* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<mpreal, mpreal>* f, derivativeCol<mpreal>& D, grid<mpreal>& Rgrid,grid<mpreal>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds, mpreal* params);

template<>
void UpdateFunctions<mpreal, ComplexMP>( void (*eom) (ComplexMP* elem, Cfunction<mpreal, ComplexMP> *F1 , grid<mpreal>& grid, mpreal* params, My_Int i, My_Int j)
										,void (*b_cond) (ComplexMP* elem, Cfunction<mpreal, ComplexMP> *F1,mpreal* params, My_Int i, My_Int j)
										,void (*l_eom)(ComplexMP **elem,  derivativeCol<mpreal>& Der,Cfunction<mpreal, ComplexMP> *F1, grid<mpreal>& grid, mpreal* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
										,void (*l_b_cond)(ComplexMP **elem,  derivativeCol<mpreal>& Der,Cfunction<mpreal, ComplexMP> *F1, mpreal* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<mpreal, ComplexMP>* f, derivativeCol<mpreal>& D, grid<mpreal>& Rgrid, grid<mpreal>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds,mpreal* params);

template<>
void UpdateFunctions<long double, long double>( void (*eom) (long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid, long double* params, My_Int i, My_Int j)
											   ,void (*b_cond) (long double* elem, Cfunction<long double, long double> *F1,long double* params, My_Int i, My_Int j)
											   ,void (*l_eom)(long double **elem, derivativeCol<long double>& Der,Cfunction<long double, long double> *F1, grid<long double>& grid, long double* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
											   ,void (*l_b_cond)(long double **elem,  derivativeCol<long double>& Der,Cfunction<long double, long double> *F1, long double* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<long double, long double>* f, derivativeCol<long double>& D, grid<long double>& Rgrid,grid<long double>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds, long double* params);

template<>
void UpdateFunctions<long double, ComplexLD>( void (*eom) (ComplexLD* elem, Cfunction<long double, ComplexLD> *F1 , grid<long double>& grid, long double* params, My_Int i, My_Int j)
											 ,void (*b_cond) (ComplexLD* elem, Cfunction<long double, ComplexLD> *F1,long double* params, My_Int i, My_Int j)
											 ,void (*l_eom)(ComplexLD **elem,  derivativeCol<long double>& Der,Cfunction<long double, ComplexLD> *F1, grid<long double>& grid, long double* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
											 ,void (*l_b_cond)(ComplexLD **elem,  derivativeCol<long double>& Der,Cfunction<long double, ComplexLD> *F1, long double* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<long double, ComplexLD>* f, derivativeCol<long double>& D, grid<long double>& Rgrid, grid<long double>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds,long double* params);

template<typename dType, typename Type>
dType MaxEOM(void (*eom) (Type* elem, Cfunction<dType, Type> *F1, grid<dType>& Rgrid, dType* params, My_Int i, My_Int j),My_Int N_eoms,Cfunction<dType, Type> *F1, grid<dType>& Rgrid, grid<dType>& Xgrid, dType* params){
	dType w2=0;
	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
	
	#pragma omp parallel
 	{Type *telem;
 	telem= new Type[N_eoms];
 	dType tw1=0, tw2=0;
	
 	#pragma omp for
 	for (My_Int i1=1; i1<Np-1; i1++) {
 		for (My_Int i2=0; i2<Nx; i2++) {
 			eom(telem,F1, Rgrid, params, i1, i2);
 			for (My_Int l=0;l<N_eoms; l++){
 				tw1=modulus<dType, Type>(telem[l]);
 				tw2=tw1>tw2?tw1:tw2;
 			}}}

	delete[] telem;
	
	#pragma omp critical
	{w2= tw2>w2?tw2:w2;}
	}
	
	return w2;
}

template<typename dType, typename Type>
dType MaxBCond(void (*b_cond) (Type* elem, Cfunction<dType, Type> *F1,dType* params, My_Int i, My_Int j), My_Int N_bconds,Cfunction<dType, Type> *F1,dType* params, grid<dType>& Xgrid, My_Int position){
	dType w2=0;
	My_Int Nx=Xgrid.TotalLength();
	
	#pragma omp parallel
	{Type *telem;
	telem= new Type[N_bconds];
	dType tw1=0, tw2=0;
	
	#pragma omp for
	for (My_Int i2=0; i2<Nx; i2++) {
		b_cond(telem, F1, params, position, i2);
		for (My_Int l=0;l<N_bconds; l++){
			tw1= modulus<dType, Type>(telem[l]);
			tw2=tw1>tw2?tw1:tw2;
		}}
	delete [] telem;
	#pragma omp critical
	w2= tw2>w2?tw2:w2;
	
	}
	return w2;
}

template<typename dType, typename Type>
void CheckEoms(void (*eom) (Type* elem, Cfunction<dType, Type> *F1, grid<dType>& Rgrid, dType* params, My_Int i, My_Int j), My_Int N_eoms, Cfunction<dType, Type>*Fsol,Cfunction<dType, Type> *Ftest, grid<dType>& Rgrid, grid<dType>& Xgrid, dType*params){

	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();

	#pragma omp parallel
	{Type *telem;
	telem=new Type[N_eoms];
	
	#pragma omp for
	for (My_Int i=1;i<Np-1;i++){
		for (My_Int j=0;j<Nx;j++){
			eom(telem,Fsol,Rgrid,params,i,j);
			for(My_Int l=0; l<N_eoms; l++){
				Ftest[l](i,j)=telem[l];
			}
		}
	}
	
	delete[] telem;}
}

#endif /* defined(__PDEs__NewtonMethod__) */
