//
//  TimeEvoLinearOp.h
//  TimeDependence
//
//  Created by Aristomenis Donos on 17/06/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#ifndef __TimeDependence__TimeEvoLinearOp__
#define __TimeDependence__TimeEvoLinearOp__

//#include <iostream>
//#include "MyTypes.h"
#include "function.h"

template<typename dType, typename Type>
void BulkLinearOp(Type **Elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, grid<dType>& grid1, dType* params, My_Int i1, My_Int j1, My_Int i, My_Int j);

template<typename dType, typename Type>
void BoundaryLinearOp(Type **Elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, dType* params, My_Int i1, My_Int j1, My_Int i, My_Int j);

template<>
void BulkLinearOp<double, std::complex<double> >(std::complex<double> **Elem, derivativeCol<double>& Der,Cfunction<double, std::complex<double> > *F1, grid<double>& grid1, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j);

template<>
void BoundaryLinearOp<double, Complex>(Complex **Elem, derivativeCol<double>& Der,Cfunction<double, Complex> *F1, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j);

template<>
void BulkLinearOp<float128, ComplexQ>(ComplexQ **Elem, derivativeCol<float128>& Der,Cfunction<float128, ComplexQ> *F1, grid<float128>& grid1, float128* params, My_Int i1, My_Int j1, My_Int i, My_Int j);

template<>
void BoundaryLinearOp<float128, ComplexQ>(ComplexQ **Elem, derivativeCol<float128>& Der,Cfunction<float128, ComplexQ> *F1, float128* params, My_Int i1, My_Int j1, My_Int i, My_Int j);

template<>
void BulkLinearOp<mpreal, ComplexMP>(ComplexMP **Elem, derivativeCol<mpreal>& Der,Cfunction<mpreal, ComplexMP> *F1, grid<mpreal>& grid1, mpreal* params, My_Int i1, My_Int j1, My_Int i, My_Int j);

template<>
void BoundaryLinearOp<mpreal, ComplexMP>(ComplexMP **Elem, derivativeCol<mpreal>& Der,Cfunction<mpreal, ComplexMP> *F1, mpreal* params, My_Int i1, My_Int j1, My_Int i, My_Int j);

template<>
void BulkLinearOp<long double, ComplexLD>(ComplexLD **Elem, derivativeCol<long double>& Der,Cfunction<long double, ComplexLD> *F1, grid<long double>& grid1, long double* params, My_Int i1, My_Int j1, My_Int i, My_Int j);

template<>
void BoundaryLinearOp<long double, ComplexLD>(ComplexLD **Elem, derivativeCol<long double>& Der,Cfunction<long double, ComplexLD> *F1, long double* params, My_Int i1, My_Int j1, My_Int i, My_Int j);

#endif /* defined(__TimeDependence__TimeEvoLinearOp__) */
