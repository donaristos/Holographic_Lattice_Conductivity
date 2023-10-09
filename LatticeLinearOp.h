//
//  InitLinearOp.h
//  TimeDependence
//
//  Created by Aristomenis Donos on 17/06/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#ifndef __TimeDependence__InitLinearOp__
#define __TimeDependence__InitLinearOp__

#include <iostream>
#include "function.h"
#include "MyTypes.h"


template<typename dType, typename Type>
void BulkLinearOp(Type **Elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, grid<dType>& grid1,dType* params,My_Int i1, My_Int j1, My_Int i, My_Int j);

template<typename dType, typename Type>
void BoundaryLinearOp(Type **Elem, derivativeCol<dType>& Der,Cfunction<dType, Type> *F1, Type* params,  My_Int i1, My_Int j1, My_Int i, My_Int j);

template<>
void BulkLinearOp<double, double>(double **elem, derivativeCol<double>& D,Cfunction<double, double> *F1, grid<double>& grid1,  double* params, My_Int i, My_Int j, My_Int i1, My_Int j1);

template<>
void BoundaryLinearOp<double, double>(double **elem, derivativeCol<double>& D,Cfunction<double, double> *F1, double* params,  My_Int i, My_Int j, My_Int i1, My_Int j1);

template<>
void BulkLinearOp<float128, float128>(float128 **elem, derivativeCol<float128>& D,Cfunction<float128, float128> *F1, grid<float128>& grid1,  float128* params, My_Int i, My_Int j, My_Int i1, My_Int j1);

template<>
void BoundaryLinearOp<float128, float128>(float128 **elem, derivativeCol<float128>& D,Cfunction<float128, float128> *F1, float128* params,  My_Int i, My_Int j, My_Int i1, My_Int j1);

template<>
void BulkLinearOp<mpreal, mpreal>(mpreal **elem, derivativeCol<mpreal>& D,Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1,  mpreal* params, My_Int i, My_Int j, My_Int i1, My_Int j1);

template<>
void BoundaryLinearOp<mpreal, mpreal>(mpreal **elem, derivativeCol<mpreal>& D,Cfunction<mpreal, mpreal> *F1, mpreal* params,  My_Int i, My_Int j, My_Int i1, My_Int j1);

template<>
void BulkLinearOp<long double, long double>(long double **elem, derivativeCol<long double>& D,Cfunction<long double, long double> *F1, grid<long double>& grid1,  long double* params, My_Int i, My_Int j, My_Int i1, My_Int j1);

template<>
void BoundaryLinearOp<long double, long double>(long double **elem, derivativeCol<long double>& D,Cfunction<long double, long double> *F1, long double* params,  My_Int i, My_Int j, My_Int i1, My_Int j1);

#endif /* defined(__TimeDependence__InitLinearOp__) */
