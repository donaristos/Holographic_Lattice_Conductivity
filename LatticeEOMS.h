//
//  EOMS.h
//  PDEs
//
//  Created by Aristomenis Donos on 05/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#ifndef __PDEs__EOMS__
#define __PDEs__EOMS__

#include <iostream>
#include "function.h"
#include "MyTypes.h"

template<typename dType, typename Type>
void equations(Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, dType* params, My_Int i, My_Int j);

template<typename dType, typename Type>
void bconditions(Type* elem, Cfunction<dType, Type> *F1,dType* params, My_Int i, My_Int j);

template<typename dType, typename Type>
void ksi2(Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, dType* params, My_Int i, My_Int j);

template<>
void equations<double, double>(double* elem, Cfunction<double, double> *F1, grid<double>& grid1,  double* params, My_Int i, My_Int j);

template<>
void bconditions<double, double>(double* elem, Cfunction<double, double> *F1,double* params,  My_Int i, My_Int j);

template<>
void ksi2<double, double>(double* elem, Cfunction<double, double> *F1, grid<double>& grid1, double* params, My_Int i, My_Int j);

template<>
void equations<float128, float128>(float128* elem, Cfunction<float128, float128> *F1, grid<float128>& grid1,  float128* params, My_Int i, My_Int j);

template<>
void bconditions<float128, float128>(float128* elem, Cfunction<float128, float128> *F1,float128* params,  My_Int i, My_Int j);

template<>
void ksi2<float128, float128>(float128* elem, Cfunction<float128, float128> *F1, grid<float128>& grid1, float128* params, My_Int i, My_Int j);

template<>
void equations<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1,  mpreal* params, My_Int i, My_Int j);

template<>
void bconditions<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1,mpreal* params,  My_Int i, My_Int j);

template<>
void ksi2<mpreal, mpreal>(mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid1, mpreal* params, My_Int i, My_Int j);

template<>
void equations<long double, long double>(long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1,  long double* params, My_Int i, My_Int j);

template<>
void bconditions<long double, long double>(long double* elem, Cfunction<long double, long double> *F1,long double* params,  My_Int i, My_Int j);

template<>
void ksi2<long double, long double>(long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid1, long double* params, My_Int i, My_Int j);

#endif /* defined(__PDEs__EOMS__) */
