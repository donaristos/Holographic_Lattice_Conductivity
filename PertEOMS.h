//
//  TimeEvoEOMS.h
//  TimeDependence
//
//  Created by Aristomenis Donos on 13/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#ifndef __TimeDependence__TimeEvoEOMS__
#define __TimeDependence__TimeEvoEOMS__

#include <iostream>
#include "function.h"
#include "MyTypes.h"

template<typename dType, typename Type>
void equations(Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, dType* params, My_Int i, My_Int j);

template<typename dType, typename Type>
void bconditions(Type* elem, Cfunction<dType, Type> *F1,dType* params, My_Int i, My_Int j);

template<typename dType, typename Type>
void SecOrderEqs(Type* elem, Cfunction<dType, Type> *F1, grid<dType>& grid1, dType* params, My_Int i, My_Int j);

template<>
void equations<double, Complex>(Complex* elem, Cfunction<double, Complex> *F1, grid<double>& grid1, double* params, My_Int i, My_Int j);

template<>
void bconditions<double, Complex>(Complex* elem, Cfunction<double, Complex> *F1,double* params, My_Int i, My_Int j);

template<>
void SecOrderEqs<double, Complex>(Complex* elem, Cfunction<double, Complex> *F1, grid<double>& grid1, double* params, My_Int i, My_Int j);

template<>
void equations<float128, ComplexQ>(ComplexQ* elem, Cfunction<float128, ComplexQ> *F1, grid<float128>& grid1, float128* params, My_Int i, My_Int j);

template<>
void bconditions<float128, ComplexQ>(ComplexQ* elem, Cfunction<float128, ComplexQ> *F1,float128* params, My_Int i, My_Int j);

template<>
void SecOrderEqs<float128, ComplexQ>(ComplexQ* elem, Cfunction<float128, ComplexQ> *F1, grid<float128>& grid1, float128* params, My_Int i, My_Int j);

template<>
void equations<mpreal, ComplexMP>(ComplexMP* elem, Cfunction<mpreal, ComplexMP> *F1, grid<mpreal>& grid1, mpreal* params, My_Int i, My_Int j);

template<>
void bconditions<mpreal, ComplexMP>(ComplexMP* elem, Cfunction<mpreal, ComplexMP> *F1,mpreal* params, My_Int i, My_Int j);

template<>
void SecOrderEqs<mpreal, ComplexMP>(ComplexMP* elem, Cfunction<mpreal, ComplexMP> *F1, grid<mpreal>& grid1, mpreal* params, My_Int i, My_Int j);

template<>
void equations<long double, ComplexLD>(ComplexLD* elem, Cfunction<long double, ComplexLD> *F1, grid<long double>& grid1, long double* params, My_Int i, My_Int j);

template<>
void bconditions<long double, ComplexLD>(ComplexLD* elem, Cfunction<long double, ComplexLD> *F1,long double* params, My_Int i, My_Int j);

template<>
void SecOrderEqs<long double, ComplexLD>(ComplexLD* elem, Cfunction<long double, ComplexLD> *F1, grid<long double>& grid1, long double* params, My_Int i, My_Int j);

#endif /* defined(__TimeDependence__TimeEvoEOMS__) */
