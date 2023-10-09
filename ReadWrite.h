//
//  ReadWrite.h
//  TimeDependence
//
//  Created by Aristomenis Donos on 13/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#ifndef __TimeDependence__ReadWrite__
#define __TimeDependence__ReadWrite__

#include <iostream>
#include <fstream>
#include <vector>
#include "function.h"
#include "MyTypes.h"


void ReadIntParameters(My_Int *parameters, My_Int length,const std::string& file_path);

//General templates

template<typename Type>
void ReadStamp(Type &t, My_Int i,const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	if(is.is_open()){
		is.seekg(i*std::streamsize(sizeof(Type)),is.beg);
		is.read(reinterpret_cast<char *>(&t), std::streamsize(sizeof(Type)));}
}


template<typename Type>
My_Int GetLength(My_Int CellLength,const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	long begin(0),end(0);
	if (is.is_open()) {
		is.seekg(0,is.beg);
		begin = is.tellg();
		is.seekg(0,is.end);
		end= is.tellg();
	}
	return (end-begin)/(CellLength*(std::streamsize(sizeof(Type))));
}

template<typename Type>
void ReadDoubleParameters(Type *parameters, My_Int length, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	if(is.is_open()){
//		std::cout << "It is open now!"<< std::endl;
		is.seekg(0,is.beg);
		is.read(reinterpret_cast<char *>(parameters),std::streamsize(length*sizeof(Type)));
	}
}


template<typename Type>
void ReadDiffFromBin(Type** diff_m,My_Int dim, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	if(is.is_open()){
//		std::cout << "Diff file is open" << std::endl;
		is.seekg(0,is.beg);
		for (My_Int ii=0; ii<dim; ++ii) {
			is.read(reinterpret_cast<char *>(diff_m[ii]), std::streamsize(dim*sizeof(Type)));
		}
	}
}


template <typename Type>
void WriteStamp(Type &t,const std::string& file_path, My_Int prec=80){
	std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out|std::ios::binary);
	os.write(reinterpret_cast<char *>(&t), std::streamsize(sizeof(Type)));
	os.close();
}


template<typename dType, typename Type>
void WriteToBin(Cfunction<dType, Type> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec=30){
	Type ***tarray;
	
	tarray=new Type**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new Type*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new Type[Nx];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			for(My_Int k=0; k<Nx;++k){
				tarray[i][j][k]=F[ind[i]](j,k);
			}
		}
	}
	
    std::ofstream os(file_path.c_str(),std::ios::out|std::ios::binary);
	
	for (My_Int i=0; i<Nf; ++i) {
		for(My_Int j=0; j<Np; ++j){
            os.write(reinterpret_cast<char *>(tarray[i][j]), std::streamsize(Nx*sizeof(Type)));
		}
	}
    
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			delete [] tarray[i][j];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
        delete [] tarray[i];
    }
    
    delete [] tarray;
    
    os.close();
}


template<typename dType, typename Type>
void AppendToBin(Cfunction<dType, Type> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec=30){
	Type ***tarray;
	
	tarray=new Type**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new Type*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new Type[Nx];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			for(My_Int k=0; k<Nx;++k){
				tarray[i][j][k]=F[ind[i]](j,k);
			}
		}
	}
	
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out|std::ios::binary);
	
	for (My_Int i=0; i<Nf; ++i) {
		for(My_Int j=0; j<Np; ++j){
            os.write(reinterpret_cast<char *>(tarray[i][j]), std::streamsize(Nx*sizeof(Type)));
		}
	}
    
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			delete [] tarray[i][j];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
        delete [] tarray[i];
    }
    
    delete [] tarray;
    
    os.close();
}

template<typename dType, typename Type>
void ReadFromBin(Cfunction<dType, Type> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int pos=-1){
	dType ***tarray;
	
	tarray=new dType**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new dType*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new dType[Nx];
		}
	}
	
    std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	
	if (is.is_open()) {
		if (pos<0) {
			is.seekg(pos*std::streamsize(Nf*Np*Nx*sizeof(Type)),is.end);
		}
		else
		{
			is.seekg(pos*std::streamsize(Nf*Np*Nx*sizeof(Type)),is.beg);
		}
		
		for (My_Int i=0; i<Nf; ++i) {
			for(My_Int j=0; j<Np; ++j){
				is.read(reinterpret_cast<char *>(tarray[i][j]), std::streamsize(Nx*sizeof(dType)));
			}
		}
		
		for (My_Int i=0; i<Nf; ++i){
			for (My_Int j=0; j<Np; ++j) {
				for(My_Int k=0; k<Nx;++k){
					F[ind[i]](j,k)=tarray[i][j][k];
				}
			}
		}}
    
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			delete [] tarray[i][j];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
        delete [] tarray[i];
    }
    
    delete [] tarray;
    
    is.close();
}

//Specialization to mpreal type

template<>
void ReadStamp<mpreal>(mpreal &t, My_Int i,const std::string& file_path);


template<>
My_Int GetLength<mpreal>(My_Int CellLength,const std::string& file_path);

template<>
void ReadDoubleParameters<mpreal>(mpreal *parameters, My_Int length, const std::string& file_path);


template<>
void ReadDiffFromBin<mpreal>(mpreal** diff_m,My_Int dim, const std::string& file_path);


template<>
void WriteStamp<mpreal>(mpreal &t,const std::string& file_path, My_Int prec);

template<>
void WriteToBin<mpreal, mpreal>(Cfunction<mpreal, mpreal> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec);

template<>
void AppendToBin<mpreal, mpreal>(Cfunction<mpreal, mpreal> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec);

template<>
void ReadFromBin<mpreal, mpreal>(Cfunction<mpreal, mpreal> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int pos);

template<>
void WriteToBin<mpreal, ComplexMP>(Cfunction<mpreal, ComplexMP> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec);

template<>
void AppendToBin<mpreal, ComplexMP>(Cfunction<mpreal, ComplexMP> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec);

template<>
void ReadFromBin<mpreal, ComplexMP>(Cfunction<mpreal, ComplexMP> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int pos);

//Specialization to long double type

template<>
void ReadStamp<long double>(long double &t, My_Int i,const std::string& file_path);


template<>
My_Int GetLength<long double>(My_Int CellLength,const std::string& file_path);

template<>
void ReadDoubleParameters<long double>(long double *parameters, My_Int length, const std::string& file_path);


template<>
void ReadDiffFromBin<long double>(long double** diff_m,My_Int dim, const std::string& file_path);


template<>
void WriteStamp<long double>(long double &t,const std::string& file_path, My_Int prec);

template<>
void WriteToBin<long double, long double>(Cfunction<long double, long double> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec);

template<>
void AppendToBin<long double, long double>(Cfunction<long double, long double> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec);

template<>
void ReadFromBin<long double, long double>(Cfunction<long double, long double> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int pos);

template<>
void WriteToBin<long double, ComplexLD>(Cfunction<long double, ComplexLD> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec);

template<>
void AppendToBin<long double, ComplexLD>(Cfunction<long double, ComplexLD> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec);

template<>
void ReadFromBin<long double, ComplexLD>(Cfunction<long double, ComplexLD> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int pos);

#endif /* defined(__TimeDependence__ReadWrite__) */
