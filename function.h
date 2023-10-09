//
//  function.h
//  PDEs
//
//  Created by Aristomenis Donos on 02/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#ifndef __PDEs__function__
#define __PDEs__function__

#include <iostream>
#include <complex>
#include "math.h"
#include "MyTypes.h"

//Oveload complex

std::complex<double> operator*(int i, std::complex<double> z);

std::complex<double> operator+(int i, std::complex<double> z);

std::complex<double> operator/(int i, std::complex<double> z);

std::complex<double> operator/(std::complex<double> z,int i);

std::complex<double> operator-(int i, std::complex<double> z);

std::complex<float128> operator*(double i, std::complex<float128> z);

std::complex<float128> operator*(std::complex<float128> z,double i);

std::complex<float128> operator+(double i, std::complex<float128> z);

std::complex<float128> operator+(std::complex<float128> z,double i);

std::complex<float128> operator/(double i, std::complex<float128> z);

std::complex<float128> operator/(std::complex<float128> z,double i);

std::complex<float128> operator-(double i, std::complex<float128> z);

std::complex<float128> operator-(std::complex<float128> z,double i);

std::complex<float128> pow(std::complex<float128> z, int i);

std::complex<mpreal> operator*(double i, std::complex<mpreal> z);

std::complex<mpreal> operator*(std::complex<mpreal> z,double i);

std::complex<mpreal> operator+(double i, std::complex<mpreal> z);

std::complex<mpreal> operator+(std::complex<mpreal> z,double i);

std::complex<mpreal> operator/(double i, std::complex<mpreal> z);

std::complex<mpreal> operator/(std::complex<mpreal> z,double i);

std::complex<mpreal> operator-(double i, std::complex<mpreal> z);

std::complex<mpreal> operator-(std::complex<mpreal> z,double i);

std::complex<mpreal> pow(std::complex<mpreal> z, int i);

std::complex<long double> operator*(double i, std::complex<long double> z);

std::complex<long double> operator*(std::complex<long double> z,double i);

std::complex<long double> operator+(double i, std::complex<long double> z);

std::complex<long double> operator+(std::complex<long double> z,double i);

std::complex<long double> operator/(double i, std::complex<long double> z);

std::complex<long double> operator/(std::complex<long double> z,double i);

std::complex<long double> operator-(double i, std::complex<long double> z);

std::complex<long double> operator-(std::complex<long double> z,double i);

std::complex<long double> operator*(int i, std::complex<long double> z);

std::complex<long double> operator+(int i, std::complex<long double> z);

std::complex<long double> operator/(int i, std::complex<long double> z);

std::complex<long double> operator/(std::complex<long double> z,int i);

std::complex<long double> operator-(int i, std::complex<long double> z);

//Overload abs function

template<typename dType, typename Type> dType modulus (Type z);

template<> double modulus<double, double> (double z);

template<> double modulus<double, std::complex<double> > (std::complex<double> z);

template<> float128 modulus<float128, float128> (float128 z);

template<> float128 modulus<float128, std::complex<float128> > (std::complex<float128> z);

template<> mpfr::mpreal modulus<mpfr::mpreal, mpfr::mpreal> (mpfr::mpreal z);

template<> mpfr::mpreal modulus<mpfr::mpreal, std::complex<mpfr::mpreal> > (std::complex<mpfr::mpreal> z);

template<> long double modulus<long double, long double> (long double z);

template<> long double modulus<long double, std::complex<long double> > (std::complex<long double> z);

//Overload long double functions

long double sqrt(long double);

long double exp(long double);

long double fabs(long double);

//long double pow(long double, int);

long double pow(long double, long double);

//long double pow(long double, double);

// Grid class definitions

template<typename Type>
class grid{
private:
	My_Int TotPatches;
	My_Int TotLength;
	My_Int* PatchLeftBoundary;
	My_Int* PatchRightBoundary;
	Type* gridarray;
	My_Int* PatchOrder;
public:
	~grid();
	
	grid(My_Int* Ends, Type *TempGrid, My_Int* Orders, My_Int NPatches=1);
	
	My_Int OrderAtPoint(My_Int i);
	
	My_Int MaxOrder();
	
	My_Int LeftPatchBoundary(My_Int i);
	
	My_Int RightPatchBoundary(My_Int i);
	
	My_Int TotalLength();
	
	bool IsLeftBoundary(My_Int i);
	
	bool IsRightBoundary(My_Int i);
	
	Type operator[](My_Int i);
};

//Derivative class definitions

template<typename Type>
class derivative{
private:
	grid<Type>* Rgrid;
	grid<Type>* Xgrid;
	Type **pr;
	Type **px;
public:
	derivative(Type **Pr=0, Type **Px=0,grid<Type>* TRgrid=0, grid<Type>* TXgrid=0 );
	
	derivative(const derivative& D);
	
	~derivative();
	
	Type operator()(const My_Int i1,const  My_Int i2,const My_Int j1,const  My_Int j2);
	
	derivative& operator=(const derivative& D);
	
	bool getpr(){return (pr==0);};
	bool getpx(){return (px==0);};
	
	My_Int getNp(){return (*Rgrid).TotalLength();}
	My_Int getNx(){return (*Xgrid).TotalLength();}
	
	My_Int Rorder(My_Int i){return (*Rgrid).OrderAtPoint(i);}
	
	My_Int Xorder(My_Int i){return (*Xgrid).OrderAtPoint(i);}
	
	My_Int RLeftBoundary(My_Int i){return (*Rgrid).LeftPatchBoundary(i);}
	
	My_Int XLeftBoundary(My_Int i){return (*Xgrid).LeftPatchBoundary(i);}
	
	My_Int RRightBoundary(My_Int i){return (*Rgrid).RightPatchBoundary(i);}
	
	My_Int XRightBoundary(My_Int i){return (*Xgrid).RightPatchBoundary(i);}
};

//Derivative collection class definitions

template<typename Type>
class derivativeCol{
private:
	grid<Type>* Rgrid;
	grid<Type>* Xgrid;
	My_Int N_ders;
	derivative<Type> *Ders;
public:
	derivativeCol(My_Int n,grid<Type>* TRgrid=0, grid<Type>* TXgrid=0);
	
	~derivativeCol();
	
	derivative<Type>& operator[](const My_Int &i);
	
	My_Int order();
	
	My_Int Rorder(My_Int i){return (*Rgrid).OrderAtPoint(i);}
	
	My_Int Xorder(My_Int i){return (*Xgrid).OrderAtPoint(i);}
	
	My_Int RLeftBoundary(My_Int i){return (*Rgrid).LeftPatchBoundary(i);}
	
	My_Int XLeftBoundary(My_Int i){return (*Xgrid).LeftPatchBoundary(i);}
	
	My_Int RRightBoundary(My_Int i){return (*Rgrid).RightPatchBoundary(i);}
	
	My_Int XRightBoundary(My_Int i){return (*Xgrid).RightPatchBoundary(i);}
	
};

//Function class definitions

template <typename  dType, typename Type>
class function{
private:
	
    My_Int Np;
    My_Int Nx;
    Type **f;
	
public:
	
	function(My_Int rad_points=1, My_Int boundary_points=1);
	
    function(Type **array,My_Int rad_points, My_Int boundary_points);

	function(Type **array,grid<dType> *TRgrid, grid<dType> *TXgrid);
	
	function(const function &f1);
    
    ~function();
    
    Type& operator()(My_Int i, My_Int j);
	
	function& operator=(const function &origf);
	
	function operator+(const function &f1);
	
	function operator-(const function &f1);
	
	template<typename dType1, typename Type1>
	friend function<dType1,Type1> operator*(const Type1& factor, function<dType1,Type1> & f);
	
	template<typename dType1, typename Type1>
	friend std::ostream& operator<<(std::ostream &out, function<dType1,Type1> &f1);
    
    Type diff(derivative<dType> &diff_operator,My_Int i, My_Int j);
	
	function diff(derivative<dType> &diff_operator);
	
};

//Cached function class definitions


template <typename  dType, typename Type>
class Cfunction{
private:
	My_Int nf;
	function<dType, Type> *fs;
public:
	
	Cfunction(My_Int Nf=0);
	
	Cfunction(function<dType,Type> &F,derivativeCol<dType> &D);
	
	~Cfunction();
	
	template<typename dType1, typename Type1>
	friend std::ostream& operator<<(std::ostream &out, Cfunction<dType1,Type1> &f1);
	
	Type& operator()(My_Int i, My_Int j);
	
	Cfunction& operator=(const Cfunction& f1);
	
	Cfunction operator-(const Cfunction& f1);
	
	Cfunction operator+(const Cfunction& f1);
	
	template<typename dType1, typename Type1>
	friend Cfunction<dType1, Type1> operator*(const Type1& factor, Cfunction<dType1, Type1> &f);

	template<typename dType1, typename Type1>
	friend Cfunction<dType1, Type1> operator*(const dType1& factor, Cfunction<dType1, Type1> &f);
	
	Type diff(My_Int di, My_Int i, My_Int j);
	
	function<dType, Type>& diff(My_Int di);
	
	void update(derivativeCol<dType> &D);

	
};

//Grid class implementation

template<typename Type>
grid<Type>::grid(My_Int* Ends, Type *grid, My_Int* Orders, My_Int NPatches){
	TotPatches=NPatches;
	My_Int length(0);
	for (My_Int ii=0; ii<TotPatches; ++ii) {
		length+=Ends[ii];
	}
	TotLength=length;
	PatchLeftBoundary= new My_Int[length];
	PatchRightBoundary= new My_Int[length];
	PatchOrder= new My_Int[length];
	gridarray=new Type[length];
	
	My_Int *tempEnds;
	tempEnds= new My_Int[TotPatches];
	tempEnds[0]=Ends[0];
	for (My_Int ii=1; ii<TotPatches; ++ii) {
		tempEnds[ii]=tempEnds[ii-1]+Ends[ii];
	}
	
	for (My_Int ii=0; ii<Ends[0]; ++ii) {
		PatchLeftBoundary[ii]=0;
		PatchRightBoundary[ii]=tempEnds[0]-1;
		PatchOrder[ii]=Orders[0];
		gridarray[ii]=grid[ii];
	}
	
	for(My_Int jj=1; jj<TotPatches; ++jj){
		for (My_Int ii=tempEnds[jj-1]; ii<tempEnds[jj]; ++ii) {
			gridarray[ii]=grid[ii];
			PatchOrder[ii]=Orders[jj];
			PatchLeftBoundary[ii]=tempEnds[jj-1];
			PatchRightBoundary[ii]=tempEnds[jj]-1;
		}
	}
	
	delete []tempEnds;
}

template<typename Type>
grid<Type>::~grid(){
	delete []PatchLeftBoundary;
	delete []PatchRightBoundary;
	delete []PatchOrder;
	delete []gridarray;
	PatchLeftBoundary=0;
	PatchRightBoundary=0;
	PatchOrder=0;
	gridarray=0;
}

template<typename Type>
My_Int grid<Type>::TotalLength(){
	return TotLength;
}

template<typename Type>
My_Int grid<Type>::LeftPatchBoundary(My_Int i){
	return PatchLeftBoundary[i];
}

template<typename Type>
My_Int grid<Type>::RightPatchBoundary(My_Int i){
	return PatchRightBoundary[i];
}

template<typename Type>
bool grid<Type>::IsLeftBoundary(My_Int i){
	return PatchLeftBoundary[i]==i;
}

template<typename Type>
bool grid<Type>::IsRightBoundary(My_Int i){
	return PatchRightBoundary[i]==i;
}

template<typename Type>
Type grid<Type>::operator[](My_Int i){
	return gridarray[i];
}

template<typename Type>
My_Int grid<Type>::OrderAtPoint(My_Int i){
	return PatchOrder[i];
}

template<typename Type>
My_Int grid<Type>::MaxOrder(){
	My_Int N(0);
	for(My_Int ii=0; ii<TotLength;++ii){
		N=(N>=PatchOrder[ii])?N:PatchOrder[ii];
	}
	return N;
}

// Function class implementation

template<typename dType, typename Type>
function<dType, Type>::function(My_Int rad_points, My_Int boundary_points){
    Np=rad_points;
    Nx=boundary_points;
    f=new Type*[Np];
    for (My_Int i=0; i<Np; i++) {
        f[i]=new Type[Nx];
    }
    for (My_Int i=0; i<Np; i++) {
        for (My_Int j=0; j<Nx; j++) {
            f[i][j]=0;
        }
    }
}

template<typename dType, typename Type>
function<dType, Type>::function(Type **array,My_Int rad_points, My_Int boundary_points){
    Np=rad_points;
    Nx=boundary_points;
    f=new Type*[Np];
    for (My_Int i=0; i<Np; i++) {
        f[i]=new Type[Nx];
    }
    for (My_Int i=0; i<Np; i++) {
        for (My_Int j=0; j<Nx; j++) {
            f[i][j]=array[i][j];
        }
    }
}

template<typename dType, typename Type>
function<dType,Type>::function(Type **array, grid<dType>* TRgrid, grid<dType>* TXgrid){
    Np=(*TRgrid).TotalLength();
    Nx=(*TXgrid).TotalLength();
    f=new Type*[Np];
    for (My_Int i=0; i<Np; i++) {
        f[i]=new Type[Nx];
    }
    for (My_Int i=0; i<Np; i++) {
        for (My_Int j=0; j<Nx; j++) {
            f[i][j]=array[i][j];
        }
    }
}



template<typename dType, typename Type>
function<dType,Type>::function(const function<dType,Type> &f1){
	Np=f1.Np;
	Nx= f1.Nx;
	f = new Type*[Np];
	for (My_Int i=0; i<Np; i++) {
		f[i]=new Type[Nx];
	}
	for (My_Int i=0; i<f1.Np; i++) {
		for (My_Int j=0; j<Nx; j++){
			f[i][j]=f1.f[i][j];
		}}
}

template<typename dType, typename Type>
function<dType,Type>::~function(){
    for (My_Int i=0; i<Np; i++) {
        delete []f[i];
    }
    delete []f;
    f=0;
}

template<typename dType, typename Type>
Type& function<dType,Type>::operator()(My_Int i, My_Int j){return f[i][j];};


template<typename dType, typename Type>
function<dType, Type> function<dType,Type>::operator+(const function<dType,Type> &f1){
	Type **temp;
    temp=new Type*[Np];
    for (My_Int i=0; i<Np; i++) {
        temp[i]=new Type[Nx];
    }
	
    for (My_Int i=0; i<Np; i++) {
        for (My_Int j=0; j<Nx; j++) {
            temp[i][j]=f[i][j]+f1.f[i][j];
        }
	}
	function<dType, Type> tempf(temp,Np,Nx);
	
	for (My_Int i=0; i<Np; i++) {
		delete[] temp[i];
	}
	
	delete [] temp;
	
	return tempf;
};

template<typename dType, typename Type>
function<dType,Type> function<dType,Type>::operator-(const function<dType,Type> &f1){
	Type **temp;
    temp=new Type*[Np];
    for (My_Int i=0; i<Np; i++) {
        temp[i]=new Type[Nx];
    }
	
    for (My_Int i=0; i<Np; i++) {
        for (My_Int j=0; j<Nx; j++) {
            temp[i][j]=f[i][j]-f1.f[i][j];
        }
	}
	function<dType,Type> tempf(temp,Np,Nx);
	
	for (My_Int i=0; i<Np; i++) {
		delete[] temp[i];
	}
	
	delete [] temp;
	
	return tempf;
};

template<typename dType1, typename Type1>
function<dType1, Type1> operator*(const Type1& x,function<dType1,Type1>& f){
	function<dType1, Type1> tempf(f.Np,f.Nx);
	for (My_Int i=0; i<tempf.Np; ++i) {
		for (My_Int j=0; j<tempf.Nx; ++j) {
			tempf(i,j)=x*f(i,j);
		}
	}
	return tempf;
}

template<typename dType, typename Type>
function<dType, Type>& function<dType, Type>::operator=(const function<dType, Type> &origf){
	if (this==&origf)
		return *this;
	
	for (My_Int i=0; i<Np; i++) {
        delete []f[i];
    }
    delete []f;
    f=0;
	
	Np=origf.Np;
	Nx=origf.Nx;
	f=new Type*[Np];
    for (My_Int i=0; i<Np; i++) {
        f[i]=new Type[Nx];
    }
	
	for (My_Int i=0; i<Np; i++){
		for (My_Int j=0; j<Nx; j++) {
			f[i][j]=(origf.f)[i][j];
		}
	}
	return *this;
}

template<typename dType, typename Type>
std::ostream& operator<<(std::ostream &out,function<dType, Type> &f1){
	out << "{ ";
	for (My_Int i=0; i<f1.Np; i++) {
		out<< "{ ";
		for (My_Int j=0; j<f1.Nx; j++) {
			out<< f1(i,j);
			if (j<f1.Nx-1) {out << ", ";}
		}
		if(i<f1.Np-1){out<< "}, ";}
		out << std::endl;}
	out<< "}}" << std::endl;
	return out;
}

template<typename dType, typename Type>
Type function<dType, Type>::diff(derivative<dType> &diff_operator, My_Int i, My_Int j){
	
	Type Sum=Type(0);

	if(diff_operator.getpr()){
		for (My_Int j1=0; j1<Nx; ++j1) {
			Sum=Sum+diff_operator(i,j,i,j1)*f[i][j1];
		}
	}
	else if(diff_operator.getpx()){
		My_Int lower=((i-diff_operator.Rorder(i))>=diff_operator.RLeftBoundary(i))?(i-diff_operator.Rorder(i)):diff_operator.RLeftBoundary(i);
		My_Int upper=((i+diff_operator.Rorder(i))<=diff_operator.RRightBoundary(i))?(i+diff_operator.Rorder(i)):diff_operator.RRightBoundary(i);
		for (My_Int i1=lower ; i1<=upper; i1++) {
			Sum=Sum+diff_operator(i,j,i1,j)*f[i1][j];
		}
		
	}
	else{
		Type* sum;
		
		sum=new Type[Np];
		
		My_Int lower=(i-diff_operator.Rorder(i))>=diff_operator.RLeftBoundary(i)?(i-diff_operator.Rorder(i)):diff_operator.RLeftBoundary(i);
		My_Int upper=(i+diff_operator.Rorder(i))<=diff_operator.RRightBoundary(i)?(i+diff_operator.Rorder(i)):diff_operator.RRightBoundary(i);
		
		for(My_Int i1=lower; i1<=upper;i1++){sum[i1]=0;}
		
		for(My_Int i1=0;i1<Np;i1++){
			for(My_Int j1=0; j1<Nx;++j1){
				sum[i1]+=diff_operator(i,j,i1,j1)*f[i1][j1];
			}
		}
		
		
		
		
		for(My_Int i1=lower;i1<upper;i1++){Sum+=sum[i1];}
		
		delete[] sum;}
	
	return Sum;
};


template<typename dType, typename Type>
function<dType, Type> function<dType, Type>::diff(derivative<dType> &diff_operator){
	Type **temp;
    temp=new Type*[Np];
    for (My_Int i=0; i<Np; i++) {
        temp[i]=new Type[Nx];
    }
	if(diff_operator.getpr()||diff_operator.getpx()){
		for (My_Int i=0; i<Np; i++) {
			for (My_Int j=0; j<Nx; j++) {
				temp[i][j]=(*this).diff(diff_operator, i, j);
			}
		}}
	else{
		
		
		Type ***temp2;
		temp2= new Type**[Np];
		for (My_Int i=0; i<Np; i++) {
			temp2[i]=new Type*[Nx];
			for (My_Int j=0; j<Nx; ++j) {
				temp2[i][j]=new Type[Np];
			}
		}

		for (My_Int i1=0; i1<Np; ++i1) {
			for (My_Int i2=0; i2<Nx; ++i2) {
				for (My_Int j1=0; j1<Np; ++j1) {
					temp2[i1][i2][j1]=0;
					if( __builtin_abs(i1-j1) <= diff_operator.Rorder(i1) && diff_operator.RLeftBoundary(i1)==diff_operator.RLeftBoundary(j1)){
					for(My_Int j2=0;j2<Nx;++j2){
						temp2[i1][i2][j1]+= diff_operator(i1,i2,j1,j2)*f[j1][j2];
					}}
				}
			}
		}
		for (My_Int i1=0; i1<Np; ++i1) {
			for (My_Int i2=0; i2<Nx; ++i2) {
				temp[i1][i2]=0;
				for (My_Int j1=((i1-diff_operator.Rorder(i1))>=diff_operator.RLeftBoundary(i1)?(i1-diff_operator.Rorder(i1)):diff_operator.RLeftBoundary(i1)); j1<=((i1+diff_operator.Rorder(i1))<=diff_operator.RRightBoundary(i1)?(i1+diff_operator.Rorder(i1)):diff_operator.RRightBoundary(i1)); ++j1) {
					temp[i1][i2]+=temp2[i1][i2][j1];
				}}}
		for (My_Int i1=0; i1<Np; ++i1) {
			for(My_Int i2=0; i2<Nx; ++i2){
				delete[] temp2[i1][i2];
			}
		}
		for (My_Int i1=0; i1<Np; ++i1) {
			delete[] temp2[i1];
		}
		delete[] temp2;
		
	}
	function tempf(temp,Np,Nx);
	
	for (My_Int i=0; i<Np; i++) {
		delete[] temp[i];
	}
	
	delete [] temp;
	
	return tempf;
};





//Cached function class implementation

template<typename dType, typename Type>
Cfunction<dType, Type>::Cfunction(My_Int Nf){
	nf=Nf;
	if(nf==0){fs=0;}
	else{
		fs=new function<dType, Type>[Nf];
	}
}

template<typename dType, typename Type>
Cfunction<dType, Type>::Cfunction(function<dType, Type> &F,derivativeCol<dType> &D){
	nf=1+D.order();
	fs=new function<dType, Type>[nf];
	fs[0]=F;
	#pragma omp parallel for
	for (My_Int i=0; i<nf-1; ++i) {
		fs[i+1]=F.diff(D[i]);
	}
}

template<typename dType, typename Type>
Cfunction<dType, Type>::~Cfunction(){
	if(fs!=0){delete [] fs;}
	fs=0;
}

template<typename dType, typename Type>
Type& Cfunction<dType, Type>::operator()(My_Int i, My_Int j){return fs[0](i,j);}

template<typename dType, typename Type>
Type Cfunction<dType, Type>::diff(My_Int di, My_Int i, My_Int j){return fs[di](i,j);}

template<typename dType, typename Type>
function<dType, Type>& Cfunction<dType, Type>::diff(My_Int di){
	return fs[di];
}

template<typename dType, typename Type>
void Cfunction<dType, Type>::update(derivativeCol<dType> &D){
	#pragma omp parallel for
	for (My_Int i=0; i<nf-1; ++i) {
		fs[i+1]=fs[0].diff(D[i]);
	}
}

template<typename dType, typename Type>
Cfunction<dType, Type>& Cfunction<dType, Type>::operator=(const Cfunction<dType,Type> &f1){
	if (fs!=0){delete [] fs;};
	nf=f1.nf;
	fs=new function<dType, Type>[nf];
	for (My_Int i=0; i<nf; ++i) {
		fs[i]=(f1.fs)[i];
	}
	return *this;
}

template<typename dType, typename Type>
Cfunction<dType, Type> Cfunction<dType, Type>::operator-(const Cfunction<dType, Type>& f1){
	Cfunction<dType, Type> tempCf(f1.nf);
	for (My_Int i=0; i<tempCf.nf; ++i) {
		tempCf.fs[i]=fs[i]-f1.fs[i];
	}
	tempCf.MaxOrder= f1.MaxOrder ;
	return tempCf;
}

template<typename dType1, typename Type1>
Cfunction<dType1, Type1> operator*(const Type1& x, Cfunction<dType1, Type1>& f1){
	Cfunction<dType1, Type1> tempCf(f1.nf);
	for (My_Int i=0; i<tempCf.nf; ++i) {
		tempCf.fs[i]=x*f1.fs[i];
	}
	
	return tempCf;
}

template<typename dType1, typename Type1>
Cfunction<dType1, Type1> operator*(const dType1& x, Cfunction<dType1, Type1>& f1){
	Cfunction<dType1, Type1> tempCf(f1.nf);
	for (My_Int i=0; i<tempCf.nf; ++i) {
		tempCf.fs[i]=Type1(x)*f1.fs[i];
	}
	
	return tempCf;
}

template<typename dType, typename Type>
std::ostream& operator<<(std::ostream &out, Cfunction<dType, Type> &f1){
	out << (f1.fs)[0];
	return out;
}

//Derivative class implementation

template<typename Type>
derivative<Type>::derivative(Type **Pr, Type **Px,grid<Type>* TRgrid, grid<Type>* TXgrid){
	Rgrid=TRgrid;
	Xgrid=TXgrid;
	My_Int Np=(Rgrid==0)?0:(*Rgrid).TotalLength();
	My_Int Nx=(Xgrid==0)?0:(*Xgrid).TotalLength();
	
	pr = new Type*[Np];
	px = new Type*[Nx];
	
	if (Pr==0&&Px==0){
		delete [] pr;
		pr=0;
		delete [] px;
		px=0;
	}
	else if (Pr==0) {
        delete [] pr;
		pr=0;
		for(My_Int i=0; i<Nx; i++) {
			px[i]=new Type[Nx];
		}
		for (My_Int i=0; i<Nx ; i++) {
			for (My_Int j=0; j<Nx; j++) {
				px[i][j]=Px[i][j];
			}}
	}
	else if (Px==0)
	{   delete [] px;
		px=0;
		for(My_Int i=0; i<Np; i++) {
			pr[i]=new Type[Np];
		}
		for (My_Int i=0; i<Np ; i++) {
			for (My_Int j=0; j<Np; j++) {
				pr[i][j]=Pr[i][j];
			}}
	}
	else{
		for(My_Int i=0; i<Nx; i++) {
			px[i]=new Type[Nx];
		}
		for (My_Int i=0; i<Nx ; i++) {
			for (My_Int j=0; j<Nx; j++) {
				px[i][j]=Px[i][j];
			}}
		for(My_Int i=0; i<Np; i++) {
			pr[i]=new Type[Np];
		}
		for (My_Int i=0; i<Np ; i++) {
			for (My_Int j=0; j<Np; j++) {
				pr[i][j]=Pr[i][j];
			}}
		
	}
}

template<typename Type>
derivative<Type>::derivative(const derivative<Type>& D){
	Rgrid=D.Rgrid;
	Xgrid=D.Xgrid;
	My_Int Np=(*Rgrid).TotalLength();
	My_Int Nx=(*Xgrid).TotalLength();
	pr = new Type*[Np];
	px = new Type*[Nx];
	
	if (D.pr==0) {
        delete [] pr;
		pr=0;
		for(My_Int i=0; i<Nx; i++) {
			px[i]=new Type[Nx];
		}
		for (My_Int i=0; i<Nx ; i++) {
			for (My_Int j=0; j<Nx; j++) {
				px[i][j]=D.px[i][j];
			}}
	}
	else if (D.px==0)
	{   delete [] px;
		px=0;
		for(My_Int i=0; i<Np; i++) {
			pr[i]=new Type[Np];
		}
		for (My_Int i=0; i<Np ; i++) {
			for (My_Int j=0; j<Np; j++) {
				pr[i][j]=D.pr[i][j];
			}}
	}
	else{
		for(My_Int i=0; i<Nx; i++) {
			px[i]=new Type[Nx];
		}
		for (My_Int i=0; i<Nx ; i++) {
			for (My_Int j=0; j<Nx; j++) {
				px[i][j]=D.px[i][j];
			}}
		for(My_Int i=0; i<Np; i++) {
			pr[i]=new Type[Np];
		}
		for (My_Int i=0; i<Np ; i++) {
			for (My_Int j=0; j<Np; j++) {
				pr[i][j]=D.pr[i][j];
			}}
		
	}
}

template<typename Type>
derivative<Type>::~derivative(){
	My_Int Np=(Rgrid==0)?0:(*Rgrid).TotalLength();
	My_Int Nx=(Xgrid==0)?0:(*Xgrid).TotalLength();
	if (pr==0&&px==0){}
	else if (pr==0) {
		for (My_Int i=0; i<Nx; i++) {
			delete [] px[i];
		}
		delete[] px;
	}
	else if (px==0){
		for (My_Int i=0; i<Np; i++) {
			delete [] pr[i];
		}
		delete[] pr;
	}
	else {
		for (My_Int i=0; i<Nx; i++) {
			delete [] px[i];
		}
		for (My_Int i=0; i<Np; i++) {
			delete [] pr[i];
		}
		delete[] px;
		delete[] pr;
	}
	pr=0;
	px=0;
	Rgrid=0;
	Xgrid=0;
}


template<typename Type>
Type derivative<Type>::operator()(const My_Int i1,const My_Int i2,const My_Int j1,const My_Int j2){
	if(pr==0&&px==0)
    {return (i1==j1?1:0)* (i2==j2?1:0);}
	else if(pr==0)
    {return (i1==j1?1:0)* px[i2][j2];}
	else if(px==0)
    {return (i2==j2?1:0)*(pr[i1][j1]);}
	else if(pr!=0&&px!=0)
    {return pr[i1][j1]*px[i2][j2];}
    else
    {std::cout << "Weird!" << std::endl;}
}

template<typename Type>
derivative<Type>& derivative<Type>::operator=(const derivative<Type>& D){
	if (this==&D)
		return *this;
	
	
	if (pr==0&&px==0){}
	else if (pr==0) {
		for (My_Int i=0; i<(*Xgrid).TotalLength(); i++) {
			delete [] px[i];
		}
		delete[] px;
	}
	else if (px==0){
		for (My_Int i=0; i<(*Rgrid).TotalLength(); i++) {
			delete [] pr[i];
		}
		delete[] pr;
	}
	else {
		for (My_Int i=0; i<(*Xgrid).TotalLength(); i++) {
			delete [] px[i];
		}
		for (My_Int i=0; i<(*Rgrid).TotalLength(); i++) {
			delete [] pr[i];
		}
		delete[] px;
		delete[] pr;
	}
	pr=0;
	px=0;
	
	Rgrid=D.Rgrid;
	Xgrid=D.Xgrid;
	
	My_Int Np=(*Rgrid).TotalLength();
	My_Int Nx=(*Xgrid).TotalLength();
	pr = new Type*[Np];
	px = new Type*[Nx];
	
	if (D.pr==0) {
        delete [] pr;
		pr=0;
		for(My_Int i=0; i<Nx; i++) {
			px[i]=new Type[Nx];
		}
		for (My_Int i=0; i<Nx ; i++) {
			for (My_Int j=0; j<Nx; j++) {
				px[i][j]=D.px[i][j];
			}}
	}
	else if (D.px==0)
	{   delete [] px;
		px=0;
		for(My_Int i=0; i<Np; i++) {
			pr[i]=new Type[Np];
		}
		for (My_Int i=0; i<Np ; i++) {
			for (My_Int j=0; j<Np; j++) {
				pr[i][j]=D.pr[i][j];
			}}
	}
	else{
		for(My_Int i=0; i<Nx; i++) {
			px[i]=new Type[Nx];
		}
		for (My_Int i=0; i<Nx ; i++) {
			for (My_Int j=0; j<Nx; j++) {
				px[i][j]=D.px[i][j];
			}}
		for(My_Int i=0; i<Np; i++) {
			pr[i]=new Type[Np];
		}
		for (My_Int i=0; i<Np ; i++) {
			for (My_Int j=0; j<Np; j++) {
				pr[i][j]=D.pr[i][j];
			}}
		
	}
	return *this;
}


//Derivative collection class implementation

template<typename Type>
derivativeCol<Type>::derivativeCol(My_Int n,grid<Type>* TRgrid, grid<Type>* TXgrid ){
	Rgrid=TRgrid;
	Xgrid=TXgrid;
	N_ders=n;
	Ders= new derivative<Type>[N_ders];
}

template<typename Type>
derivativeCol<Type>::~derivativeCol(){
	delete[] Ders;
	Ders=0;
	Rgrid=0;
	Xgrid=0;
}

template<typename Type>
derivative<Type>& derivativeCol<Type>::operator[](const My_Int &i){return Ders[i];}

template<typename Type>
My_Int derivativeCol<Type>::order(){
	return N_ders;
}


#endif /* defined(__PDEs__function__) */