//
//  ReadWrite.cpp
//  TimeDependence
//
//  Created by Aristomenis Donos on 13/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "ReadWrite.h"

void ReadIntParameters(My_Int *parameters, My_Int length, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	if(is.is_open()){
//		std::cout << "It is open now!"<< std::endl;
		is.seekg(0,is.beg);
		is.read(reinterpret_cast<char *>(parameters),std::streamsize(length*sizeof(My_Int)));
	}
}

//Implementation of full mpreal specialization

template<>
void ReadStamp<mpreal>(mpreal &t, My_Int i,const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	std::vector<mpreal> tempstamp;
	if(is.is_open()){
		mpreal temp;
		while(is >> temp){
			tempstamp.push_back(temp);
		}
		t=tempstamp[i];
	}
}


template<>
My_Int GetLength<mpreal>(My_Int CellLength,const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	My_Int i(0);
	if(is.is_open()){
		mpreal temp;
		while(is >> temp){
			++i;}
	}
	return i;
}

template<>
void ReadDoubleParameters<mpreal>(mpreal *parameters, My_Int length, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in);
	if(is.is_open()){
//		std::cout << "It is open now!"<< std::endl;
		is.seekg(0,is.beg);
		for (My_Int ii=0; ii<length; ++ii) {
			is >> parameters[ii];
		}
	}
}


template<>
void ReadDiffFromBin<mpreal>(mpreal** diff_m,My_Int dim, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in);
	if(is.is_open()){
//		std::cout << "Diff file is open" << std::endl;
		is.seekg(0,is.beg);
		for (My_Int ii=0; ii<dim; ++ii) {
			for (My_Int jj=0; jj<dim; ++jj){
				is >> diff_m[ii][jj];
			}}
	}
}


template<>
void WriteStamp<mpreal>(mpreal &t,const std::string& file_path, My_Int prec){
	std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out);
	os << std::fixed;
	os.precision(prec);
	os <<  t << "\t" ;
	os.close();
}

template<>
void WriteToBin<mpreal, mpreal>(Cfunction<mpreal, mpreal> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec){
	mpreal ***tarray;
	
	tarray=new mpreal**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new mpreal*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new mpreal[Nx];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			for(My_Int k=0; k<Nx;++k){
				tarray[i][j][k]=F[ind[i]](j,k);
			}
		}
	}
	
    std::ofstream os(file_path.c_str(),std::ios::out);
	
	os << std::fixed;
	os.precision(prec);
	
	for (My_Int i=0; i<Nf; ++i) {
		for(My_Int j=0; j<Np; ++j){
			for(My_Int k=0; k<Nx; ++k){
				os << tarray[i][j][k] << "\t";
			}
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

template<>
void AppendToBin<mpreal, mpreal>(Cfunction<mpreal, mpreal> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec){
	mpreal ***tarray;
	
	tarray=new mpreal**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new mpreal*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new mpreal[Nx];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			for(My_Int k=0; k<Nx;++k){
				tarray[i][j][k]=F[ind[i]](j,k);
			}
		}
	}
	
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out);
	
	os << std::fixed;
	os.precision(prec);
	
	
	for (My_Int i=0; i<Nf; ++i) {
		for(My_Int j=0; j<Np; ++j){
			for(My_Int k=0; k<Nx; ++k){
				os << tarray[i][j][k] << "\t";
			}
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

template<>
void ReadFromBin<mpreal, mpreal>(Cfunction<mpreal, mpreal> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int pos){
	mpreal ***tarray;
	
	tarray=new mpreal**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new mpreal*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new mpreal[Nx];
		}
	}
	
    std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	
	if(is.is_open()){
		
		is.seekg(0, is.beg);
		
		mpreal tempvar;
		std::vector<mpreal> tempvec;
		My_Int length;
		My_Int position;
		
		while(is >> tempvar){
			tempvec.push_back(tempvar);
		}
		
		length=tempvec.size();
		
		if (pos<0) {
			position=length+ pos*Nx*Np*Nf;
		}
		else
		{
			position=pos*Nx*Np*Nf;
		}
		
		for (My_Int i=0; i<Nf; ++i) {
			for(My_Int j=0; j<Np; ++j){
				for(My_Int k=0; k<Nx; ++k){
					tarray[i][j][k]=tempvec[position+i*Np*Nx+j*Nx+k];
				}
			}
		}
		
		for (My_Int i=0; i<Nf; ++i){
			for (My_Int j=0; j<Np; ++j) {
				for(My_Int k=0; k<Nx;++k){
					F[ind[i]](j,k)=tarray[i][j][k];
				}
			}
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
    
    is.close();
}

template<>
void WriteToBin<mpreal, ComplexMP>(Cfunction<mpreal, ComplexMP> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec){
	ComplexMP ***tarray;
	
	tarray=new ComplexMP**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new ComplexMP*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new ComplexMP[Nx];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			for(My_Int k=0; k<Nx;++k){
				tarray[i][j][k]=F[ind[i]](j,k);
			}
		}
	}
	
    std::ofstream os(file_path.c_str(),std::ios::out);
	
	os << std::fixed;
	os.precision(prec);
	
	for (My_Int i=0; i<Nf; ++i) {
		for(My_Int j=0; j<Np; ++j){
			for(My_Int k=0; k<Nx; ++k){
				os << tarray[i][j][k] << "\t";
			}
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

template<>
void AppendToBin<mpreal, ComplexMP>(Cfunction<mpreal, ComplexMP> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec){
	ComplexMP ***tarray;
	
	tarray=new ComplexMP**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new ComplexMP*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new ComplexMP[Nx];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			for(My_Int k=0; k<Nx;++k){
				tarray[i][j][k]=F[ind[i]](j,k);
			}
		}
	}
	
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out);
	
	os << std::fixed;
	os.precision(prec);
	
	for (My_Int i=0; i<Nf; ++i) {
		for(My_Int j=0; j<Np; ++j){
			for(My_Int k=0; k<Nx; ++k){
				os << tarray[i][j][k] << "\t";
			}
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

template<>
void ReadFromBin<mpreal, ComplexMP>(Cfunction<mpreal, ComplexMP> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int pos){
	ComplexMP ***tarray;
	
	tarray=new ComplexMP**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new ComplexMP*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new ComplexMP[Nx];
		}
	}
	
    std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	
	if(is.is_open()){
		
		is.seekg(0, is.beg);
		
		ComplexMP tempvar;
		std::vector<ComplexMP> tempvec;
		My_Int length;
		My_Int position;
		
		while(is >> tempvar){
			tempvec.push_back(tempvar);
		}
		
		length=tempvec.size();
		
		if (pos<0) {
			position=length+ pos*Nx*Np*Nf;
		}
		else
		{
			position=pos*Nx*Np*Nf;
		}
		
		for (My_Int i=0; i<Nf; ++i) {
			for(My_Int j=0; j<Np; ++j){
				for(My_Int k=0; k<Nx; ++k){
					tarray[i][j][k]=tempvec[position+i*Np*Nx+j*Nx+k];
				}
			}
		}
		
		for (My_Int i=0; i<Nf; ++i){
			for (My_Int j=0; j<Np; ++j) {
				for(My_Int k=0; k<Nx;++k){
					F[ind[i]](j,k)=tarray[i][j][k];
				}
			}
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
    
    is.close();
}

//Implementation of full long double specialization

template<>
void ReadStamp<long double>(long double &t, My_Int i,const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	std::vector<long double> tempstamp;
	if(is.is_open()){
		long double temp;
		while(is >> temp){
			tempstamp.push_back(temp);
		}
		t=tempstamp[i];
	}
}


template<>
My_Int GetLength<long double>(My_Int CellLength,const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	My_Int i(0);
	if(is.is_open()){
		long double temp;
		while(is >> temp){
			++i;}
	}
	return i;
}

template<>
void ReadDoubleParameters<long double>(long double *parameters, My_Int length, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in);
	if(is.is_open()){
//		std::cout << "It is open now!"<< std::endl;
		is.seekg(0,is.beg);
		for (My_Int ii=0; ii<length; ++ii) {
			is >> parameters[ii];
		}
	}
}


template<>
void ReadDiffFromBin<long double>(long double** diff_m,My_Int dim, const std::string& file_path){
	std::ifstream is(file_path.c_str(),std::ios::in);
	if(is.is_open()){
//		std::cout << "Diff file is open" << std::endl;
		is.seekg(0,is.beg);
		for (My_Int ii=0; ii<dim; ++ii) {
			for (My_Int jj=0; jj<dim; ++jj){
				is >> diff_m[ii][jj];
			}}
	}
}


template<>
void WriteStamp<long double>(long double &t,const std::string& file_path, My_Int prec){
	std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out);
	os << std::fixed;
	os.precision(prec);
	os <<  t << "\t" ;
	os.close();
}

template<>
void WriteToBin<long double, long double>(Cfunction<long double, long double> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec){
	long double ***tarray;
	
	tarray=new long double**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new long double*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new long double[Nx];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			for(My_Int k=0; k<Nx;++k){
				tarray[i][j][k]=F[ind[i]](j,k);
			}
		}
	}
	
    std::ofstream os(file_path.c_str(),std::ios::out);
	
	os << std::fixed;
	os.precision(prec);
	
	for (My_Int i=0; i<Nf; ++i) {
		for(My_Int j=0; j<Np; ++j){
			for(My_Int k=0; k<Nx; ++k){
				os << tarray[i][j][k] << "\t";
			}
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

template<>
void AppendToBin<long double, long double>(Cfunction<long double, long double> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec){
	long double ***tarray;
	
	tarray=new long double**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new long double*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new long double[Nx];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			for(My_Int k=0; k<Nx;++k){
				tarray[i][j][k]=F[ind[i]](j,k);
			}
		}
	}
	
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out);
	
	os << std::fixed;
	os.precision(prec);
	
	
	for (My_Int i=0; i<Nf; ++i) {
		for(My_Int j=0; j<Np; ++j){
			for(My_Int k=0; k<Nx; ++k){
				os << tarray[i][j][k] << "\t";
			}
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

template<>
void ReadFromBin<long double, long double>(Cfunction<long double, long double> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int pos){
	long double ***tarray;
	
	tarray=new long double**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new long double*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new long double[Nx];
		}
	}
	
    std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	
	if(is.is_open()){
		
		is.seekg(0, is.beg);
		
		long double tempvar;
		std::vector<long double> tempvec;
		My_Int length;
		My_Int position;
		
		while(is >> tempvar){
			tempvec.push_back(tempvar);
		}
		
		length=tempvec.size();
		
		if (pos<0) {
			position=length+ pos*Nx*Np*Nf;
		}
		else
		{
			position=pos*Nx*Np*Nf;
		}
		
		for (My_Int i=0; i<Nf; ++i) {
			for(My_Int j=0; j<Np; ++j){
				for(My_Int k=0; k<Nx; ++k){
					tarray[i][j][k]=tempvec[position+i*Np*Nx+j*Nx+k];
				}
			}
		}
		
		for (My_Int i=0; i<Nf; ++i){
			for (My_Int j=0; j<Np; ++j) {
				for(My_Int k=0; k<Nx;++k){
					F[ind[i]](j,k)=tarray[i][j][k];
				}
			}
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
    
    is.close();
}

template<>
void WriteToBin<long double, ComplexLD>(Cfunction<long double, ComplexLD> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec){
	ComplexLD ***tarray;
	
	tarray=new ComplexLD**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new ComplexLD*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new ComplexLD[Nx];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			for(My_Int k=0; k<Nx;++k){
				tarray[i][j][k]=F[ind[i]](j,k);
			}
		}
	}
	
    std::ofstream os(file_path.c_str(),std::ios::out);
	
	os << std::fixed;
	os.precision(prec);
	
	for (My_Int i=0; i<Nf; ++i) {
		for(My_Int j=0; j<Np; ++j){
			for(My_Int k=0; k<Nx; ++k){
				os << tarray[i][j][k] << "\t";
			}
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

template<>
void AppendToBin<long double, ComplexLD>(Cfunction<long double, ComplexLD> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int prec){
	ComplexLD ***tarray;
	
	tarray=new ComplexLD**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new ComplexLD*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new ComplexLD[Nx];
		}
	}
	
	for (My_Int i=0; i<Nf; ++i){
		for (My_Int j=0; j<Np; ++j) {
			for(My_Int k=0; k<Nx;++k){
				tarray[i][j][k]=F[ind[i]](j,k);
			}
		}
	}
	
    std::ofstream os(file_path.c_str(),std::ios::app|std::ios::out);
	
	os << std::fixed;
	os.precision(prec);
	
	for (My_Int i=0; i<Nf; ++i) {
		for(My_Int j=0; j<Np; ++j){
			for(My_Int k=0; k<Nx; ++k){
				os << tarray[i][j][k] << "\t";
			}
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

template<>
void ReadFromBin<long double, ComplexLD>(Cfunction<long double, ComplexLD> * F,My_Int * ind, My_Int Nf, My_Int Np, My_Int Nx,const std::string& file_path, My_Int pos){
	ComplexLD ***tarray;
	
	tarray=new ComplexLD**[Nf];
	
	for (My_Int i=0; i<Nf; ++i) {
		tarray[i]=new ComplexLD*[Np];
		for (My_Int j=0; j<Np; ++j) {
			tarray[i][j]=new ComplexLD[Nx];
		}
	}
	
    std::ifstream is(file_path.c_str(),std::ios::in|std::ios::binary);
	
	if(is.is_open()){
		
		is.seekg(0, is.beg);
		
		ComplexLD tempvar;
		std::vector<ComplexLD> tempvec;
		My_Int length;
		My_Int position;
		
		while(is >> tempvar){
			tempvec.push_back(tempvar);
		}
		
		length=tempvec.size();
		
		if (pos<0) {
			position=length+ pos*Nx*Np*Nf;
		}
		else
		{
			position=pos*Nx*Np*Nf;
		}
		
		for (My_Int i=0; i<Nf; ++i) {
			for(My_Int j=0; j<Np; ++j){
				for(My_Int k=0; k<Nx; ++k){
					tarray[i][j][k]=tempvec[position+i*Np*Nx+j*Nx+k];
				}
			}
		}
		
		for (My_Int i=0; i<Nf; ++i){
			for (My_Int j=0; j<Np; ++j) {
				for(My_Int k=0; k<Nx;++k){
					F[ind[i]](j,k)=tarray[i][j][k];
				}
			}
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
    
    is.close();
}