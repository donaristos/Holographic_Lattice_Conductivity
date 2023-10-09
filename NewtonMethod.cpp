//
//  NewtonMethod.cpp
//  PDEs
//
//  Created by Aristomenis Donos on 04/05/2013.
//  Copyright (c) 2013 Aristomenis Donos. All rights reserved.
//

#include "NewtonMethod.h"

//using namespace mpfr;

//template<>
//void UpdateFunctions<double, double>( void (*eom) (double* elem, Cfunction<double, double> *F1, grid<double>& grid, double* params, My_Int i, My_Int j)
//									 ,void (*b_cond) (double* elem, Cfunction<double, double> *F1,double* params, My_Int i, My_Int j)
//									 ,void (*l_eom)(double **elem, derivativeCol<double>& Der,Cfunction<double, double> *F1, grid<double>& grid, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
//									 ,void (*l_b_cond)(double **elem,  derivativeCol<double>& Der,Cfunction<double, double> *F1, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<double, double>* f, derivativeCol<double>& D, grid<double>& Rgrid,grid<double>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds, double* params){
//	
//	
//	My_Int Np=Rgrid.TotalLength();
//	My_Int Nx=Xgrid.TotalLength();
//	My_Int nbdof = (Np-1)*Nx*N_eoms+(N_bconds-N_eoms)*Nx; // number of degrees of freedom.
//	My_Int nvars = (Rgrid.MaxOrder()*N_eoms+(N_bconds-N_eoms))*Xgrid.MaxOrder();
//	
//	
//	std::vector<Eigen::Triplet<double, My_Int> > TriplVec;
//	Eigen::Matrix<double, Eigen::Dynamic, 1> B(nbdof);
//	
//	TriplVec.reserve(nvars*nbdof);
//	
//	ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid, N_functions,N_eoms,N_bconds,params, 1.E-13);
//	
//	std::cout << "Reserved: " << nvars*nbdof << ", Needed: " << TriplVec.size() << std::endl;
//	
//	long nz=TriplVec.size();
//	
//	void *Symbolic, *Numeric;
//	double *Info, *Control;
//	
//	Info=new double[UMFPACK_INFO];
//	Control=new double[UMFPACK_CONTROL];
//	
//	/* get the default control parameters */
//    umfpack_dl_defaults (Control) ;
//	
//    /* change the default print level for this demo */
//    /* (otherwise, nothing will print) */
//    Control [UMFPACK_PRL] = 6 ;
//	
//    /* print the license agreement */
//	//    umfpack_zi_report_status (Control, UMFPACK_OK) ;
//    Control [UMFPACK_PRL] = 5 ;
//	
//	double *x;
//	x= new double[nbdof];
//	
//	long *Ap;
//	Ap= new long[nbdof+1];
//	
//	long *Ai;
//	Ai= new long[nz];
//	
//	double *Ax;
//	Ax= new double[nz];
//	
//	long *Ti;
//	Ti= new long[nz];
//	
//	long *Tj;
//	Tj= new long[nz];
//	
//	double *Tx;
//	Tx= new double[nz];
//	
//	double *b;
//	b= new double[nbdof];
//	
//	long status;
//	
//	for(My_Int ii=0;ii<nbdof;++ii){
//		b[ii]=B[ii];
//		x[ii]=0;
//	}
//	
//	
//	for(My_Int ii=0;ii<nz;++ii){
//		Ti[ii]=TriplVec[ii].row();
//		Tj[ii]=TriplVec[ii].col();
//		Tx[ii]=TriplVec[ii].value();
//	}
//	
//	TriplVec.erase(TriplVec.begin(), TriplVec.end());
//	
//	std::cout << "Starting solver with dimension " << nbdof << " x " << nbdof << " matrices and with " << nz << " non-zero elements"<< std::endl;
//	
//	status=umfpack_dl_triplet_to_col(nbdof, nbdof, nz, Ti, Tj, Tx,Ap, Ai, Ax, NULL);
//	
//	if (status!=0) {
//		std::cout<< status <<std::endl;
//	}
//	
//	umfpack_dl_symbolic(nbdof, nbdof, Ap, Ai, Ax, &Symbolic, Control, Info);
//	
//	status=umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
//	
//	if (status!=0) {
//		std::cout<< status <<std::endl;
//	}
//	
//	umfpack_dl_free_symbolic(&Symbolic);
//	
//	umfpack_dl_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL);
//	
//	umfpack_dl_free_numeric(&Numeric);
//	
//	for (My_Int l3=0; l3<N_eoms; l3++) {
//		for (My_Int i1=1; i1 < Np; i1++) {
//			for (My_Int i2=0; i2<Nx; i2++) {
//				
//				f[l3](i1,i2)=f[l3](i1,i2)-x[i2+(i1-1)*Nx+l3*(Np-1)*Nx];
//			}}}
//    
//	for (My_Int l3=0; l3<N_bconds-N_eoms;++l3){
//		for(My_Int i2=0; i2<Nx;++i2){
//			
//			f[l3](0,i2)=f[l3](0,i2)-x[N_eoms*(Np-1)*Nx+l3*Nx+i2];
//		}
//	}
//	
//	
//	for (My_Int l3=0; l3<N_functions; ++l3) {
//		f[l3].update(D);
//	}
//	
//	std::cout << "Called the <double, double> function!" << std::endl;
//	std::cout << "Matrix Condition Number=" << Info[UMFPACK_RCOND] << std::endl;
//	
//	for(My_Int ii=0; ii<nz;++ii){
//		b[Ti[ii]]-=Tx[ii]*x[Tj[ii]];
//	}
//	
//	double res(0);
//	
//	for(My_Int ii=0;ii<nbdof;++ii){
//		res=(fabs(b[ii])>res)?fabs(b[ii]):res;
//	}
//	
//	std::cout << "Residue estimate: " << res << std::endl;
//	
//	delete []x;
//	delete []Ap;
//	delete []Ai;
//	delete []Ax;
//	delete []Ti;
//	delete []Tj;
//	delete []Tx;
//	delete []b;
//	delete []Info;
//	delete []Control;
//
//}
//
//template<>
//void UpdateFunctions<double, Complex>( void (*eom) (Complex* elem, Cfunction<double, Complex> *F1 , grid<double>& grid, double* params, My_Int i, My_Int j)
//									  ,void (*b_cond) (Complex* elem, Cfunction<double, Complex> *F1,double* params, My_Int i, My_Int j)
//									  ,void (*l_eom)(Complex **elem,  derivativeCol<double>& Der,Cfunction<double, Complex> *F1, grid<double>& grid, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
//									  ,void (*l_b_cond)(Complex **elem,  derivativeCol<double>& Der,Cfunction<double, Complex> *F1, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<double, Complex>* f, derivativeCol<double>& D, grid<double>& Rgrid, grid<double>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds,double* params){
//	
//	My_Int Np=Rgrid.TotalLength();
//	My_Int Nx=Xgrid.TotalLength();
//	My_Int nbdof = (Np-1)*Nx*N_eoms+(N_bconds-N_eoms)*Nx; // number of degrees of freedom.
//	My_Int nvars = (Rgrid.MaxOrder()*N_eoms+(N_bconds-N_eoms))*Xgrid.MaxOrder();
//	
//	
//	std::vector<Eigen::Triplet<Complex, My_Int> > TriplVec;
//	Eigen::Matrix<Complex, Eigen::Dynamic, 1> B(nbdof);
//	
//	TriplVec.reserve(nvars*nbdof);
//	
//	ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid,N_functions,N_eoms,N_bconds,params, 1.E-13);
//	
//	std::cout << "Reserved: " << nvars*nbdof << ", Needed: " << TriplVec.size() << std::endl;
//	
//	My_Int nz=TriplVec.size();
//	
//	void *Symbolic, *Numeric;
//	long status;
//	double *Info, *Control;
//	
//	Info=new double[UMFPACK_INFO];
//	Control=new double[UMFPACK_CONTROL];
//	
//	/* get the default control parameters */
//    umfpack_zi_defaults (Control) ;
//	
//    /* change the default print level for this demo */
//    /* (otherwise, nothing will print) */
//    Control [UMFPACK_PRL] = 6 ;
//	
//    /* print the license agreement */
//	//    umfpack_zi_report_status (Control, UMFPACK_OK) ;
//    Control [UMFPACK_PRL] = 5 ;
//	
//    /* print the control parameters */
//	//    umfpack_di_report_control (Control) ;
//	
//	//	Control[UMFPACK_PIVOT_TOLERANCE]=0.8;
//	//	Control[UMFPACK_SYM_PIVOT_TOLERANCE]=0.8;
//	//	Control[UMFPACK_ORDERING]=UMFPACK_ORDERING_BEST;
//	//	Control[UMFPACK_IRSTEP]=20;
//	
//	double *xr;
//	xr= new double[nbdof];
//	
//	double *xi;
//	xi= new double[nbdof];
//	
//	long *Ap;
//	Ap= new long[nbdof+1];
//	
//	long *Ai;
//	Ai= new long[nz];
//	
//	double *Ax;
//	Ax= new double[nz];
//	
//	double *Az;
//	Az= new double[nz];
//	
//	long *Ti;
//	Ti= new long[nz];
//	
//	long *Tj;
//	Tj= new long[nz];
//	
//	double *Tx;
//	Tx= new double[nz];
//	
//	double *Tz;
//	Tz= new double[nz];
//	
//	double *br;
//	br= new double[nbdof];
//	
//	double *bi;
//	bi= new double[nbdof];
//	
//	
//	for(My_Int ii=0;ii<nbdof;++ii){
//		br[ii]=B[ii].real();
//		bi[ii]=B[ii].imag();
//		xr[ii]=0;
//	}
//	
//	
//	for(My_Int ii=0;ii<nz;++ii){
//		Ti[ii]=TriplVec[ii].row();
//		Tj[ii]=TriplVec[ii].col();
//		Tx[ii]=TriplVec[ii].value().real();
//		Tz[ii]=TriplVec[ii].value().imag();
//	}
//	
//	TriplVec.erase(TriplVec.begin(), TriplVec.end());
//	
//	std::cout << "Starting solver with dimension "<< nbdof << " x " << nbdof << " matrices and with " <<nz << " non-zero elements"<< std::endl;
//	
//	status=umfpack_zl_triplet_to_col(nbdof, nbdof, nz, Ti, Tj, Tx, Tz, Ap, Ai, Ax, Az, NULL);
//	
//	if (status!=0) {
//		std::cout<< status <<std::endl;
//	}
//	
//	status=umfpack_zl_symbolic(nbdof, nbdof, Ap, Ai, Ax,Az, &Symbolic, Control, Info);
//	
//	if (status!=0) {
//		std::cout<< status <<std::endl;
//	}
//	
//	status=umfpack_zl_numeric(Ap, Ai, Ax,Az, Symbolic, &Numeric, Control, Info);
//	
//	
//	if (status!=0) {
//		std::cout << status << std::endl;
//		umfpack_zl_report_info(Control, Info) ;
//		umfpack_zl_report_status(Control, status) ;}
//	
//	umfpack_zl_free_symbolic(&Symbolic);
//	
//	umfpack_zl_solve(UMFPACK_A, Ap, Ai, Ax, Az, xr, xi, br, bi, Numeric, Control, Info);
//	
//	umfpack_zl_free_numeric(&Numeric);
//	
//	
//	for (My_Int l3=0; l3<N_eoms; l3++) {
//		for (My_Int i1=1; i1 < Np; i1++) {
//			for (My_Int i2=0; i2<Nx; i2++) {
//				
//				f[l3](i1,i2)=f[l3](i1,i2)-Complex(xr[i2+(i1-1)*Nx+l3*(Np-1)*Nx],xi[i2+(i1-1)*Nx+l3*(Np-1)*Nx]);
//			}}}
//    
//	for (My_Int l3=0; l3<N_bconds-N_eoms;++l3){
//		for(My_Int i2=0; i2<Nx;++i2){
//			
//			f[l3](0,i2)=f[l3](0,i2)-Complex(xr[N_eoms*(Np-1)*Nx+l3*Nx+i2],xi[N_eoms*(Np-1)*Nx+l3*Nx+i2]);
//		}
//	}
//	
//	
//	for (My_Int l3=0; l3<N_functions; ++l3) {
//		f[l3].update(D);
//	}
//	
//	std::cout << "Called the <double, complex> function!" << std::endl;
//	std::cout << "Matrix Condition Number=" << Info[UMFPACK_RCOND] << std::endl;
//	
//	delete []xr;
//	delete []xi;
//	delete []Ap;
//	delete []Ai;
//	delete []Ax;
//	delete []Az;
//	delete []Ti;
//	delete []Tj;
//	delete []Tx;
//	delete []Tz;
//	delete []br;
//	delete []bi;
//	delete []Info;
//	delete []Control;
//}


template<>
void UpdateFunctions<double, double>( void (*eom) (double* elem, Cfunction<double, double> *F1, grid<double>& grid, double* params, My_Int i, My_Int j)
									 ,void (*b_cond) (double* elem, Cfunction<double, double> *F1,double* params, My_Int i, My_Int j)
									 ,void (*l_eom)(double **elem, derivativeCol<double>& Der,Cfunction<double, double> *F1, grid<double>& grid, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
									 ,void (*l_b_cond)(double **elem,  derivativeCol<double>& Der,Cfunction<double, double> *F1, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<double, double>* f, derivativeCol<double>& D, grid<double>& Rgrid,grid<double>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds, double* params){
	
	
	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
    My_Int nbdof = (Np-1)*Nx*N_eoms+(N_bconds-N_eoms)*Nx; // number of degrees of freedom.
	My_Int nvars = (Rgrid.MaxOrder()*N_eoms+(N_bconds-N_eoms))*Xgrid.MaxOrder();
	
	
	std::vector<Eigen::Triplet<double, My_Int> > TriplVec;
	Eigen::SparseMatrix<double, Eigen::ColMajor, MKL_INT> MSp(nbdof,nbdof);
	Eigen::Matrix<double, Eigen::Dynamic, 1> B(nbdof);
	Eigen::Matrix<double, Eigen::Dynamic, 1> X(nbdof);
	Eigen::Matrix<double, Eigen::Dynamic, 1> Test(nbdof);
	
	TriplVec.reserve(nvars*nbdof);
	
	ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid, N_functions,N_eoms,N_bconds,params, 1.E-13);
	
    
	std::cout << "Reserved: " << nvars*nbdof << ", Needed: " << TriplVec.size() << std::endl;
	
	std::cout << "Starting solver with dimension "<< nbdof << " x " << nbdof << " matrices and with " <<TriplVec.size() << " non-zero elements"<< std::endl;
    
	MSp.setFromTriplets(TriplVec.begin(),TriplVec.end());
    
	TriplVec.erase(TriplVec.begin(), TriplVec.end());
    
    MKL_INT mtype = 11;       /* Real unsymmetric matrix */
    
    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: My_Int pt[64]; 64-bit: long My_Int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    
    /* PARDISO control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    double ddum;          /* Double dummy */
    MKL_INT idum(0);         /* Integer dummy. */
    double res, res0;
    MKL_INT nn=nbdof;
    
    for (My_Int i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 3;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 2;        /* Conjugate transposed/transpose solve */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = -1;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = 0;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = 0;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    
    iparm[20] = 0;  // 1x1 pivoting
    iparm[26] = 0;  // No matrix checker
    //    iparm[27] = (sizeof(RealScalar) == 4) ? 1 : 0;
    iparm[34] = 1;         // C indexing
    iparm[59] = 0;          // InCore Calculation
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;         /* Which factorization to use. */
    msglvl = 0;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */
    
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for (My_Int i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
    
    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    phase = 12;
    

    
    PARDISO_64 (pt, &maxfct, &mnum, &mtype, &phase,
             &nn, MSp.valuePtr() , MSp.outerIndexPtr() , MSp.innerIndexPtr() , &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    

    
    if ( error != 0 )
    {
        std::cout << "Error number" << error << " during analysis and numerical factorization" << std::endl;
    }
    
    std::cout << "Peak memory during symbolic factorization: "<< iparm[14] << " KBytes" <<std::endl;
    
    phase = 33;
    
    
    PARDISO_64 (pt, &maxfct, &mnum, &mtype, &phase,
             &nn, MSp.valuePtr() , MSp.outerIndexPtr()  , MSp.innerIndexPtr(), &idum, &nrhs, iparm, &msglvl, B.data(), X.data() ,&error);
    
    if ( error != 0 )
    {
        std::cout << "Error number" << error << " during solving and refinement" << std::endl;
    }
    
    
    phase = -1;
    
    PARDISO_64 (pt, &maxfct, &mnum, &mtype, &phase,
             &nn, MSp.valuePtr()  , MSp.outerIndexPtr() , MSp.innerIndexPtr(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    
	for (My_Int l3=0; l3<N_eoms; l3++) {
		for (My_Int i1=1; i1 < Np; i1++) {
			for (My_Int i2=0; i2<Nx; i2++) {
				
				f[l3](i1,i2)=f[l3](i1,i2)-X[i2+(i1-1)*Nx+l3*(Np-1)*Nx];
			}}}
    
	for (My_Int l3=0; l3<N_bconds-N_eoms;++l3){
		for(My_Int i2=0; i2<Nx;++i2){
			
			f[l3](0,i2)=f[l3](0,i2)-X[N_eoms*(Np-1)*Nx+l3*Nx+i2];
		}
	}
	
	
	for (My_Int l3=0; l3<N_functions; ++l3) {
		f[l3].update(D);
	}
	
	std::cout << "Called the <double, double> function!" << std::endl;
	
	Test=B-MSp*X;
	
    res=0;
	
	for(My_Int ii=0;ii<nbdof;++ii){
		res=(fabs(Test[ii])>res)?fabs(Test[ii]):res;
	}
	
	std::cout << "Residue estimate: " << res << std::endl;
    
}

template<>
void UpdateFunctions<double, Complex>( void (*eom) (Complex* elem, Cfunction<double, Complex> *F1 , grid<double>& grid, double* params, My_Int i, My_Int j)
									  ,void (*b_cond) (Complex* elem, Cfunction<double, Complex> *F1,double* params, My_Int i, My_Int j)
									  ,void (*l_eom)(Complex **elem,  derivativeCol<double>& Der,Cfunction<double, Complex> *F1, grid<double>& grid, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
									  ,void (*l_b_cond)(Complex **elem,  derivativeCol<double>& Der,Cfunction<double, Complex> *F1, double* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<double, Complex>* f, derivativeCol<double>& D, grid<double>& Rgrid, grid<double>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds,double* params){
	
	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
	My_Int nbdof = (Np-1)*Nx*N_eoms+(N_bconds-N_eoms)*Nx; // number of degrees of freedom.
	My_Int nvars = (Rgrid.MaxOrder()*N_eoms+(N_bconds-N_eoms))*Xgrid.MaxOrder();
	
	
	std::vector<Eigen::Triplet< Complex , My_Int> > TriplVec;
	Eigen::SparseMatrix<Complex , Eigen::ColMajor, MKL_INT> MSp(nbdof,nbdof);
	Eigen::Matrix<Complex, Eigen::Dynamic, 1> B(nbdof);
	Eigen::Matrix<Complex, Eigen::Dynamic, 1> Test(nbdof);
    Eigen::Matrix<Complex, Eigen::Dynamic, 1> X(nbdof);
    
	TriplVec.reserve(nvars*nbdof);
	
	ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid,N_functions,N_eoms,N_bconds,params, 1.E-13);

	
    MSp.setFromTriplets(TriplVec.begin(), TriplVec.end());
    
	std::cout << "Reserved: " << nvars*nbdof << ", Needed: " << TriplVec.size() << std::endl;
	
	std::cout << "Starting solver with dimension "<< nbdof << " x " << nbdof << " matrices and with " << TriplVec.size() << " non-zero elements"<< std::endl;
    
    TriplVec.erase(TriplVec.begin(), TriplVec.end());
    
    MSp.makeCompressed();
    
    MKL_INT mtype = 13;       /* Complex unsymmetric matrix */
    
    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: My_Int pt[64]; 64-bit: long My_Int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    
    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    Complex ddum;          /* Double dummy */
    MKL_INT idum(0);         /* Integer dummy. */
    double res, res0;
    MKL_INT nn=nbdof;
    
    for (My_Int i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 0;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 20;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 2;        /* Conjugate transposed/transpose solve */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = -1;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = 0;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = 0;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    
    iparm[20] = 1;  // 1x1 pivoting
    iparm[23] = 1;  // Parallel factorization control
    iparm[26] = 0;  // No matrix checker
    iparm[34] = 1;         // C indexing
    iparm[59] = 0;          // InCore Calculation
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;         /* Which factorization to use. */
    msglvl = 0;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */
    
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for (My_Int i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
    
    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    phase = 12;
    
    PARDISO_64 (pt, &maxfct, &mnum, &mtype, &phase,
             &nn, MSp.valuePtr() , MSp.outerIndexPtr() , MSp.innerIndexPtr() , &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    
    if ( error != 0 )
    {
        std::cout << "Error number" << error << " during analysis and numerical factorization" << std::endl;
    }
    
    std::cout << "Peak memory during symbolic factorization: "<< iparm[14] << " KBytes" <<std::endl;

    
    phase = 33;
    
    
    PARDISO_64 (pt, &maxfct, &mnum, &mtype, &phase,
             &nn, MSp.valuePtr() , MSp.outerIndexPtr()  , MSp.innerIndexPtr(), &idum, &nrhs, iparm, &msglvl, B.data(), X.data() ,&error);
    
    if ( error != 0 )
    {
        std::cout << "Error number" << error << " during solving and refinement" << std::endl;
    }
    
    std::cout << "Iterative steps: "<< iparm[6] <<std::endl;
    
    phase = -1;
    
    PARDISO_64 (pt, &maxfct, &mnum, &mtype, &phase,
             &nn, MSp.valuePtr()  , MSp.outerIndexPtr() , MSp.innerIndexPtr(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

    
	for (My_Int l3=0; l3<N_eoms; l3++) {
		for (My_Int i1=1; i1 < Np; i1++) {
			for (My_Int i2=0; i2<Nx; i2++) {
				
				f[l3](i1,i2)=f[l3](i1,i2)-X[i2+(i1-1)*Nx+l3*(Np-1)*Nx];
			}}}
    
	for (My_Int l3=0; l3<N_bconds-N_eoms;++l3){
		for(My_Int i2=0; i2<Nx;++i2){
			
			f[l3](0,i2)=f[l3](0,i2)- X[N_eoms*(Np-1)*Nx+l3*Nx+i2];
		}
	}
	
	
	for (My_Int l3=0; l3<N_functions; ++l3) {
		f[l3].update(D);
	}
	
	std::cout << "Called the <double, double> function!" << std::endl;
	
	Test=B-MSp*X;
	
    res=0;
	
	for(My_Int ii=0;ii<nbdof;++ii){
		res=(abs(Test[ii])>res)?abs(Test[ii]):res;
	}
    
    std::cout << "Residue estimate: " << res << std::endl;
}


template<>
void UpdateFunctions<float128,float128>( void (*eom) (float128* elem, Cfunction<float128, float128> *F1, grid<float128>& grid, float128* params, My_Int i, My_Int j)
										,void (*b_cond) (float128* elem, Cfunction<float128, float128> *F1,float128* params, My_Int i, My_Int j)
										,void (*l_eom)(float128 **elem, derivativeCol<float128>& Der,Cfunction<float128, float128> *F1, grid<float128>& grid, float128* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
										,void (*l_b_cond)(float128 **elem,  derivativeCol<float128>& Der,Cfunction<float128, float128> *F1, float128* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<float128, float128>* f, derivativeCol<float128>& D, grid<float128>& Rgrid,grid<float128>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds, float128* params){
	
	
	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
	My_Int nbdof = (Np-1)*Nx*N_eoms+(N_bconds-N_eoms)*Nx; // number of degrees of freedom.
	My_Int nvars = (Rgrid.MaxOrder()*N_eoms+(N_bconds-N_eoms))*Xgrid.MaxOrder();
	
	std::vector<Eigen::Triplet<float128, My_Int> > TriplVec;
	Eigen::SparseMatrix<float128, Eigen::ColMajor, My_Int> MSp(nbdof,nbdof);
	Eigen::Matrix<float128, Eigen::Dynamic, 1> B(nbdof);
	Eigen::Matrix<float128, Eigen::Dynamic, 1> X(nbdof);
	Eigen::Matrix<float128, Eigen::Dynamic, 1> Test(nbdof);
	
	TriplVec.reserve(nvars*nbdof);
    
    float128 tol=5.E-30;
	
	ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid, N_functions,N_eoms,N_bconds,params, tol );
	
	std::cout << "Reserved: " << nvars*nbdof << ", Needed: " << TriplVec.size() << std::endl;
	
	std::cout << "Starting solver with dimension "<< nbdof << " x " << nbdof << " matrices and with " <<TriplVec.size() << " non-zero elements"<< std::endl;
	
	MSp.setFromTriplets(TriplVec.begin(),TriplVec.end());
	
	TriplVec.erase(TriplVec.begin(), TriplVec.end());
	
	MSp.makeCompressed();
	
	//	Eigen::UmfPackLU<Eigen::SparseMatrix<float128> > solver;
	//	Eigen::BiCGSTAB<Eigen::SparseMatrix<float128>,Eigen::IncompleteLUT<float128>> solver(MSp);
	Eigen::SparseLU<Eigen::SparseMatrix<float128, Eigen::ColMajor, My_Int>,Eigen::COLAMDOrdering<My_Int> > solver(MSp);
	//	solver.compute(MSp);
	std::cout<< solver.info() << std::endl;
	X=solver.solve(B);
	std::cout<< solver.info() << std::endl;
	
	for (My_Int l3=0; l3<N_eoms; l3++) {
		for (My_Int i1=1; i1 < Np; i1++) {
			for (My_Int i2=0; i2<Nx; i2++) {
				
				f[l3](i1,i2)=f[l3](i1,i2)-X[i2+(i1-1)*Nx+l3*(Np-1)*Nx];
			}}}
    
	for (My_Int l3=0; l3<N_bconds-N_eoms;++l3){
		for(My_Int i2=0; i2<Nx;++i2){
			
			f[l3](0,i2)=f[l3](0,i2)-X[N_eoms*(Np-1)*Nx+l3*Nx+i2];
		}
	}
	
	
	for (My_Int l3=0; l3<N_functions; ++l3) {
		f[l3].update(D);
	}
	
	std::cout << "Called the <float128, float128> function!" << std::endl;
	
	Test=B-MSp*X;
	
	float128 res(0);
	
	for(My_Int ii=0;ii<nbdof;++ii){
		res=(fabs(Test[ii])>res)?fabs(Test[ii]):res;
	}
	
	std::cout << "Residue estimate: " << res << std::endl;

	
}

template<>
void UpdateFunctions<float128,ComplexQ>( void (*eom) (ComplexQ* elem, Cfunction<float128, ComplexQ> *F1 , grid<float128>& grid, float128* params, My_Int i, My_Int j)
										,void (*b_cond) (ComplexQ* elem, Cfunction<float128, ComplexQ> *F1,float128* params, My_Int i, My_Int j)
										,void (*l_eom)(ComplexQ **elem,  derivativeCol<float128>& Der,Cfunction<float128, ComplexQ> *F1, grid<float128>& grid, float128* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
										,void (*l_b_cond)(ComplexQ **elem,  derivativeCol<float128>& Der,Cfunction<float128, ComplexQ> *F1, float128* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<float128, ComplexQ>* f, derivativeCol<float128>& D, grid<float128>& Rgrid, grid<float128>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds,float128* params){
	
	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
	My_Int nbdof = (Np-1)*Nx*N_eoms+(N_bconds-N_eoms)*Nx; // number of degrees of freedom.
	My_Int nvars = (Rgrid.MaxOrder()*N_eoms+(N_bconds-N_eoms))*Xgrid.MaxOrder();
	
	
	std::vector<Eigen::Triplet<ComplexQ, My_Int> > TriplVec;
	Eigen::SparseMatrix<ComplexQ, Eigen::ColMajor, My_Int> MSp(nbdof,nbdof);
	Eigen::Matrix<ComplexQ, Eigen::Dynamic, 1> B(nbdof);
	Eigen::Matrix<ComplexQ, Eigen::Dynamic, 1> X(nbdof);
	
	TriplVec.reserve(nvars*nbdof);
    
    float128 tol=5.E-30;
	
	ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid,N_functions,N_eoms,N_bconds,params, tol);
	
	std::cout << "Reserved: " << nvars*nbdof << ", Needed: " << TriplVec.size() << std::endl;
	
	std::cout << "Starting solver with dimension "<< nbdof << " x " << nbdof << " matrices and with " <<TriplVec.size() << " non-zero elements"<< std::endl;
	
	MSp.setFromTriplets(TriplVec.begin(),TriplVec.end());
	MSp.makeCompressed();
	
	TriplVec.erase(TriplVec.begin(), TriplVec.end());
    
	
	Eigen::SparseLU<Eigen::SparseMatrix<ComplexQ, Eigen::ColMajor, My_Int>,Eigen::COLAMDOrdering<My_Int> > solver;
	//	solver.analyzePattern(MSp);
	//	solver.factorize(MSp);
	solver.compute(MSp);
	std::cout<< solver.info() << std::endl;
	X=solver.solve(B);
	std::cout<< solver.info() << std::endl;
	
	for (My_Int l3=0; l3<N_eoms; l3++) {
		for (My_Int i1=1; i1 < Np; i1++) {
			for (My_Int i2=0; i2<Nx; i2++) {
				
				f[l3](i1,i2)=f[l3](i1,i2)-X[i2+(i1-1)*Nx+l3*(Np-1)*Nx];
			}}}
    
	for (My_Int l3=0; l3<N_bconds-N_eoms;++l3){
		for(My_Int i2=0; i2<Nx;++i2){
			
			f[l3](0,i2)=f[l3](0,i2)-X[N_eoms*(Np-1)*Nx+l3*Nx+i2];
		}
	}
	
	
	for (My_Int l3=0; l3<N_functions; ++l3) {
		f[l3].update(D);
	}
	
	std::cout << "Called the <float128, complexQ> function!" << std::endl;
}

template<>
void UpdateFunctions<mpreal,mpreal>( void (*eom) (mpreal* elem, Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid, mpreal* params, My_Int i, My_Int j)
									,void (*b_cond) (mpreal* elem, Cfunction<mpreal, mpreal> *F1,mpreal* params, My_Int i, My_Int j)
									,void (*l_eom)(mpreal **elem, derivativeCol<mpreal>& Der,Cfunction<mpreal, mpreal> *F1, grid<mpreal>& grid, mpreal* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
									,void (*l_b_cond)(mpreal **elem,  derivativeCol<mpreal>& Der,Cfunction<mpreal, mpreal> *F1, mpreal* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<mpreal, mpreal>* f, derivativeCol<mpreal>& D, grid<mpreal>& Rgrid,grid<mpreal>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds, mpreal* params){
	
	
		
	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
	My_Int nbdof = (Np-1)*Nx*N_eoms+(N_bconds-N_eoms)*Nx; // number of degrees of freedom.
	My_Int nvars = (Rgrid.MaxOrder()*N_eoms+(N_bconds-N_eoms))*Xgrid.MaxOrder();
	
	std::vector<Eigen::Triplet<mpreal, My_Int> > TriplVec;
	Eigen::SparseMatrix<mpreal, Eigen::ColMajor, My_Int> MSp(nbdof,nbdof);
	Eigen::Matrix<mpreal, Eigen::Dynamic, 1> B(nbdof);
	Eigen::Matrix<mpreal, Eigen::Dynamic, 1> X(nbdof);
	Eigen::Matrix<mpreal, Eigen::Dynamic, 1> Test(nbdof);
	
	TriplVec.reserve(nvars*nbdof);
	
	ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid, N_functions,N_eoms,N_bconds,params, mpreal(0));
	
	std::cout << "Reserved: " << nvars*nbdof << ", Needed: " << TriplVec.size() << std::endl;
	
	std::cout << "Starting solver with dimension "<< nbdof << " x " << nbdof << " matrices and with " <<TriplVec.size() << " non-zero elements"<< std::endl;
	
	MSp.setFromTriplets(TriplVec.begin(),TriplVec.end());
	
	TriplVec.erase(TriplVec.begin(), TriplVec.end());
	
	MSp.makeCompressed();
	
	//	Eigen::UmfPackLU<Eigen::SparseMatrix<mpreal> > solver;
	//	Eigen::BiCGSTAB<Eigen::SparseMatrix<mpreal>,Eigen::IncompleteLUT<mpreal>> solver(MSp);
	Eigen::SparseLU<Eigen::SparseMatrix<mpreal, Eigen::ColMajor, My_Int>,Eigen::COLAMDOrdering<My_Int> > solver(MSp);
	//	solver.compute(MSp);
	std::cout<< solver.info() << std::endl;
	X=solver.solve(B);
	std::cout<< solver.info() << std::endl;
	
	for (My_Int l3=0; l3<N_eoms; l3++) {
		for (My_Int i1=1; i1 < Np; i1++) {
			for (My_Int i2=0; i2<Nx; i2++) {
				
				f[l3](i1,i2)=f[l3](i1,i2)-X[i2+(i1-1)*Nx+l3*(Np-1)*Nx];
			}}}
    
	for (My_Int l3=0; l3<N_bconds-N_eoms;++l3){
		for(My_Int i2=0; i2<Nx;++i2){
			
			f[l3](0,i2)=f[l3](0,i2)-X[N_eoms*(Np-1)*Nx+l3*Nx+i2];
		}
	}
	
	
	for (My_Int l3=0; l3<N_functions; ++l3) {
		f[l3].update(D);
	}
	
	std::cout << "Called the <mpreal, mpreal> function!" << std::endl;
	
	Test=B-MSp*X;
	
	mpreal res(0);
	
	for(My_Int ii=0;ii<nbdof;++ii){
		res=(fabs(Test[ii])>res)?fabs(Test[ii]):res;
	}
	
	std::cout << "Residue estimate: " << res << std::endl;
	
}

template<>
void UpdateFunctions<mpreal,ComplexMP>( void (*eom) (ComplexMP* elem, Cfunction<mpreal, ComplexMP> *F1 , grid<mpreal>& grid, mpreal* params, My_Int i, My_Int j)
									   ,void (*b_cond) (ComplexMP* elem, Cfunction<mpreal, ComplexMP> *F1,mpreal* params, My_Int i, My_Int j)
									   ,void (*l_eom)(ComplexMP **elem,  derivativeCol<mpreal>& Der,Cfunction<mpreal, ComplexMP> *F1, grid<mpreal>& grid, mpreal* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
									   ,void (*l_b_cond)(ComplexMP **elem,  derivativeCol<mpreal>& Der,Cfunction<mpreal, ComplexMP> *F1, mpreal* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<mpreal, ComplexMP>* f, derivativeCol<mpreal>& D, grid<mpreal>& Rgrid, grid<mpreal>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds,mpreal* params){
	
	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
	My_Int nbdof = (Np-1)*Nx*N_eoms+(N_bconds-N_eoms)*Nx; // number of degrees of freedom.
	My_Int nvars = (Rgrid.MaxOrder()*N_eoms+(N_bconds-N_eoms))*Xgrid.MaxOrder();
	
	
	std::vector<Eigen::Triplet<ComplexMP, My_Int> > TriplVec;
	Eigen::SparseMatrix<ComplexMP, Eigen::ColMajor, My_Int> MSp(nbdof,nbdof);
	Eigen::Matrix<ComplexMP, Eigen::Dynamic, 1> B(nbdof);
	Eigen::Matrix<ComplexMP, Eigen::Dynamic, 1> X(nbdof);
	
	TriplVec.reserve(nvars*nbdof);
	
	ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid,N_functions,N_eoms,N_bconds,params,mpreal(0));
	
	std::cout << "Reserved: " << nvars*nbdof << ", Needed: " << TriplVec.size() << std::endl;
	
	std::cout << "Starting solver with dimension "<< nbdof << " x " << nbdof << " matrices and with " <<TriplVec.size() << " non-zero elements"<< std::endl;
	
	MSp.setFromTriplets(TriplVec.begin(),TriplVec.end());
	MSp.makeCompressed();
	
	TriplVec.erase(TriplVec.begin(), TriplVec.end());
	
	Eigen::SparseLU<Eigen::SparseMatrix<ComplexMP, Eigen::ColMajor, My_Int>, Eigen::COLAMDOrdering<My_Int> > solver(MSp);
	//	solver.analyzePattern(MSp);
	//	solver.factorize(MSp);
	solver.compute(MSp);
	std::cout<< solver.info() << std::endl;
	X=solver.solve(B);
	std::cout<< solver.info() << std::endl;
	
	//	X=scal.RightScaling().cwiseProduct(X);
	
	for (My_Int l3=0; l3<N_eoms; l3++) {
		for (My_Int i1=1; i1 < Np; i1++) {
			for (My_Int i2=0; i2<Nx; i2++) {
				
				f[l3](i1,i2)=f[l3](i1,i2)-X[i2+(i1-1)*Nx+l3*(Np-1)*Nx];
			}}}
    
	for (My_Int l3=0; l3<N_bconds-N_eoms;++l3){
		for(My_Int i2=0; i2<Nx;++i2){
			
			f[l3](0,i2)=f[l3](0,i2)-X[N_eoms*(Np-1)*Nx+l3*Nx+i2];
		}
	}
	
	
	for (My_Int l3=0; l3<N_functions; ++l3) {
		f[l3].update(D);
	}
	
	std::cout << "Called the <mpreal, ComplexMP> function!" << std::endl;
}

template<>
void UpdateFunctions<long double,long double>( void (*eom) (long double* elem, Cfunction<long double, long double> *F1, grid<long double>& grid, long double* params, My_Int i, My_Int j)
											  ,void (*b_cond) (long double* elem, Cfunction<long double, long double> *F1,long double* params, My_Int i, My_Int j)
											  ,void (*l_eom)(long double **elem, derivativeCol<long double>& Der,Cfunction<long double, long double> *F1, grid<long double>& grid, long double* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
											  ,void (*l_b_cond)(long double **elem,  derivativeCol<long double>& Der,Cfunction<long double, long double> *F1, long double* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<long double, long double>* f, derivativeCol<long double>& D, grid<long double>& Rgrid,grid<long double>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds, long double* params){
	
	
	
	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
	My_Int nbdof = (Np-1)*Nx*N_eoms+(N_bconds-N_eoms)*Nx; // number of degrees of freedom.
	My_Int nvars = (Rgrid.MaxOrder()*N_eoms+(N_bconds-N_eoms))*Xgrid.MaxOrder();
	
	std::vector<Eigen::Triplet<long double, My_Int> > TriplVec;
	Eigen::SparseMatrix<long double, Eigen::ColMajor, My_Int> MSp(nbdof,nbdof);
	Eigen::Matrix<long double, Eigen::Dynamic, 1> B(nbdof);
	Eigen::Matrix<long double, Eigen::Dynamic, 1> X(nbdof);
	Eigen::Matrix<long double, Eigen::Dynamic, 1> Test(nbdof);
	
	TriplVec.reserve(nvars*nbdof);
	
	ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid, N_functions,N_eoms,N_bconds,params, 1.E-15L);
	
	std::cout << "Reserved: " << nvars*nbdof << ", Needed: " << TriplVec.size() << std::endl;
	
	std::cout << "Starting solver with dimension "<< nbdof << " x " << nbdof << " matrices and with " <<TriplVec.size() << " non-zero elements"<< std::endl;
	
	MSp.setFromTriplets(TriplVec.begin(),TriplVec.end());
	
	TriplVec.erase(TriplVec.begin(), TriplVec.end());
	
	MSp.makeCompressed();
	
	//	Eigen::UmfPackLU<Eigen::SparseMatrix<long double> > solver;
	//	Eigen::BiCGSTAB<Eigen::SparseMatrix<long double>,Eigen::IncompleteLUT<long double>> solver(MSp);
	Eigen::SparseLU<Eigen::SparseMatrix<long double, Eigen::ColMajor, My_Int>,Eigen::COLAMDOrdering<My_Int> > solver(MSp);
	//	solver.compute(MSp);
	std::cout<< solver.info() << std::endl;
	X=solver.solve(B);
	std::cout<< solver.info() << std::endl;
	
	for (My_Int l3=0; l3<N_eoms; l3++) {
		for (My_Int i1=1; i1 < Np; i1++) {
			for (My_Int i2=0; i2<Nx; i2++) {
				
				f[l3](i1,i2)=f[l3](i1,i2)-X[i2+(i1-1)*Nx+l3*(Np-1)*Nx];
			}}}
    
	for (My_Int l3=0; l3<N_bconds-N_eoms;++l3){
		for(My_Int i2=0; i2<Nx;++i2){
			
			f[l3](0,i2)=f[l3](0,i2)-X[N_eoms*(Np-1)*Nx+l3*Nx+i2];
		}
	}
	
	
	for (My_Int l3=0; l3<N_functions; ++l3) {
		f[l3].update(D);
	}
	
	std::cout << "Called the <long double, long double> function!" << std::endl;
	
	Test=B-MSp*X;
	
	long double res(0);
	
	for(My_Int ii=0;ii<nbdof;++ii){
		res=(fabs(Test[ii])>res)?fabs(Test[ii]):res;
	}
	
	std::cout << "Residue estimate: " << res << std::endl;
	
}

template<>
void UpdateFunctions<long double,ComplexLD>( void (*eom) (ComplexLD* elem, Cfunction<long double, ComplexLD> *F1 , grid<long double>& grid, long double* params, My_Int i, My_Int j)
											,void (*b_cond) (ComplexLD* elem, Cfunction<long double, ComplexLD> *F1,long double* params, My_Int i, My_Int j)
											,void (*l_eom)(ComplexLD **elem,  derivativeCol<long double>& Der,Cfunction<long double, ComplexLD> *F1, grid<long double>& grid, long double* params, My_Int i1, My_Int j1, My_Int i, My_Int j)
											,void (*l_b_cond)(ComplexLD **elem,  derivativeCol<long double>& Der,Cfunction<long double, ComplexLD> *F1, long double* params, My_Int i1, My_Int j1, My_Int i, My_Int j),Cfunction<long double, ComplexLD>* f, derivativeCol<long double>& D, grid<long double>& Rgrid, grid<long double>& Xgrid,My_Int N_functions,My_Int N_eoms, My_Int N_bconds,long double* params){
	
	My_Int Np=Rgrid.TotalLength();
	My_Int Nx=Xgrid.TotalLength();
	My_Int nbdof = (Np-1)*Nx*N_eoms+(N_bconds-N_eoms)*Nx; // number of degrees of freedom.
	My_Int nvars = (Rgrid.MaxOrder()*N_eoms+(N_bconds-N_eoms))*Xgrid.MaxOrder();
	
	
	std::vector<Eigen::Triplet<ComplexLD, My_Int> > TriplVec;
	Eigen::SparseMatrix<ComplexLD, Eigen::ColMajor, My_Int> MSp(nbdof,nbdof);
	Eigen::Matrix<ComplexLD, Eigen::Dynamic, 1> B(nbdof);
	Eigen::Matrix<ComplexLD, Eigen::Dynamic, 1> X(nbdof);
	
	TriplVec.reserve(nvars*nbdof);
	
	ConstructLinearOp(TriplVec, B, eom, b_cond, l_eom, l_b_cond,f,D,Rgrid, Xgrid,N_functions,N_eoms,N_bconds,params, 1.E-15L);
	
	std::cout << "Reserved: " << nvars*nbdof << ", Needed: " << TriplVec.size() << std::endl;
	
	std::cout << "Starting solver with dimension "<< nbdof << " x " << nbdof << " matrices and with " <<TriplVec.size() << " non-zero elements"<< std::endl;
	
	MSp.setFromTriplets(TriplVec.begin(),TriplVec.end());
	MSp.makeCompressed();
	
	TriplVec.erase(TriplVec.begin(), TriplVec.end());
	
	Eigen::SparseLU<Eigen::SparseMatrix<ComplexLD, Eigen::ColMajor, My_Int>, Eigen::COLAMDOrdering<My_Int> > solver(MSp);
	//	solver.analyzePattern(MSp);
	//	solver.factorize(MSp);
	solver.compute(MSp);
	std::cout<< solver.info() << std::endl;
	X=solver.solve(B);
	std::cout<< solver.info() << std::endl;
	
	
	for (My_Int l3=0; l3<N_eoms; l3++) {
		for (My_Int i1=1; i1 < Np; i1++) {
			for (My_Int i2=0; i2<Nx; i2++) {
				
				f[l3](i1,i2)=f[l3](i1,i2)-X[i2+(i1-1)*Nx+l3*(Np-1)*Nx];
			}}}
    
	for (My_Int l3=0; l3<N_bconds-N_eoms;++l3){
		for(My_Int i2=0; i2<Nx;++i2){
			
			f[l3](0,i2)=f[l3](0,i2)-X[N_eoms*(Np-1)*Nx+l3*Nx+i2];
		}
	}
	
	
	for (My_Int l3=0; l3<N_functions; ++l3) {
		f[l3].update(D);
	}
	
	std::cout << "Called the <long double, ComplexLD> function!" << std::endl;
}