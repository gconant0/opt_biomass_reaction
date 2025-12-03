//c++ include file for the library lbprsh_lapack
//Computes the inverse of a matrix--requires LU factored matrix
//from dgetrf
extern "C" int dgetri_( int *n, double *a, int *lda,
		        int *ipiv, double *work, int *lwork,
		       int *info);

//Computes the LU factorization of a matrix
extern "C" int dgetrf_ (int *m, int *n, double *a,  int *lda,
			int *ipiv, int *info);



//Computes the eigen values and vectors of a matrix 
extern "C" int dgeev_(char *jobvl, char *jobvr,  int *n,
		  double *a, int *lda, double *wr,
		  double *wi, double *vl, int *ldvl,
		  double *vr,  int *ldvr, double *work,
		  int *lwork, int *info, int *jobvl_len,
		  int *jobvr_len);

//Performs a matrix/matrix multiplication
extern "C" int dgemm_ (char *transa, char *transb, int *m, int *n,
		        int *k, double *alpha, double *a,
		        int *lda, double *b,  int *ldb,
		       double *beta, double *c__,  int *ldc,
		        int *transa_len, int *transb_len);



//Solves over or underdetermined system
extern "C" int dgels_(char *trans, int *m, int *n, int *nrhs, double *a, int* lda, double *b, int* ldb, double *work, int* lwork, int* info);


//LQ factorization of underdetermined system
extern "C" int dgelqf_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* lwork, int* info);

//Singular-value decomposition
extern "C" void dgesdd_(char *JOBZ,  int *M, int *N, double *A, int *LDA,double *S,          double *U, int *LDU,   double *VT, int *LDVT,  double *WORK, int *LWORK,
          int *IWORK, int *INFO);


extern "C" void dgesvd_(char *JOBU, char *JOBVT, int *M, int *N, double *A, int *LDA,
                        double *S,  double *U, int *LDU,double *VT, int *LDVT,
                        double *WORK, int *LWORK,   int *INFO);
