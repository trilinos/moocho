// ///////////////////////////////////////////////////////////////////////////////////////
// BLAS_Cpp.cpp

// ////////////////////////////////////////////////////////////
// When using the Intel BLAS you have to use lower case names

#ifdef _WINDOWS
#define _INTEL_BLAS
#endif

#ifdef _INTEL_BLAS
#define FORTRAN_NAMES_UPPERCASE 0
#endif

#include <stdexcept>

#include "../include/BLAS_Cpp.h"

// /////////////////////////////////////
// Fortran function declarations.

namespace BLAS_C_Decl {

using namespace std;
typedef FortranTypes::f_int  		f_int;
typedef FortranTypes::f_real		f_real;
typedef FortranTypes::f_dbl_prec	f_dbl_prec;

// function declarations 
// (don't use these directly, use the later overloaded wrappers in namespace BLAS)
extern "C" {

// ////////////////////////////////////////
// Level 1 BLAS (vector-vector operations)

// Generate plane rotation
FORTRAN_FUNC_DECL_UL(void,DROTG,drotg)(f_dbl_prec& A, f_dbl_prec& B, f_dbl_prec& C, f_dbl_prec& S);

// Generate modified plane rotation
FORTRAN_FUNC_DECL_UL(void,DROTMG,drotmg)(f_dbl_prec& D1, f_dbl_prec& D2, f_dbl_prec& A, const f_dbl_prec& B, f_dbl_prec* PARAM);
 
// Apply plane rotation
FORTRAN_FUNC_DECL_UL(void,DROT,drot)(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY,
		  const f_dbl_prec& C, const f_dbl_prec& S);

// Apply modified plane rotation
FORTRAN_FUNC_DECL_UL(void,DROTM,drotm)(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY,
		  const f_dbl_prec* PARAM);

// Interchange vectors
FORTRAN_FUNC_DECL_UL(void,DSWAP,dswap)(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY);

// Vector scaling
FORTRAN_FUNC_DECL_UL(void,DSCAL,dscal)(const f_int& N, const f_dbl_prec& ALPHA, f_dbl_prec* X, const f_int& INCX);

// Vector copy 
FORTRAN_FUNC_DECL_UL(void,DCOPY,dcopy)(const f_int& N, const f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY);

// y = a*x + y
FORTRAN_FUNC_DECL_UL(void,DAXPY,daxpy)(const f_int& N, const f_dbl_prec& A, const f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y,
			 const f_int& INCY);

// Dot product
FORTRAN_FUNC_DECL_UL(f_dbl_prec,DDOT,ddot)(const f_int& N, const f_dbl_prec* X, const f_int& INCX, const f_dbl_prec* Y, const f_int& INCY);
FORTRAN_FUNC_DECL_UL(f_dbl_prec,DSDOT,dsdot)(const f_int& N, const f_real* X, const f_int& INCX, const f_real* Y, const f_int& INCY);

// 2-Norm
FORTRAN_FUNC_DECL_UL(f_dbl_prec,DNRM2,dnrm2)(const f_int& N, const f_dbl_prec* X, const f_int& INCX);

// Sum of magnitudes
FORTRAN_FUNC_DECL_UL(f_dbl_prec,DASUM,dasum)(const f_int& N, const f_dbl_prec* X, const f_int& INCX);

// Largest component of vector
FORTRAN_FUNC_DECL_UL(f_dbl_prec,IDAMAX,idamax)(const f_int& N, const f_dbl_prec* X, const f_int& INCX);

// ////////////////////////////////////////
// Level 2 BLAS (matrix-vector operations)

// General rectangular matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DGEMV,dgemv)(const char* TRANSA, const f_int& M, const f_int& N, const f_dbl_prec& ALPHA,
			const f_dbl_prec* A, const f_int& LDA, const f_dbl_prec* X, const f_int& INCX,
			const f_dbl_prec& BETA, f_dbl_prec* Y, const f_int& INCY);

// General band matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DGBMV,dgbmv)(const char* TRANSA, const f_int& M, const f_int& N, const f_int& KL, const f_int& KU,
		   const f_dbl_prec& ALPHA,	const f_dbl_prec* A, const f_int& LDA, const f_dbl_prec* X,
		   const f_int& INCX,	const f_dbl_prec& BETA, f_dbl_prec* Y, const f_int& INCY);

// Hermitian matrix-vector products

// Hermitian band matrix-vector products

// Hermitian packed matrix-vector products

// Symmetric matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DSYMV,dsymv)(const char* UPLO, const f_int& N,
		   const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA,
		   const f_dbl_prec* X, const f_int& INCX, const f_dbl_prec& BETA,
		   f_dbl_prec* Y, const f_int& INCY);

// Symmetric band matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DSBMV,dsbmv)(const char* UPLO, const f_int& N, const f_int& K,
		   const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA,
		   const f_dbl_prec* X, const f_int& INCX, const f_dbl_prec& BETA,
		   f_dbl_prec* Y, const f_int& INCY);

// Symmetric packed matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DSPMV,dspmv)(const char* UPLO, const f_int& N,
		   const f_dbl_prec& ALPHA, const f_dbl_prec* AP,
		   const f_dbl_prec* X, const f_int& INCX, const f_dbl_prec& BETA,
		   f_dbl_prec* Y, const f_int& INCY);

// Triangular matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DTRMV,dtrmv)(const char* UPLO, const char* TRANS, const char* DIAG, const f_int& N,
		   const f_dbl_prec* A, const f_int& LDA, f_dbl_prec* X, const f_int& INCX);

// Triangular band matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DTBMV,dtbmv)(const char* UPLO, const char* TRANS, const char* DIAG, const f_int& N, const f_int& K,
		   const f_dbl_prec* A, const f_int& LDA, f_dbl_prec* X, const f_int& INCX);

// Triangular packed matrix-vector products
FORTRAN_FUNC_DECL_UL(void,DTPMV,dtpmv)(const char* UPLO, const char* TRANS, const char* DIAG, const f_int& N,
		   const f_dbl_prec* AP, f_dbl_prec* X, const f_int& INCX);

// Triangular equation solve
FORTRAN_FUNC_DECL_UL(void,DTRSV,dtrsv)(const char* UPLO, const char* TRANS, const char* DIAG, const f_int& N,
		   const f_dbl_prec* A, const f_int& LDA, f_dbl_prec* X, const f_int& INCX);

// Triangular band equation solve
FORTRAN_FUNC_DECL_UL(void,DTBSV,dtbsv)(const char* UPLO, const char* TRANS, const char* DIAG, const f_int& N, const f_int& K,
		   const f_dbl_prec* A, const f_int& LDA, f_dbl_prec* X, const f_int& INCX);

// Triangular packed equation solve
FORTRAN_FUNC_DECL_UL(void,DTPSV,dtpsv)(const char* UPLO, const char* TRANS, const char* DIAG, const f_int& N,
		   const f_dbl_prec* AP, f_dbl_prec* X, const f_int& INCX);

// General rank-1 update
FORTRAN_FUNC_DECL_UL(void,DGER,dger)(const f_int& M, const f_int& N, const f_dbl_prec& ALPHA, const f_dbl_prec* X, const f_int& INCX,
		  const f_dbl_prec* Y, const f_int& INCY, f_dbl_prec* A, const f_int& LDA);

// Hermitian rank-1 update

// Hermitian packed rank-1 update

// Hermitian rank-2 update

// Hermitian packed rank-2 update

// Symmetric rank-1 update
FORTRAN_FUNC_DECL_UL(void,DSYR,dsyr)(const char* UPLO, const f_int& N, const f_dbl_prec& ALPHA, const f_dbl_prec* X, const f_int& INCX,
		  f_dbl_prec* A, const f_int& LDA);

// Symmetric packed rank-1 update
FORTRAN_FUNC_DECL_UL(void,DSPR,dspr)(const char* UPLO, const f_int& N, const f_dbl_prec& ALPHA, const f_dbl_prec* X, const f_int& INCX,
		  f_dbl_prec* AP);

// Symmetric rank-2 update
FORTRAN_FUNC_DECL_UL(void,DSYR2,dsyr2)(const char* UPLO, const f_int& N, const f_dbl_prec& ALPHA, const f_dbl_prec* X, const f_int& INCX,
		   const f_dbl_prec* Y, const f_int& INCY, f_dbl_prec* A, const f_int& LDA);

// Symmetric packed rank-2 update
FORTRAN_FUNC_DECL_UL(void,DSPR2,dspr2)(const char* UPLO, const f_int& N, const f_dbl_prec& ALPHA, const f_dbl_prec* X, const f_int& INCX,
		   const f_dbl_prec* Y, const f_int& INCY, f_dbl_prec* AP);

// ////////////////////////////////////////
// Level 3 BLAS (matrix-matrix operations)

// General rectangular matrix-matrix product
FORTRAN_FUNC_DECL_UL(void,DGEMM,dgemm)(const char* TRANSA, const char* TRANSB, const f_int& M, const f_int& N, const f_int& K,  
		const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA, const f_dbl_prec* B, const f_int& LDB,
		const f_dbl_prec& BETA, f_dbl_prec* C, const f_int& LDC);

// Symmetric matrix-matrix product
FORTRAN_FUNC_DECL_UL(void,DSYMM,dsymm)(const char* SIDE, const char* UPLO, const f_int& M, const f_int& N,  
		const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA, const f_dbl_prec* B, const f_int& LDB,
		const f_dbl_prec& BETA, f_dbl_prec* C, const f_int& LDC);

// Hermitian matrix-matrix product

// Symmetric rank-k update
FORTRAN_FUNC_DECL_UL(void,DSYRK,dsyrk)(const char* UPLO, const char* TRANS, const f_int& N, const f_int& K,  
		const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA,
		const f_dbl_prec& BETA, f_dbl_prec* C, const f_int& LDC);

// Hermitian rank-k update

// Symmetric rank-2k update
FORTRAN_FUNC_DECL_UL(void,DSYR2K,dsyr2k)(const char* UPLO, const char* TRANS, const f_int& N, const f_int& K,  
		const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA, const f_dbl_prec* B, const f_int& LDB,
		const f_dbl_prec& BETA, f_dbl_prec* C, const f_int& LDC);

// Hermitian rank-2k update

// Triangular matrix-matrix product
FORTRAN_FUNC_DECL_UL(void,DTRMM,dtrmm)(const char* SIDE, const char* UPLO, const char* TRANSA, const char* DIAG,
			const f_int& M, const f_int& N,	const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA,
			f_dbl_prec* B, const f_int& LDB);

// Solution of triangular system 
FORTRAN_FUNC_DECL_UL(void,DTRSM,dtrsm)(const char* SIDE, const char* UPLO, const char* TRANSA, const char* DIAG,
			const f_int& M, const f_int& N,	const f_dbl_prec& ALPHA, const f_dbl_prec* A, const f_int& LDA,
			f_dbl_prec* B, const f_int& LDB);

} // end extern "C"

} // end namespace BLAS_C_Decl

// ////////////////////////////////////////////////////////
// C++ BLAS Function Declarations

// Level 1 BLAS (vector-vector operations)

void BLAS_Cpp::rotg(f_dbl_prec& a, f_dbl_prec& b, f_dbl_prec& c, f_dbl_prec& s) {
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DROTG,drotg)(a,b,c,s);
}

// Generate modified plane rotation

void BLAS_Cpp::rotmg(f_dbl_prec& d1, f_dbl_prec& d2, f_dbl_prec& a, const f_dbl_prec& b, f_dbl_prec* param) {
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DROTMG,drotmg)(d1, d2, a, b, param);
}
 
// Apply plane rotation

void BLAS_Cpp::rot(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY
	, const f_dbl_prec& C, const f_dbl_prec& S)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DROT,drot)(N, X, INCX, Y, INCY, C, S);
}


// Apply modified plane rotation

void BLAS_Cpp::rot(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY
	, const f_dbl_prec* PARAM)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DROTM,drotm)(N, X, INCX, Y, INCY, PARAM);
}

// Interchange vectors

void BLAS_Cpp::swap(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSWAP,dswap)(N, X, INCX, Y, INCY);
}		 	

// Vector scaling

void BLAS_Cpp::scal(const f_int& N, const f_dbl_prec& ALPHA, f_dbl_prec* X, const f_int& INCX)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSCAL,dscal)(N, ALPHA, X, INCX);
}

// Vector copy

void BLAS_Cpp::copy(const f_int& N, const f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DCOPY,dcopy)(N, X, INCX, Y, INCY);
}

// y = a*x + y

void BLAS_Cpp::axpy(const f_int& N, const f_dbl_prec& A, const f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y
	, const f_int& INCY)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DAXPY,daxpy)(N, A, X, INCX, Y, INCY);
}

// Dot product

BLAS_Cpp::f_dbl_prec BLAS_Cpp::dot(const f_int& N, const f_dbl_prec* X, const f_int& INCX, const f_dbl_prec* Y, const f_int& INCY)
{
	return BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DDOT,ddot)(N, X, INCX, Y, INCY);
}

// 2-Norm

BLAS_Cpp::f_dbl_prec BLAS_Cpp::nrm2(const f_int& N, const f_dbl_prec* X, const f_int& INCX)
{
	return BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DNRM2,dnrm2)(N, X, INCX);
}

// 1-Norm

BLAS_Cpp::f_dbl_prec BLAS_Cpp::asum(const f_int& N, const f_dbl_prec* X, const f_int& INCX)
{
	return BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DASUM,dasum)(N, X, INCX);
}

// Inifinity-Norm

BLAS_Cpp::f_dbl_prec BLAS_Cpp::iamax(const f_int& N, const f_dbl_prec* X, const f_int& INCX)
{
	return BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(IDAMAX,idamax)(N, X, INCX);
}

// Level-2 BLAS (matrix-vector operations)

// General rectangular matrix-vector products

void BLAS_Cpp::gemv(Transp transa, f_int m, f_int n, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DGEMV,dgemv)(&TransChar[transa], m, n, alpha, pa, lda, x, incx, beta, py, incy);
}

// General band matrix-vector products

void BLAS_Cpp::gbmv(Transp transa, f_int m, f_int n, f_int kl, f_int ku, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DGBMV,dgbmv)(&TransChar[transa], m, n, kl, ku, alpha, pa, lda, x, incx, beta, py, incy);
}
			 	
// Hermitian matrix-vector products

// Hermitian band matrix-vector products

// Hermitian packed matrix-vector products

// Symmetric matrix-vector products

void BLAS_Cpp::symv(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSYMV,dsymv)(&UploChar[uplo], n, alpha, pa, lda, x, incx, beta, py, incy);
}
			
// Symmetric band matrix-vector products

void BLAS_Cpp::sbmv(Uplo uplo, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSBMV,dsbmv)(&UploChar[uplo], n, k, alpha, pa, lda, x, incx, beta, py, incy);
}

// Symmetric packed matrix-vector products

void BLAS_Cpp::spmv(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* pap
	, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSPMV,dspmv)(&UploChar[uplo], n, alpha, pap, x, incx, beta, py, incy);
}

// Triangular matrix-vector products

void BLAS_Cpp::trmv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pa
	, f_int lda, f_dbl_prec* px, f_int incx)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTRMV,dtrmv)(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pa, lda, px, incx);
}

// Triangular band matrix-vector products

void BLAS_Cpp::tbmv(Uplo uplo, Transp trans, Diag diag, f_int n, f_int k, const f_dbl_prec* pa
	, f_int lda, f_dbl_prec* px, f_int incx)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTBMV,dtbmv)(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, k, pa, lda, px, incx);
}

// Triangular packed matrix-vector products

void BLAS_Cpp::tpmv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pap
	, f_dbl_prec* px, f_int incx)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTPMV,dtpmv)(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pap, px, incx);
}

// Triangular equation solve

void BLAS_Cpp::trsv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pa
	, f_int lda, f_dbl_prec* px, f_int incx)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTRSV,dtrsv)(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pa, lda, px, incx);
}

// Triangular band equation solve

void BLAS_Cpp::tbsv(Uplo uplo, Transp trans, Diag diag, f_int n, f_int k, const f_dbl_prec* pa
	, f_int lda, f_dbl_prec* px, f_int incx)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTBSV,dtbsv)(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, k, pa, lda, px, incx);
}

// Triangular packed equation solve

void BLAS_Cpp::tpsv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pap
	, f_dbl_prec* px, f_int incx)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTPSV,dtpsv)(&UploChar[uplo], &TransChar[trans], &DiagChar[diag], n, pap, px, incx);
}

// General rank-1 update

void BLAS_Cpp::ger(f_int m, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
	, f_int incx, const f_dbl_prec* py, f_int incy, f_dbl_prec* pa, f_int lda)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DGER,dger)(m, n, alpha, px, incx, py, incy, pa, lda);
}

// Hermitian rank-1 update

// Hermitian packed rank-1 update

// Hermitian rank-2 update

// Hermitian packed rank-2 update

// Symmetric rank-1 update

void BLAS_Cpp::syr(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
	, f_int incx, f_dbl_prec* pa, f_int lda)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSYR,dsyr)(&UploChar[uplo], n, alpha, px, incx, pa, lda);
}

// Symmetric packed rank-1 update

void BLAS_Cpp::spr(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
	, f_int incx, f_dbl_prec* pap)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSPR,dspr)(&UploChar[uplo], n, alpha, px, incx, pap);
}

// Symmetric rank-2 update

void BLAS_Cpp::syr2(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
	, f_int incx, const f_dbl_prec* py, f_int incy, f_dbl_prec* pa, f_int lda)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSYR2,dsyr2)(&UploChar[uplo], n, alpha, px, incx, py, incy, pa, lda);
}

// Symmetric packed rank-2 update

void BLAS_Cpp::spr2(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
	, f_int incx, const f_dbl_prec* py, f_int incy, f_dbl_prec* pap)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSPR2,dspr2)(&UploChar[uplo], n, alpha, px, incx, py, incy, pap);
}

// Level 3 BLAS (matrix-matrix operations)	

// General rectangular matrix-matrix product

void BLAS_Cpp::gemm(Transp transa, Transp transb, f_int m, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* pb, f_int ldb, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DGEMM,dgemm)(&TransChar[transa], &TransChar[transb], m, n, k, alpha, pa, lda, pb, ldb
		, beta, pc, ldc);
}

// Symmetric matrix-matrix product

void BLAS_Cpp::symm(Side side, Uplo uplo, f_int m, f_int n, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* pb, f_int ldb, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSYMM,dsymm)(&SideChar[side], &UploChar[uplo], m, n, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}

// Hermitian matrix-matrix product

// Symmetric rank-k update

void BLAS_Cpp::syrk(Uplo uplo, Transp trans, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSYRK,dsyrk)(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, beta, pc, ldc);
}

// Hermitian rank-k update

// Symmetric rank-2k update

void BLAS_Cpp::syr2k(Uplo uplo, Transp trans, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* pb, f_int ldb, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DSYR2K,dsyr2k)(&UploChar[uplo], &TransChar[trans], n, k, alpha, pa, lda, pb, ldb, beta, pc, ldc);
}

// Hermitian rank-2k update

// Triangular matrix-matrix product

void BLAS_Cpp::trmm(Side side, Uplo uplo, Transp transa, Diag diag, f_int m, f_int n, f_dbl_prec alpha
	, const f_dbl_prec* pa, f_int lda, f_dbl_prec* pb, f_int ldb)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTRMM,dtrmm)(&SideChar[side], &UploChar[uplo], &TransChar[transa], &DiagChar[diag]
		, m, n, alpha, pa, lda, pb, ldb);
}

// Solution of triangular system

void BLAS_Cpp::trsm(Side side, Uplo uplo, Transp transa, Diag diag, f_int m, f_int n, f_dbl_prec alpha
	, const f_dbl_prec* pa, f_int lda, f_dbl_prec* pb, f_int ldb)
{
	BLAS_C_Decl::FORTRAN_FUNC_CALL_UL(DTRSM,dtrsm)(&SideChar[side], &UploChar[uplo], &TransChar[transa], &DiagChar[diag]
		, m, n, alpha, pa, lda, pb, ldb);
}