// ///////////////////////////////////////////////////////////////////////////////////////
// BLAS_Cpp.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.
//
// C++ overloads for BLAS kernals (element type removed from name and enum for operations)

#ifndef BLAS_CPP_OVERLOADS_DECLARATIONS_H
#define BLAS_CPP_OVERLOADS_DECLARATIONS_H

#include "MoochoMoreUtilities/src/fortran_types.hpp"
#include "MoochoMoreUtilities/src/BLAS_CppTypes.hpp"

// Overloaded BLAS wrappers.
// The naming convention is the Fortran BLAS name minus the type prefix.
namespace BLAS_Cpp {

typedef FortranTypes::f_int			f_int;
typedef FortranTypes::f_real		f_real;
typedef FortranTypes::f_dbl_prec	f_dbl_prec;

/* * @name Option Arguments
 * These are emumerations that are used with the overloaded C++ BLAS declarations to replace the
 * error prone use of characters for specifying options
 * @memo enumerations (enum)
 */
// @{

///
const char SideChar[]	= {'L'	, 'R'			};
///
const char TransChar[]	= {'N'	, 'T'	, 'C'	};
///
const char UploChar[]	= {'U'	, 'L'			};
///
const char DiagChar[]	= {'U'	, 'N'			};

// @}

/* * @name C++ BLAS Function Declarations
 * These are overloaded C++ functions that have removed the element type from the name
 * of the BLAS functions and use enumerations for the options arguments.
 */
// @{	

// ///////////////////////////////////////////////////////////////////////////////////////////
/* * @name Level 1 BLAS (vector-vector operations) */
// @{	

/* * @name Generate plane rotation */
// @{

///
void rotg( f_dbl_prec* a, f_dbl_prec* b, f_dbl_prec* c, f_dbl_prec* s );

// @}	 	
 
/* * @name Apply plane rotation */
// @{

///
void rot(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY
		 , const f_dbl_prec& C, const f_dbl_prec& S);
// @}

/* * @name  Interchange vectors */
// @{

///
void swap(const f_int& N, f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY);

// @}

/* * @name  DVector scaling */
// @{

/// 
void scal(const f_int& N, const f_dbl_prec& ALPHA, f_dbl_prec* X, const f_int& INCX);

// @}

/* * @name DVector copy */
// @{

/// 
void copy(const f_int& N, const f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y, const f_int& INCY);

// @}

/* * @name  y = a*x + y */
// @{

///
void axpy(const f_int& N, const f_dbl_prec& A, const f_dbl_prec* X, const f_int& INCX, f_dbl_prec* Y
	, const f_int& INCY);
	
// @}

/* * @name  Dot product */
// @{

///
f_dbl_prec dot(const f_int& N, const f_dbl_prec* X, const f_int& INCX, const f_dbl_prec* Y, const f_int& INCY);

// @}

/* * @name  2-Norm */
// @{

///
f_dbl_prec nrm2(const f_int& N, const f_dbl_prec* X, const f_int& INCX);

// @}

/* * @name  1-Norm */
// @{

///
f_dbl_prec asum(const f_int& N, const f_dbl_prec* X, const f_int& INCX);

// @}

/* * @name  Inifinity-Norm */
// @{

///
f_dbl_prec iamax(const f_int& N, const f_dbl_prec* X, const f_int& INCX);
// @}

//		end Level-1 BLAS
// @}	

// /////////////////////////////////////////////////
/* * @name Level-2 BLAS (matrix-vector operations) */
// @{	

/* * @name General rectangular matrix-vector products */
// @{

///
void gemv(Transp transa, f_int m, f_int n, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy);
			 
// @}	

/* * @name General band matrix-vector products */
// @{

///
void gbmv(Transp transa, f_int m, f_int n, f_int kl, f_int ku, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy);
			 	
// @}

/* * @name Hermitian matrix-vector products */
// @{

			 	
// @}

/* * @name Hermitian band matrix-vector products */
// @{

// @}

/* * @name Hermitian packed matrix-vector products */
// @{


// @}

/* * @name Symmetric matrix-vector products */
// @{

///
void symv(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy);
			
// @}

/* * @name Symmetric band matrix-vector products */
// @{

///
void sbmv(Uplo uplo, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy);

// @}

/* * @name Symmetric packed matrix-vector products */
// @{

///
void spmv(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* pap
	, const f_dbl_prec* x, f_int incx, f_dbl_prec beta, f_dbl_prec* py, f_int incy);

// @}

/* * @name Triangular matrix-vector products */
// @{

///
void trmv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pa
	, f_int lda, f_dbl_prec* px, f_int incx);

// @}

/* * @name Triangular band matrix-vector products */
// @{

///
void tbmv(Uplo uplo, Transp trans, Diag diag, f_int n, f_int k, const f_dbl_prec* pa
	, f_int lda, f_dbl_prec* px, f_int incx);

// @}

/* * @name Triangular packed matrix-vector products */
// @{

///
void tpmv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pap
	, f_dbl_prec* px, f_int incx);

// @}

/* * @name Triangular equation solve */
// @{

///
void trsv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pa
	, f_int lda, f_dbl_prec* px, f_int incx);

// @}

/* * @name Triangular band equation solve */
// @{

///
void tbsv(Uplo uplo, Transp trans, Diag diag, f_int n, f_int k, const f_dbl_prec* pa
	, f_int lda, f_dbl_prec* px, f_int incx);

// @}

/* * @name Triangular packed equation solve */
// @{

///
void tpsv(Uplo uplo, Transp trans, Diag diag, f_int n, const f_dbl_prec* pap
	, f_dbl_prec* px, f_int incx);

// @}

/* * @name General rank-1 update */
// @{

///
void ger(f_int m, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
	, f_int incx, const f_dbl_prec* py, f_int incy, f_dbl_prec* pa, f_int lda);

// @}

/* * @name Hermitian rank-1 update */
// @{

// @}

/* * @name Hermitian packed rank-1 update */
// @{

// @}

/* * @name Hermitian rank-2 update */
// @{

// @}

/* * @name Hermitian packed rank-2 update */
// @{

// @}

/* * @name Symmetric rank-1 update */
// @{

///
void syr(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
	, f_int incx, f_dbl_prec* pa, f_int lda);

// @}

/* * @name Symmetric packed rank-1 update */
// @{

///
void spr(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
	, f_int incx, f_dbl_prec* pap);

// @}

/* * @name Symmetric rank-2 update */
// @{

///
void syr2(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
	, f_int incx, const f_dbl_prec* py, f_int incy, f_dbl_prec* pa, f_int lda);

// @}

/* * @name Symmetric packed rank-2 update */
// @{

///
void spr2(Uplo uplo, f_int n, f_dbl_prec alpha, const f_dbl_prec* px
	, f_int incx, const f_dbl_prec* py, f_int incy, f_dbl_prec* pap);

// @}

//		end Level 2 BLAS
// @}
	
// /////////////////////////////////////////
/* * @name Level 3 BLAS (matrix-matrix operations) */
// @{	

/* * @name General rectangular matrix-matrix product */
// @{

///
void gemm(Transp transa, Transp transb, f_int m, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* pb, f_int ldb, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc);

// @}

/* * @name Symmetric matrix-matrix product */
// @{

///
void symm(Side side, Uplo uplo, f_int m, f_int n, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* pb, f_int ldb, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc);

// @}

/* * @name Hermitian matrix-matrix product */
// @{

// @}

/* * @name Symmetric rank-k update */
// @{

///
void syrk(Uplo uplo, Transp trans, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc);

// @}

/* * @name Hermitian rank-k update */
// @{

// @}

/* * @name Symmetric rank-2k update */
// @{

///
void syr2k(Uplo uplo, Transp trans, f_int n, f_int k, f_dbl_prec alpha, const f_dbl_prec* pa
	, f_int lda, const f_dbl_prec* pb, f_int ldb, f_dbl_prec beta, f_dbl_prec* pc, f_int ldc);

// @}

/* * @name Hermitian rank-2k update */
// @{

// @}

/* * @name Triangular matrix-matrix product */
// @{

///
void trmm(Side side, Uplo uplo, Transp transa, Diag diag, f_int m, f_int n, f_dbl_prec alpha
	, const f_dbl_prec* pa, f_int lda, f_dbl_prec* pb, f_int ldb);

// @}

/* * @name Solution of triangular system */
// @{

///
void trsm(Side side, Uplo uplo, Transp transa, Diag diag, f_int m, f_int n, f_dbl_prec alpha
	, const f_dbl_prec* pa, f_int lda, f_dbl_prec* pb, f_int ldb);
		 	 	
// @}

//		end Level 3 BLAS
// @}	

//		end overloaded functions
// @}

}	// end namespace BLAS_Cpp

#endif // BLAS_CPP_OVERLOADS_DECLARATIONS_H
