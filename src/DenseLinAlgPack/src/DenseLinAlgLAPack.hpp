// ///////////////////////////////////////////////////////
// DenseLinAlgLAPack.hpp
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

#ifndef LIN_ALG_LA_PACK_H
#define LIN_ALG_LA_PACK_H

#include "DenseLinAlgLAPackTypes.hpp"

namespace DenseLinAlgLAPack {

///
/** Calls xPOTRF to compute the cholesky factorization of a
 * symmetric positive definte matrix.
 *
 * On input A contains the upper or lower trianagular elements of
 * a symmetric positive definite matrix.  On output it contians
 * the cholesky factor U (A = U'*U) if #A.uplo() == upper# and
 * L (A=L*L') if #A.uplo() == lower# of the matrix.
 * 
 * If the matrix is not positive definite a #FactorizationException#
 * will be thrown with an error message.
 */
void potrf( DMatrixSliceTriEle* A );

///
/** Calls xGEQRF to compute the QR factorization of a matrix A.
 *
 * The factors A = Q*R are stored in the storage space for A
 * and in an extra vector tau of length k = min(A.rows(),A.cols()).
 * The matrix R is stored in the upper diagonal of A on
 * exit and Q is stored as a set of elementary reflectors
 * in the elements of A below the diagonal and in the extra
 * vector tau.
 *
 * See LAPACK documentation of xGEQRF for more detail of these
 * parameters.
 *
 * If there is a problem with one of the arguments
 * a std::invalid_argument exception will be thrown.
 *
 * @param	A	 [in/out]	The mxn matrix to be factorized on input.
 *						On output contains R and part of Q.
 *	@param	tau	 [out]	DVector of length min(m,n) for the other part
 *						of Q.
 *	@param	work [work]	Workspace vector of length > N.
 */
void geqrf( DMatrixSlice* A, DVectorSlice* tau, DVectorSlice* work );

///
/** Calls xORMRQ to compute a matrix matrix op(C)=op(Q)*op(C) where
 * Q is stored in A and tau computed from geqrf(..).
 *
 * See the LAPACK documentation for xQRMRQ to see a description
 * of the parameters.
 *
 *	@param	side	[in]	Determines if op(Q) is on the right
 *						(C=op(Q)*C) or on the left (C=C*op(Q)).
 *	@param	trans	[in]	Determines if op(Q) = Q or Q'.
 *	@param	A		[in]	Contains part of Q stored in the lower
 *						diagonal elements as set by geqrf(...).
 *	@param	tau		[in]	Contains the other part of Q as set by
 *						geqrf(...).
 *	@param	C		[in/out] On input contains C and on output contains
 *						  op(Q)*C (left) or C*op(Q) (right).
 *	@param	work	[work] Workspace vector of length > N if side==left
 *						or > M if side==right.
 *	
 */
void ormrq(
	BLAS_Cpp::Side side, BLAS_Cpp::Transp trans
	, const DMatrixSlice& A, const DVectorSlice& tau
	, DMatrixSlice* C, DVectorSlice* work
	);

///
/** Calls xSYTRF to compute the P*A*P' = L'*D*L factorization of a symmetric
 * possibly indefinite matrix.
 *
 * See LAPACK documentation for xSYTRF so see a description of
 * the function.
 *
 * This factorization can be used to solve linear systems by calling
 * sytrs(...).  This function will throw #FactorizationException# if
 * the matrix is singular.
 *
 *	@param	A		[in/out]	On input contains the elements of the symmetric
 *						matrix to be factorized.  On output contains
 *						the factorization.
 *	@param	ipiv	[out]	Array of length >= A.rows().  On output contains
 *						permutation matrix P.
 *	@param	work	[work]	Workspace used to compute factoization.  The length
 *						must be greater than 0 but optimal is n * NB.  NB
 *						is the optimal block size (64?), n = A.rows().
 */
void sytrf( DMatrixSliceTriEle* A, FortranTypes::f_int ipiv[], DVectorSlice* work );

///
/** Calls xSYTRS(...) to compute the solution of the factorized system
 * A * X = B where A was factorized by xSYTRF(...).
 *
 * See LAPACK documentation for xSYTRS so see a description of
 * the function.
 *
 * This factorization can be used to solve linear systems by calling
 * sytrs(...).
 *
 *	@param	A		[in]	Factorization of A computed by sytrf(...).
 *	@param	ipiv	[out]	Array of length >= A.rows() computed by sytrf(...).
 *	@param	B		[in/out] On input contains the rhs matrix B.  On output
 *						contiains the solution matrix X.
 *	@param	work	[work]	Workspace used to compute factoization.  The length
 *						must be greater than 0 but optimal is n * NB.  NB
 *						is the optimal block size (64?), n = A.rows().
 */
void sytrs(
	const DMatrixSliceTriEle& A, FortranTypes::f_int ipiv[]
	,DMatrixSlice* B, DVectorSlice* work
	);

///
/** Calls xGETRF to compute the P'*A = L*U factorization of a general
 * rectuangular matrix.
 *
 * See LAPACK documentation for xGETRF so see a description of
 * the function.
 *
 * This factorization can be used to solve linear systems by calling
 * getrs(...).  This function will throw #FactorizationException# if
 * the matrix is singular.
 *
 * @param  A    [in/out] On input contains the elements of the general
 *              matrix to be factorized.  On output contains
 *              the factorization.
 * @param  ipiv	[out] Array of length >= A.rows().  On output contains
 *              permutation matrix P.
 * @param  rank [out] Returns the rank of the matrix.
 */
void getrf( DMatrixSlice* A, FortranTypes::f_int ipiv[], FortranTypes::f_int* rank );

///
/** Calls xGETRS to solve linear systems with the factors of P'*A = L*U
 * generated by xGETRF.
 *
 * See LAPACK documentation for xGETRS so see a description of
 * the function.
 *
 * This factorization can be used to solve linear systems by calling
 * getrs(...).  This function will throw #FactorizationException# if
 * the matrix is singular.
 *
 * @param  LU    [in] On input contains the elements of the LU factors.
 * @param  ipiv  [out] Array of length >= LU.rows() that contains
 *               permutation matrix P.
 * @param  tranp [in] Determines of (P*L*U) * X = Y (no_trans) is solved or
 *               (P*L*U)' * X = Y (trans) is solved.
 * @param  B     [in/out] On input, contains the rhs matrix Y, on output, contains
 *               the solution matix X.
 */
void getrs(
	const DMatrixSlice& LU, const FortranTypes::f_int ipiv[], BLAS_Cpp::Transp transp
	,DMatrixSlice* B
	);

}	// end namespace DenseLinAlgLAPack

#endif	// LIN_ALG_LA_PACK_H