// ///////////////////////////////////////////////////////
// MatrixSymPosDefBandedChol.h
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

#ifndef MATRIX_SYM_POS_DEF_BANDED_CHOL_H
#define MATRIX_SYM_POS_DEF_BANDED_CHOL_H

#include "ConstrainedOptimizationPackTypes.h"
#include "SparseLinAlgPack/include/MatrixSymWithOpFactorized.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "Misc/include/ref_count_ptr.h"
#include "Misc/include/ReleaseResource.h"

namespace ConstrainedOptimizationPack {
///
/** Matrix subclass for banded symmetric positive definite matrices and their
 * Cholesky factors.
 * 
 * This class is designed to support the LAPACK routines for banded symmetric positive
 * definite matrices.  The banded matrix and/or its cholesky factor are stored in
 * simple flat rectangular matrices compatible with the LAPACK routines.
 *
 * For example, for #n = 8, kd = 3# the original matrix #M# (if set) is stored in the
 * following format #MB#:
 \begin{verbatim}

         M                                   MB

 [ x x x x         ]                  [ x x x x x x x x ]
 [ x x x x x       ]  lower triangle  [ x x x x x x x o ]
 [ x x x x x x     ]         =>       [ x x x x x x o o ]
 [ x x x x x x x   ]                  [ x x x x x o o o ]
 [   x x x x x x x ]
 [     x x x x x x ]                  [ o o o x x x x x ]
 [       x x x x x ]  upper triangle  [ o o x x x x x x ]
 [         x x x x ]         =>       [ o x x x x x x x ]
                                      [ x x x x x x x x ]
 \end{verbatim}
 * The Cholesky factor #U# is sorted in a similar format #UB#.  Technically, the matrix
 * is #M = scale * U'*U# so that #M# may be negative definite as well.
 */
class MatrixSymPosDefBandedChol	: public MatrixSymWithOpFactorized
{
public:
	
	///
	typedef ReferenceCountingPack::ref_count_ptr<
		ResourceManagementPack::ReleaseResource>  release_resource_ptr_t;

	// //////////////
    // Constructors

	///
	/** Construct and Initialize.
	 *
	 * This constructor just calls #this->initialize(...)#.
	 */
	MatrixSymPosDefBandedChol(
		size_type                         n                       = 0
		,size_type                        kd                      = 0
		,GenMatrixSlice                   *MB                     = NULL
		,const release_resource_ptr_t&    MB_release_resource_ptr = NULL
		,BLAS_Cpp::Uplo                   MB_uplo                 = BLAS_Cpp::lower
		,GenMatrixSlice                   *UB                     = NULL
		,const release_resource_ptr_t&    UB_release_resource_ptr = NULL
		,BLAS_Cpp::Uplo                   UB_uplo                 = BLAS_Cpp::lower
		,bool                             update_factor           = false
		,value_type                       scale                   = 1.0
		);

	// ///////////////////////////
    // Access representation

	///
	/** Initialize
	 *
	 * If called with all of the default arguments then #this# will become uninitialized.
	 *
	 * ToDo: Finish pre and post conditions!
	 *
	 * @param  n        [in] Determines the size of the banded matrix (n x n).
	 *                  If #n == 0# then all of the following arguments should be left at
	 *                  their defaults and #this# will become uninitialized.
	 * @param  kd       [in] Determines the band width of the matrix as defined by xPBTRF(...).
	 * @param  MB       [in/state] If #MB != NULL# then this matrix (size (kd+1) x n) is used to store
	 *                  the original banded matrix #MB# in the format of xPBTRF(...).  This matrix must
	 *                  be initialized on input.
	 * @param  MB_release_resource_ptr
	 *                  [in] Only significant if #MB != NULL#.  Points to a resource to
	 *                  be released when #MB# is no longer needed.
	 * @param  MB_uplo  [in] Determines if #MB# is stores the upper or lower triangular elements.
	 * @param  UB       [in/state] If #UB != NULL# then this matrix (size (kd+1) x n) is used to store
	 *                  the Cholesky factor of the banded matrix in the format of xPBTRF(...).
	 *                  This matrix may or may not be initialized on input.
	 *                  If #update_factor == false# this this matrix must  already be initialized.
	 *                  If #update_factor == true# then this matrix will be computed.
2	 *                  If #UB == NULL# then storage for the Cholesky factor will be computed
	 *                  on the fly and will be factored.
	 * @param  UB_release_resource_ptr
	 *                  [in] Only significant if #UB != NULL#.  Points to a resource to
	 *                  be released when #UB# is no longer needed.
	 * @param  UB_uplo  [in] Determines if #UB# is stores the upper or lower triangular elements.
	 * @param  update_factor
	 *                  [in] If true then the factor will be computed within this function call.
	 * @param  scale    [in] Only significant if #MB != NULL# or #UB != NULL# (see intro).
	 */
	void initialize(
		size_type                         n                       = 0
		,size_type                        kd                      = 0
		,GenMatrixSlice                   *MB                     = NULL
		,const release_resource_ptr_t&    MB_release_resource_ptr = NULL
		,BLAS_Cpp::Uplo                   MB_uplo                 = BLAS_Cpp::lower
		,GenMatrixSlice                   *UB                     = NULL
		,const release_resource_ptr_t&    UB_release_resource_ptr = NULL
		,BLAS_Cpp::Uplo                   UB_uplo                 = BLAS_Cpp::lower
		,bool                             update_factor           = false
		,value_type                       scale                   = 1.0
		);

	///
	size_type kd() const;
	///
	/** Get view of MB.
	 */
	GenMatrixSlice& MB();
	///
	const GenMatrixSlice& MB() const;
	///
	BLAS_Cpp::Uplo MB_uplo() const;
	///
	/** Get view of UB.
	 */
	GenMatrixSlice& UB();
	///
	const GenMatrixSlice& UB() const;
	///
	BLAS_Cpp::Uplo UB_uplo() const;

	// /////////////////////////////
	// Overridden from MatrixWithOp

	///
	size_type rows() const;
	///
	size_type nz() const;
	///
	std::ostream& output(std::ostream& out) const;
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const VectorSlice& vs_rhs3, value_type beta) const;
	///
	void Vp_StPtMtV(VectorSlice* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const SpVectorSlice& sv_rhs3, value_type beta) const;

	// //////////////////////////////////
	// Overridden from MatrixFactorized

	/// With throw exception if factorization is not allowed.
	void V_InvMtV(VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2) const;

private:
	
	// /////////////////////////////
	// Private data members

	size_type                       n_;
	size_type                       kd_;
	GenMatrixSlice                  MB_;
	release_resource_ptr_t          MB_release_resource_ptr_;
	BLAS_Cpp::Uplo                  MB_uplo_;
	mutable GenMatrixSlice          UB_;
	mutable release_resource_ptr_t  UB_release_resource_ptr_;
	mutable BLAS_Cpp::Uplo          UB_uplo_;
	mutable bool                    factor_updated_;
	value_type                      scale_;

	// /////////////////////////////
	// Private member functions

	void assert_initialized() const;
	void update_factorization() const;

}; // end class MatrixSymPosDefBandedChol

// ///////////////////////////////////////////////////////
// Inline members for MatrixSymPosDefBandedChol

inline
size_type MatrixSymPosDefBandedChol::kd() const
{
	return kd_;
}

inline
GenMatrixSlice& MatrixSymPosDefBandedChol::MB()
{
	return MB_;
}

inline
const GenMatrixSlice& MatrixSymPosDefBandedChol::MB() const
{
	return MB_;
}

inline
BLAS_Cpp::Uplo MatrixSymPosDefBandedChol::MB_uplo() const
{
	return MB_uplo_;
}

inline
GenMatrixSlice& MatrixSymPosDefBandedChol::UB()
{
	return UB_;
}

inline
const GenMatrixSlice& MatrixSymPosDefBandedChol::UB() const
{
	return UB_;
}

inline
BLAS_Cpp::Uplo MatrixSymPosDefBandedChol::UB_uplo() const
{
	return UB_uplo_;
}

} // end namespace ConstrainedOptimizationPack

#endif // MATRIX_SYM_POS_DEF_BANDED_CHOL_H
