// /////////////////////////////////////////////////////////////////
// MatrixSymIdentitySerial.h
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

#ifndef COP_MATRIX_SYM_IDENTITY_SERIAL_H
#define COP_MATRIX_SYM_IDENTITY_SERIAL_H

#include "MatrixExtractInvCholFactor.h"
#include "SparseLinAlgPack/include/MatrixSymWithOpNonsingularSerial.h"

namespace ConstrainedOptimizationPack {

///
/** Matrix class for a serial scaled identity matrix.
 *
 * More operations will be overridden as they are needed by various applications.
 */
class MatrixSymIdentitySerial
	: public MatrixSymWithOpNonsingularSerial
	, public MatrixExtractInvCholFactor
{
public:
	
    /** @name Constructors/initalizers */
	//@{

	/// Calls <tt>this->initalize()</tt>
	MatrixSymIdentitySerial( size_type size = 1, value_type scale = 1.0 );

	///
	void initalize( size_type size, value_type scale );

	//@}

	/** @name Access */
	//@{

	///
	value_type scale() const;

	//@}

	/** Overridden from MatrixBase */
	//@{

	///
	size_type rows() const;
	///
	size_type nz() const;

	//@}

	/** Overridden from MatrixWithOp */
	//@{

	///
	std::ostream& output(std::ostream& out) const;

	//@}

	/** Overridden from MatrixWithOpSerial */
	//@{

	///
	void Vp_StMtV(
		VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		,const VectorSlice& vs_rhs2, value_type beta) const;

	//@}

	/** @name Overridden from MatrixNonsingularSerial */
	//@{

	///
	void V_InvMtV(
		VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1,const VectorSlice& vs_rhs2 ) const;

	//@}

	/** @name Overridden from MatrixSymNonsingular */
	//@{

	///
	void M_StMtInvMtM(
		sym_gms* sym_gms_lhs, value_type alpha
		, const MatrixWithOp& mwo, BLAS_Cpp::Transp mwo_trans, EMatrixDummyArg
		) const;

	//@}

	/** @name Overridden from MatrixExtractInvCholFactor */
	//@{

	///
	void extract_inv_chol( tri_ele_gms* InvChol ) const;

	//@}

private:

	size_type   size_;
	value_type  scale_;

	
}; // end class MatrixSymIdentitySerial

// //////////////////////////////////////
// Inline members

inline
value_type MatrixSymIdentitySerial::scale() const
{
	return scale_;
}

} // end namespace ConstrainedOptimizationPack

#endif // COP_MATRIX_SYM_IDENTITY_SERIAL_H