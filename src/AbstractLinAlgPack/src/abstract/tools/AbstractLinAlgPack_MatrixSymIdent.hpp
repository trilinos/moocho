// ///////////////////////////////////////////////////////////////////////
// MatrixSymIdentity.h
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

#ifndef ALAP_MATRIX_SYM_IDENTITY_H
#define ALAP_MATRIX_SYM_IDENTITY_H

#include "MatrixSymWithOpNonsingular.h"
#include "VectorSpace.h"

namespace AbstractLinAlgPack {

///
/** Matrix subclass for a scaled identity matrix.
 *
 * More operations will be overridden as they are needed by various applications.
 */
class MatrixSymIdentity : virtual public MatrixSymWithOpNonsingular {
public:
	
	/** @name Constructors/initializers */
	//@{

	/// Calls <tt>this->initialize()</tt>.
	MatrixSymIdentity(
		const VectorSpace::space_ptr_t&          vec_space = NULL
		,const value_type                        scale     = 1.0
		);

	///
	void initialize(
		const VectorSpace::space_ptr_t&          vec_space
		,const value_type                        scale
	);

	//@}

	/** @name Access */
	//@{

	///
	value_type scale() const;

	//@}

	/** @name Overridden from MatrixBase */
	//@{

	/// Returns 0 if not initalized.
	size_type rows() const;
	/// Returns <tt>this->rows()</tt>
	size_type nz() const;

	//@}

	/** @name Overridden from MatrixWithOp */
	//@{

	///
	const VectorSpace& space_cols() const;
	///
	std::ostream& output(std::ostream& out) const;
	///
	void Vp_StMtV(
		VectorWithOpMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		,const VectorWithOp& v_rhs2, value_type beta ) const;

	//@}

	/** @name Overridden from MatrixNonsingular */
	//@{

	///
	void V_InvMtV(
		VectorWithOpMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		,const VectorWithOp& v_rhs2 ) const;

	//@}

private:

	VectorSpace::space_ptr_t  vec_space_;
	value_type                scale_;
	
}; // end class MatrixSymIdentity

// ///////////////////////////////////////////
// Inline members

inline
value_type MatrixSymIdentity::scale() const
{
	return scale_;
}

} // end namespace AbstractLinAlgPack

#endif // ALAP_MATRIX_SYM_IDENTITY_H
