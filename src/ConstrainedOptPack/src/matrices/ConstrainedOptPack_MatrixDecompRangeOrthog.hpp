// //////////////////////////////////////////////////////////////
// MatrixDecompRangeOrthog.hpp
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

#ifndef MATRIX_DECOMP_RANGE_ORTHOG_H
#define MATRIX_DECOMP_RANGE_ORTHOG_H

#include "ConstrainedOptimizationPackTypes.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/MatrixSymOpNonsing.hpp"

namespace ConstrainedOptimizationPack {

///
/** Matrix subclass for variable reduction orthogonal matrix </tt>R = Gc(:,con_decomp)'*Y</tt>.
 *
 * This matrix class is used to represent the matrix:
 \verbatim

 R = C*(I + D*D')

 inv(R) =  inv(I + D*D') * inv(C)

        = (I - D * inv(I + D'*D) * D') * inv(C)
                       \______/
                          S
 \endverbatim
 * Above, the expresion for <tt>inv(R)</tt> is derived using the Sherman-Morrison-Woodbury formula.
 * The nonsingular matrix <tt>S = I + D'*D</tt> is setup by the client, along with the basis matrix
 * \c C and the direct sensitivity matrix \c D.
 */
class MatrixDecompRangeOrthog
	: public AbstractLinAlgPack::MatrixOpNonsing // full path needed for doxygen
{
public:

	/** @name Public types */
	//@{
#ifndef DOXYGEN_COMPILE
	///
	typedef MemMngPack::ref_count_ptr<const MatrixOpNonsing>     C_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<const MatrixOp>                D_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<const MatrixSymOpNonsing>  S_ptr_t;
#endif
	//@}

	/** @name Constructors/initializers */
	//@{

	///
	/** Constructs uninitialized.
	 *
	 * Postconditions:<ul>
	 * <li> Same as for <tt>this->set_uninitialized()</tt>.
	 * </ul>
	 */
	MatrixDecompRangeOrthog();

	/// Calls <tt>this->initialize()</tt>.
	MatrixDecompRangeOrthog(
		const C_ptr_t   &C_ptr
		,const D_ptr_t  &D_ptr
		,const S_ptr_t  &S_ptr
		);

	///
	/** Initialize the matrix object.
	 *
	 * @param  C_ptr  [in]
	 * @param  D_ptr  [in]
	 * @param  S_ptr  [in]
	 *
	 * Preconditions:<ul>
	 * <li> <tt>C_ptr.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>D_ptr.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>S_ptr.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>C_ptr->space_rows().is_compatible(D_ptr->space_cols()) == true</tt>
	 *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>S_ptr->space_cols().is_compatible(D_ptr->space_rows()) == true</tt>
	 *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->C_ptr().get() == C_ptr.get()</tt>
	 * <li> <tt>this->D_ptr().get() == D_ptr.get()</tt>
	 * <li> <tt>this->S_ptr().get() == S_ptr.get()</tt>
	 * <li> <tt>this->space_cols().is_compatible(C_ptr->space_cols())</tt>
	 * <li> <tt>this->space_rows().is_compatible(C_ptr->space_rows())</tt>
	 * </ul>
	 */
	void initialize(
		const C_ptr_t   &C_ptr
		,const D_ptr_t  &D_ptr
		,const S_ptr_t  &S_ptr
		);

	///
	/** Make uninitialized.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->C_ptr().get() == NULL</tt>
	 * <li> <tt>this->D_ptr().get() == NULL</tt>
	 * <li> <tt>this->S_ptr().get() == NULL</tt>
	 * <li> <tt>this->rows() == 0</tt>
	 * <li> <tt>this->cols() == 0</tt>
	 * </ul>
	 */
	void set_uninitialized();

	//@}

	/** @name Access */
	//@{

	///
	const C_ptr_t& C_ptr() const;
	///
	const D_ptr_t& D_ptr() const;
	///
	const S_ptr_t& S_ptr() const;

	//@}
	
	/** @name Overridden from MatrixOp */
	//@{

	///
	size_type rows() const;
	///
	size_type cols() const;
	///
	const VectorSpace& space_cols() const;
	///
	const VectorSpace& space_rows() const;
	///
	std::ostream& output(std::ostream& out) const;
	///
	void Vp_StMtV(
		VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const Vector& v_rhs2, value_type beta ) const;

	//@}

	/** @name Overridden from MatrixOpNonsing */
	//@{

	///
	void V_InvMtV(
		VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const Vector& v_rhs2 ) const;

	//@}

private:

#ifdef DOXYGEN_COMPILE
	AbstractLinAlgPack::MatrixOpNonsing      *C;
	AbstractLinAlgPack::MatrixOp                 *D;
	AbstractLinAlgPack::MatrixSymOpNonsing   *S;
#else
	C_ptr_t    C_ptr_;
	D_ptr_t    D_ptr_;
	S_ptr_t    S_ptr_;
#endif

	//
	void assert_initialized(const char func_name[]) const;

	// not defined and not to be called
	MatrixDecompRangeOrthog(const MatrixDecompRangeOrthog&);
	MatrixDecompRangeOrthog& operator=(const MatrixDecompRangeOrthog&);

}; // end class MatrixDecompRangeOrthog

// /////////////////////////////////
// Inline members

inline
const MatrixDecompRangeOrthog::C_ptr_t&
MatrixDecompRangeOrthog::C_ptr() const
{
	return C_ptr_;
}

inline
const MatrixDecompRangeOrthog::D_ptr_t&
MatrixDecompRangeOrthog::D_ptr() const
{
	return D_ptr_;
}

inline
const MatrixDecompRangeOrthog::S_ptr_t&
MatrixDecompRangeOrthog::S_ptr() const
{
	return S_ptr_;
}

} // end namespace ConstrainedOptimizationPack

#endif // MATRIX_DECOMP_RANGE_ORTHOG_H
