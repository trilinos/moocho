// //////////////////////////////////////////////////////////////////////////////////
// MatrixVarReductImplicit.hpp
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

#ifndef MATRIX_VAR_REDUCT_IMPLICIT_H
#define MATRIX_VAR_REDUCT_IMPLICIT_H

#include <vector>
#include <list>

#include "ConstrainedOptimizationPackTypes.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/MatrixOp.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorSpace.hpp"
#include "ref_count_ptr.hpp"

namespace ConstrainedOptimizationPack {

///
/** Implements <tt>D = - inv(C) * N</tt> for a variable reduction projection.
 *
 * This class is used to implement an implicit matrix defined as
 \verbatim
 
 D = - inv(C) * N
 \endverbatim
 * The operations <tt>y = op(D)*x</tt> are implemented as:
 \verbatim

	y = D * x
	  = -inv(C) * (N * x)

	y = D' * x
	  = - N' * (inv(C') * x)
 \endverbatim
 * This class also allows the client to set a precomputed matrix \c D_direct
 * that represents <tt>D = -inv(C)*N</tt> that will be used to extract rows or columns
 * or to implement operations involving <tt>GenPermMatrixSlice</tt>
 * when convenient.  One might ask why this subclass would even be used if
 * \c D_direct was even available.  The reason is that it may be cheaper to perform
 * the sparse solve and matrix-vector multiplication with \c C and \c N than it is
 * to use a dense precomputed matrix \c D_direct.  Determining if this class is
 * even useful when \c D_direct is availible must be determined at runtime
 * using timing data (which can be very hard to do well in general).
 *
 * This implementation is designed to deal efficiently with the case where matrix-
 * vector multiplications will only be performed with subsets of rows of 
 * <tt>inv(C)*N</tt> or columns of <tt>N'*inv(C')</tt>.  This primarily affects
 * two types of operations:<ul>
 * <li> <tt>y = b*y + a*[-N'*inv(C')]*x</tt> (\c x is a <tt>SpVectorSlice</tt> object)
 * <li> <tt>y = b*y + a*op(P)*[-inv(C)*N]*x</tt> (\c P has few nonzeros, x is any vector)
 * </ul>
 * when \c D_direct is not set then needed rows of <tt>inv(C)*N</tt> are generated on the fly
 * (as abstract vectors) and stored away for later use.  When <tt>this->initialize()</tt> is
 * called then all of these computed rows are discarded and they must be generated again.
 */
class MatrixVarReductImplicit
	: public AbstractLinAlgPack::MatrixOp
{
public:

	/** @name Public types */
	//@{
	
	///
	typedef MemMngPack::ref_count_ptr<const MatrixOpNonsing>   mat_nonsing_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<const MatrixOp>              mat_ptr_t;

	//@}

	/** @name Constructors / initializers */
	//@{

	///
	/** Initialize \c this matrix object.
	 *
	 * @param  C  [in] Nonsingular basis matrix object.  This matrix object must
	 *            not be altered while \c this is in use or until
	 *            <tt>this->initialize()</tt> is called again.
	 * @param  N  [in] Genaral nonbasis matrix object.  This matrix object must
	 *            not be altered while \c this is in use or until
	 *            <tt>this->initialize()</tt> is called again.
	 * @param  D_direct
	 *            [in] Matrix object for <tt>D = -inv(C)*N</tt> already
	 *            computed.  The matrix object \c D_direct will not be modifed
	 *            by \c this and must not be altered while c\ this matrix object
	 *            is in use or until <tt>this->initialize()</tt> is called again.
	 *            <tt>D_direct == NULL</tt> is allowed and \c this matrix object
	 *            will just have to do without.  For most applications (except
	 *            those using direct linear solvers for \c C and when \c N has
	 *            many columns) <tt>D_direct</tt> should be set to <tt>NULL</tt>
	 *            and should not be computed by the client (that is the whole
	 *            purpose for this matrix class).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>C.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>N.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> [<tt>D_direct != NULL</tt>]
	 *      <tt>D_direct->space_cols().is_compatible(C->space_cols()) == true
	 *      && D_direct->space_rows().is_compatible(N->space_rows()) == true</tt>
	 *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->C_ptr().get() == C.get()</tt>
	 * <li> <tt>this->N_ptr().get() == N.get()</tt>
	 * <li> <tt>this->D_direct_ptr().get() == N.get()</tt>
	 * </ul>
	 */
	virtual void initialize(
		const mat_nonsing_ptr_t          &C
		,const mat_ptr_t                 &N
		,const mat_ptr_t                 &D_direct
		);
	///
	/** Set the matrix to uninitialized.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->C_ptr().get() == NULL</tt>
	 * <li> <tt>this->N_ptr().get() == NULL</tt>
	 * <li> <tt>this->D_direct_ptr().get() == NULL</tt>
	 * </ul>
	 */
	virtual void set_uninitialized();

	//@}

	/** @name Access */
	//@{

	///
	/** Return the smart pointer to the aggregate basis matrix object \c C.
	 *
	 * If <tt>this</tt> is the only reference to this matrix object, then
	 * <tt>return.count() == 1</tt> will be true.
	 */
	const mat_nonsing_ptr_t& C_ptr() const;
	
	///
	/** Return the smart pointer to the aggregate nonbasis matrix object \c N.
	 *
	 * If <tt>this</tt> is the only reference to this matrix object, then
	 * <tt>return.count() == 1</tt> will be true.
	 */
	const mat_ptr_t& N_ptr() const;
	
	///
	/** Return the smart pointer to the aggregate precomputed matrix object \c D_direct (if set).
	 *
	 * If <tt>this</tt> is the only reference to this matrix object, then
	 * <tt>return.count() == 1</tt> will be true.
	 */
	const mat_ptr_t& D_direct_ptr() const;
	
	//@}

	/** @name Overridden from MatrixBase. */
	//@{
	///
	size_type rows() const;
	///
	size_type cols() const;
	//@}

	/** @name Overridden from MatrixOp. */
	//@{

	///
	const VectorSpace& space_cols() const;
	///
	const VectorSpace& space_rows() const;
	///
	MatrixOp& operator=(const MatrixOp& M);
	///
	std::ostream& output(std::ostream&) const;
	///
	void Vp_StMtV(
		VectorMutable* v_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const Vector& v_rhs2, value_type beta
		) const;
	///
	void Vp_StMtV(
		VectorMutable* v_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const SpVectorSlice& sv_rhs2, value_type beta
		) const;
	///
	void Vp_StPtMtV(
		VectorMutable* v_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,BLAS_Cpp::Transp M_rhs2_trans
		,const Vector& v_rhs3, value_type beta
		) const;
	///
	void Vp_StPtMtV(
		VectorMutable* v_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,BLAS_Cpp::Transp M_rhs2_trans
		,const SpVectorSlice& sv_rhs3, value_type beta
		) const;
	
	//@}
	
private:
	
	// //////////////////////////
	// Private types
	
	typedef std::vector<VectorSpace::vec_mut_ptr_t>    InvCtN_rows_t;
	typedef std::list<index_type>                      InvCtN_rows_set_list_t;

	// //////////////////////////
	// Private data members

#ifdef DOXYGEN_COMPILE
	AbstractLinAlgPack::MatrixOpNonsing  *C;
	AbstractLinAlgPack::MatrixOp             *N;
	AbstractLinAlgPack::MatrixOp             *D_direct;
#else
	mat_nonsing_ptr_t                   C_;
	mat_ptr_t                           N_;
	mat_ptr_t                           D_direct_;

	mutable InvCtN_rows_t               InvCtN_rows_;
	// InvCtN_rows_ keeps track of a set pointers of computed rows of inv(C)*N.
	// If D_direct_ is setup then InvCtN_rows_ is not necessary.  However, if
	// not, then InvCtN_rows_[j-1] will be !=NULL if it points to the precomputed
	// jth row of inv(C)*N and will be ==NULL if this row has not been computed
	// yet.  Each time initialize(...) is called, these rows are deallocated and
	// InvCtN_rows_[j-1], j=1...this->rows() is set to NULL.

	mutable InvCtN_rows_set_list_t      InvCtN_rows_set_list_;
	// InvCtN_rows_set_list_ keeps a unorderd=ed list of the row indexes that are
	// currently updated.  Keeping this list allows the vectors allocated in 
	// InvCtN_rows_[] to be deallocated much faster than if the whole array
	// had to be searched everytime that this->initialize() was called again.

#endif

	// //////////////////////////////////
	// Private member functions

	///
	void assert_initialized() const;

};	// end class MatrixVarReductImplicit

// ////////////////////////////////
// Inline members

inline
const MatrixVarReductImplicit::mat_nonsing_ptr_t&
MatrixVarReductImplicit::C_ptr() const
{
	return C_;
}

inline
const MatrixVarReductImplicit::mat_ptr_t&
MatrixVarReductImplicit::N_ptr() const
{
	return N_;
}

inline
const MatrixVarReductImplicit::mat_ptr_t&
MatrixVarReductImplicit::D_direct_ptr() const
{
	return D_direct_;
}

}	// end namespace ConstrainedOptimizationPack 

#endif	// MATRIX_VAR_REDUCT_IMPLICIT_H
