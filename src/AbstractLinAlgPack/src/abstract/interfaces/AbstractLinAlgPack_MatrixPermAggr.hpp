// /////////////////////////////////////////////////////////////////////
// MatrixPermAggr.hpp
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

#ifndef MATRIX_PERM_AGGR_H
#define MATRIX_PERM_AGGR_H

#include "AbstractLinAlgPack/src/abstract/interfaces/MatrixOp.hpp"

namespace AbstractLinAlgPack {

///
/** Aggregate matrix class for a matrix and its permuted view.
 *
 * <tt>mat_perm = row_perm * mat_orig * col_perm'</tt>.
 *
 */
class MatrixPermAggr
	: virtual public MatrixOp
{
public:

	/** @name Public types */
	//@{

	///
	typedef MemMngPack::ref_count_ptr<const Permutation>   perm_ptr_t;

	//@}

	/** @name Constructors / initializers */
	//@{

	/// Construct to uninitialized
	MatrixPermAggr();

	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	MatrixPermAggr(
		const mat_ptr_t      &mat_orig
		,const perm_ptr_t    &row_perm
		,const perm_ptr_t    &col_perm
		,const mat_ptr_t     &mat_perm
		);

	///
	/** Initialize.
	 *
	 * <tt>mat_perm = row_perm' * mat_orig * col_perm</tt>.
	 *
	 * @param  mat_orig  [in] Smart pointer to original unpermuted matrix.
	 * @param  row_perm  [in] Smart pointer to row permutation.  If <tt>row_perm.get() == NULL</tt>
	 *                   then the identity permutation is assumed.
	 * @param  col_perm  [in] Smart pointer to column permutation.  If <tt>col_perm.get() == NULL</tt>
	 *                   then the identity permutation is assumed.
	 * @param  mat_perm  [in] Smart pointer to permuted matrix.  It is allowed for
	 *                   <tt>mat_perm.get() == NULL</tt> in which case all of the linear algebra
	 *                   methods are implemented in terms of \c mat_orig, \c row_perm and \c col_perm.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>mat_perm.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> [<tt>row_perm.get() != NULL</tt>] <tt>mat_orig->space_cols().is_compatible(row_perm->space()) == true</tt>
	 *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * <li> [<tt>col_perm.get() != NULL</tt>] <tt>mat_orig->space_rows().is_compatible(col_perm->space()) == true</tt>
	 *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->mat_orig().get() == mat_orig.get()</tt>
	 * <li> <tt>this->row_perm().get() == row_perm.get()</tt>
	 * <li> <tt>this->col_perm().get() == col_perm.get()</tt>
	 * <li> <tt>this->mat_perm().get() == mat_perm.get()</tt>
	 * </ul>
	 */
	void initialize(
		const mat_ptr_t      &mat_orig
		,const perm_ptr_t    &row_perm
		,const perm_ptr_t    &col_perm
		,const mat_ptr_t     &mat_perm
		);

	///
	/** Set uninitialized.
	 *
	 * ToDo: Finish documentation.
	 */
	void set_uninitialized();

	//@}
	
	/** @name Access */
	//@{

	///
	const mat_ptr_t& mat_orig() const;
	///
	const perm_ptr_t& row_perm() const;
	///
	const perm_ptr_t& col_perm() const;
	///
	const mat_ptr_t& mat_perm() const;

	//@}

	/** @name Overridden from MatrixBase */
	//@{

	///
	size_type rows() const;
	///
	size_type cols() const;
	///
	size_type nz() const;

	//@}

	/** @name Overridden from MatrixOp */
	//@{

	///
	const VectorSpace& space_cols() const;
	///
	const VectorSpace& space_rows() const;
	///
	MatrixOp::mat_ptr_t sub_view(const Range1D& row_rng, const Range1D& col_rng) const;
	///
	MatrixOp& operator=(const MatrixOp& M);
	///
	std::ostream& output(std::ostream& out) const;

protected:

	///
	bool Mp_StM(
		MatrixOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs) const;
	///
	bool Mp_StMtP(
		MatrixOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp M_trans
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		) const;
	///
	bool Mp_StPtM(
		MatrixOp* mwo_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		, BLAS_Cpp::Transp M_trans
		) const;
	///
	bool Mp_StPtMtP(
		MatrixOp* mwo_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
		) const;
	///
	void Vp_StMtV(
		VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const Vector& v_rhs2, value_type beta) const;
	///
	void Vp_StMtV(
		VectorMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	void Vp_StPtMtV(
		VectorMutable* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const Vector& v_rhs3, value_type beta) const;
	///
	void Vp_StPtMtV(
		VectorMutable* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const SpVectorSlice& sv_rhs3, value_type beta) const;
	///
	value_type transVtMtV(
		const Vector& v_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const Vector& v_rhs3) const;
	///
	value_type transVtMtV(
		const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const SpVectorSlice& sv_rhs3) const;
	///
	void syr2k(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
		, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
		, value_type beta, MatrixSymOp* symwo_lhs ) const;
	///
	bool Mp_StMtM(
		MatrixOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs1, const MatrixOp& mwo_rhs2
		, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const;
	///
	bool Mp_StMtM(
		MatrixOp* mwo_lhs, value_type alpha
		, const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const;
	///
	bool syrk(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		, value_type beta, MatrixSymOp* sym_lhs ) const;
	
	//@}

private:

#ifdef DOXYGEN_COMPILE
	MatrixOp         *mat_orig;
	Permutation          *row_perm;
	Permutation          *col_perm;
	MatrixOp         *mat_perm;
#else
	mat_ptr_t            mat_orig_;
	perm_ptr_t           row_perm_;
	perm_ptr_t           col_perm_;
	mat_ptr_t            mat_perm_;
#endif

}; // end class MatrixPermAggr

// ////////////////////////////////////
// Inline members

inline
const MatrixOp::mat_ptr_t&
MatrixPermAggr::mat_orig() const
{
	return mat_orig_;
}

inline
const MatrixPermAggr::perm_ptr_t&
MatrixPermAggr::row_perm() const
{
	return row_perm_;
}

inline
const MatrixPermAggr::perm_ptr_t&
MatrixPermAggr::col_perm() const
{
	return col_perm_;
}

inline
const MatrixOp::mat_ptr_t& MatrixPermAggr::mat_perm() const
{
	return mat_perm_;
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_PERM_AGGR_H