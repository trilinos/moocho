// //////////////////////////////////////////////////////////////////////////////////
// MatrixWithOp.h
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

#ifndef ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_SUB_VIEW_H
#define ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_SUB_VIEW_H

#include <iosfwd>

#include "MatrixWithOp.h"
#include "ref_count_ptr.h"
#include "VectorSpace.h"

namespace AbstractLinAlgPack {

///
/** Standard subclass for representing a sub, possibly transposed, view of a matrix
 * 
 * The matrix \c M_view represented by \c this is:
 \verbatim
 
 M_view = op(M_full(rng_rows,rng_cols))
 \endverbatim
 *
 * ToDo: Finish Documentation!
 */
class MatrixWithOpSubView : public virtual MatrixWithOp {
public:
	
	///
	typedef ReferenceCountingPack::ref_count_ptr<MatrixWithOp>   mat_ptr_t;

	/** @name Constructors/initalizers */
	//@{

	///
	/** Calls <tt>this->initialize(...)</tt>
	 */
	MatrixWithOpSubView(
		const mat_ptr_t&   M_full    = NULL
		,const Range1D&    rng_rows  = Range1D()
		,const Range1D&    rng_cols  = Range1D()
		,BLAS_Cpp::Transp   M_trans  = BLAS_Cpp::no_trans
		);
		
	///
	/** Initialize the view of a matrix.
	 *
	 * @param  M_full   [in] Smart pointer for the matrix to provide a view of.
	 *                  It is allowed for <tt>M_full.get() == NULL</tt> in which case
	 *                  \c this will become uninitialized and none of the rest of the
	 *                  arguments matter and any value will do (i.e. the default values).
	 * @param  rng_rows [in] Range in the rows of <tt>*M_full</tt> that \c this will represent.
	 *                  Only significant if <tt>M_full.get() != NULL</tt>.
	 * @param  rng_cols [in] Range in the columns of <tt>*M_full</tt> that \c this will represent.
	 *                  Only significant if <tt>M_full.get() != NULL</tt>.
	 * @param  M_trans  [in] If <tt>M_trans == no_trans</tt> then \c this will represent
	 *                  <tt>M_full(rng_rows,rng_cols)</tt>.  If If <tt>M_trans == trans</tt>
	 *                  then \c this will represent the transpose <tt>M_full(rng_rows,rng_cols)'</tt>.
	 *
	 * Preconditions:<ul>
	 * <li>[<tt>M_full.get()!=NULL && !rng_rows.full_range()</tt>]
	 *     <tt>rng_rows.ubound() <= M_full->rows()</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li>[<tt>M_full.get()!=NULL && !rng_cols.full_range()</tt>]
	 *     <tt>rng_cols.ubound() <= M_full->cols()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->M_full_ptr().get() == M_full.get()</tt>
	 * <li>[<tt>M_full.get()!=NULL</tt>] <tt>&this->M_full() == M_full.get()</tt>
	 * <li>[<tt>M_full.get()!=NULL && (rng_rows.full_range() || (rng_rows.lbound() == 1 && rng_rows.ubound() == M_full->rows()))</tt>]
	 *     <tt>this->rng_rows().full_range() == true</tt>
	 * <li>[<tt>M_full.get()!=NULL && (rng_cols.full_range() || (rng_cols.lbound() == 1 && rng_cols.ubound() == M_full->cols()))</tt>]
	 *     <tt>this->rng_cols().full_range() == true</tt>
	 * <li>[<tt>M_full.get()!=NULL</tt>] <tt>this->M_trans == M_trans</tt>
	 * <li>[<tt>M_full.get()==NULL</tt>] <tt>this->rng_rows() == Range1D::Invalid</tt>
	 * <li>[<tt>M_full.get()==NULL</tt>] <tt>this->rng_cols() == Range1D::Invalid</tt>
	 * </ul>
	 */
	void initialize(
		const mat_ptr_t&   M_full
		,const Range1D&    rng_rows  = Range1D()
		,const Range1D&    rng_cols  = Range1D()
		,BLAS_Cpp::Transp  M_trans = BLAS_Cpp::no_trans
		);

	//@}

	/** @name Representation access */
	//@{

	///
	const mat_ptr_t& M_full_ptr();
	///
	MatrixWithOp& M_full();
	///
	const MatrixWithOp& M_full() const;
	///
	Range1D rng_rows() const;
	///
	Range1D rng_cols() const;
	///
	BLAS_Cpp::Transp M_trans();

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

	/** @name Overridden from MatrixWithOp */
	//@{

	///
	const VectorSpace& space_cols() const;
	///
	const VectorSpace& space_rows() const;
	///
	MatrixWithOp::mat_ptr_t sub_view(const Range1D& row_rng, const Range1D& col_rng) const;
	///
	void zero_out();
	///
	void Mt_S( value_type alpha );
	///
	MatrixWithOp& operator=(const MatrixWithOp& M);
	///
	std::ostream& output(std::ostream& out) const;
	///
	bool Mp_StM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs) const;
	///
	bool Mp_StMtP(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp M_trans
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		) const;
	///
	bool Mp_StPtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		, BLAS_Cpp::Transp M_trans
		) const;
	///
	bool Mp_StPtMtP(
		MatrixWithOp* mwo_lhs, value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
		) const;
	///
	bool Mp_StM(
		value_type alpha,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs);
	///
	bool Mp_StMtP(
		value_type alpha
		,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		);
	///
	bool Mp_StPtM(
		value_type alpha
		,const GenPermMatrixSlice& P_rhs, BLAS_Cpp::Transp P_rhs_trans
		,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
		);
	///
	bool Mp_StPtMtP(
		value_type alpha
		,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		,const MatrixWithOp& M_rhs, BLAS_Cpp::Transp M_trans
		,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
		);
	///
	void Vp_StMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorWithOp& v_rhs2, value_type beta) const;
	///
	void Vp_StMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;
	///
	void Vp_StPtMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const VectorWithOp& v_rhs3, value_type beta) const;
	///
	void Vp_StPtMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha
		, const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
		, BLAS_Cpp::Transp M_rhs2_trans
		, const SpVectorSlice& sv_rhs3, value_type beta) const;
	///
	value_type transVtMtV(
		const VectorWithOp& v_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const VectorWithOp& v_rhs3) const;
	///
	value_type transVtMtV(
		const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const SpVectorSlice& sv_rhs3) const;
	///
	void syr2k(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		, const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
		, const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
		, value_type beta, MatrixSymWithOp* symwo_lhs ) const;
	///
	bool Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& mwo_rhs2
		, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const;
	///
	bool Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		, const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		, BLAS_Cpp::Transp trans_rhs2, value_type beta ) const;
	///
	bool Mp_StMtM(
		value_type alpha
		,const MatrixWithOp& mvw_rhs1, BLAS_Cpp::Transp trans_rhs1
		,const MatrixWithOp& mwo_rhs2,BLAS_Cpp::Transp trans_rhs2
		,value_type beta );
	///
	void syrk(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		, value_type beta, MatrixSymWithOp* sym_lhs ) const;
	
	//@}

private:
	
#ifdef DOXYGEN_COMPILE	
	MatrixWithOp              *M_full;
	Range1D                   rng_rows;
	Range1D                   rng_cols;
#else
	mat_ptr_t                 M_full_;
	Range1D                   rng_rows_;
	Range1D                   rng_cols_;
	BLAS_Cpp::Transp          M_trans_;
	VectorSpace::space_ptr_t  space_cols_;
	VectorSpace::space_ptr_t  space_rows_;
#endif

	//
	void assert_initialized() const;
	
};	// end class MatrixWithOpSubView

// //////////////////////////////////
// Inline members

inline
const MatrixWithOpSubView::mat_ptr_t&
MatrixWithOpSubView::M_full_ptr()
{
	return M_full_;
}

inline
MatrixWithOp& MatrixWithOpSubView::M_full()
{
	return *M_full_;
}

inline
const MatrixWithOp& MatrixWithOpSubView::M_full() const
{
	return *M_full_;
}

inline
Range1D MatrixWithOpSubView::rng_rows() const
{
	return rng_rows_;
}

inline
Range1D MatrixWithOpSubView::rng_cols() const
{
	return rng_rows_;
}

inline
BLAS_Cpp::Transp MatrixWithOpSubView::M_trans()
{
	return M_trans_;
}

}	// end namespace AbstractLinAlgPack

#endif	// ABSTRACT_LINALG_PACK_MATRIX_WITH_OP_SUB_VIEW_H
