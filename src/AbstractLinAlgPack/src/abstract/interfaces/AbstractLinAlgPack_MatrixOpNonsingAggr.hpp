// /////////////////////////////////////////////////////////////////////
// MatrixOpNonsingAggr.hpp
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

#ifndef MATRIX_WITH_OP_NONSINGULAR_AGGR_H
#define MATRIX_WITH_OP_NONSINGULAR_AGGR_H

#include "MatrixOpNonsing.hpp"

namespace AbstractLinAlgPack {

///
/** Aggregate matrix class pulling together a \c MatrixOp object and a
 * \c MatrixNonsing object into a unified matrix object.
 *
 * ToDo: Finish documentation!
 */
class MatrixOpNonsingAggr
	: virtual public MatrixOpNonsing
{
public:

	/** @name Public types */
	//@{

	///
	typedef MemMngPack::ref_count_ptr<const MatrixOp>        mwo_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<const MatrixNonsing>   mns_ptr_t;

	//@}

	/** @name Constructors / initializers */
	//@{

	/// Construct to uninitialized
	MatrixOpNonsingAggr();

	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	MatrixOpNonsingAggr(
		const mwo_ptr_t       &mwo
		,BLAS_Cpp::Transp     mwo_trans
		,const mns_ptr_t      &mns
		,BLAS_Cpp::Transp     mns_trans
		);

	///
	/** Initialize.
	 *
	 * ToDo: Finish documentation.
	 */
	void initialize(
		const mwo_ptr_t       &mwo
		,BLAS_Cpp::Transp     mwo_trans
		,const mns_ptr_t      &mns
		,BLAS_Cpp::Transp     mns_trans
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
	const mwo_ptr_t& mwo() const;
	///
	BLAS_Cpp::Transp mwo_trans() const;
	///
	const mns_ptr_t& mns() const;
	///
	BLAS_Cpp::Transp mns_trans() const;

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

	/** @name Overridden from MatrixNonsing */
	//@{

	///
	void V_InvMtV(
		VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		,const Vector& v_rhs2) const;
	///
	void V_InvMtV(
		VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2) const;
	///
	value_type transVtInvMtV(
		const Vector& v_rhs1
		,BLAS_Cpp::Transp trans_rhs2, const Vector& v_rhs3) const;
	///
	value_type transVtInvMtV(
		const SpVectorSlice& sv_rhs1
		,BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const;
	///
	void M_StInvMtM(
		MatrixOp* m_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
		) const;
	///
	void M_StMtInvM(
		MatrixOp* m_lhs, value_type alpha
		,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		,BLAS_Cpp::Transp trans_rhs2
		) const;

	//@}

private:

#ifdef DOXYGEN_COMPILE
	MatrixOp         *mwo;
	MatrixNonsing    *mns;
#else
	mwo_ptr_t            mwo_;
	BLAS_Cpp::Transp     mwo_trans_;
	mns_ptr_t            mns_;
	BLAS_Cpp::Transp     mns_trans_;
#endif

}; // end class MatrixOpNonsingAggr

// ////////////////////////////////////
// Inline members

inline
const MatrixOpNonsingAggr::mwo_ptr_t&
MatrixOpNonsingAggr::mwo() const
{
	return mwo_;
}

inline
BLAS_Cpp::Transp MatrixOpNonsingAggr::mwo_trans() const
{
	return mwo_trans_;
}

inline
const MatrixOpNonsingAggr::mns_ptr_t&
MatrixOpNonsingAggr::mns() const
{
	return mns_;
}

inline
BLAS_Cpp::Transp MatrixOpNonsingAggr::mns_trans() const
{
	return mns_trans_;
}

} // end namespace AbstractLinAlgPack

#endif // MATRIX_WITH_OP_NONSINGULAR_AGGR_H
