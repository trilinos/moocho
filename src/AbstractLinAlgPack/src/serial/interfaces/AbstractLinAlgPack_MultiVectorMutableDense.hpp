// /////////////////////////////////////////////////////////////////////////
// MultiVectorMutableDense.hpp
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

#ifndef MULTI_VECTOR_MUTABLE_DENSE_H
#define MULTI_VECTOR_MUTABLE_DENSE_H

#include "AbstractLinAlgPack/src/serial/interfaces/MatrixOpSerial.hpp"
#include "AbstractLinAlgPack/src/serial/interfaces/MatrixOpGetGMSMutable.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/MultiVectorMutable.hpp"
#include "DenseLinAlgPack/src/DMatrixClass.hpp"
#include "ReleaseResource.hpp"

namespace AbstractLinAlgPack {

///
/** <tt>MultiVectorMutable</tt> "Adapter" subclass for <tt>DenseLinAlgPack::DMatrixSlice</tt>
 * or <tt>DenseLinAlgPack::DMatrix</tt> object.
 *
 * This class can be used either as a view of a <tt>DenseLinAlgPack::DMatrixSlice</tt> object
 * or as a storage type for a <tt>DenseLinAlgPack::DMatrix</tt> object.
 *
 * To create a storage type with the dimensions of <tt>rows x cols</tt> just call the
 * constructor <tt>MatrixWithOpMutableDense(rows,cols)</tt> or after construction you
 * can call <tt>this->initialize(rows,cols)</tt>.
 *
 * To simply create a view of a matrix or its transpose, say \c M, without ownership
 * just call <tt>MatrixWithOpMutableDense(M(),M_trans,NULL)</tt> or after construction call
 * <tt>this->initialize(M(),M_trans,NULL)</tt>.
 *
 * Alternately, \c this can be given a matrix with the responsibility to delete any associated
 * memory by calling <tt>this->initialize()</tt> with a <tt>ReleaseResource</tt> object to
 * perform the deallocation.
 *
 * If \c this has been initialized by <tt>this->initialize(rows,cols)</tt> and
 * if the client really needs to get at the <tt>DenseLinAlgPack::DMatrix</tt> object
 * itself, then it can be obtained as:
 \code
 void f( MatrixWithOpMutableDense* M )
     namespace rmp = MemMngPack;
     DMatrix &_M = *dynamic_cast<rmp::ReleaseResource_ref_count_ptr<DMatrix&> >(*M.gms_release()).ptr;
 \endcode
 * This is not pretty but it is not supposed to be.  Of course the above function will throw
 * an exception if the <tt>dynamic_cast<></tt> fails.
 */
class MultiVectorMutableDense
	: public AbstractLinAlgPack::MultiVectorMutable   // doxygen needs the full path
	, public MatrixOpSerial
	, public MatrixOpGetGMS
{
public:

	///
	typedef MemMngPack::ref_count_ptr<
		MemMngPack::ReleaseResource>  release_resource_ptr_t;

	/** @name Constructors / initializers */
	//@{

	///
	/** Calls <tt>this->initialize(rows,cols)</tt>.
	 */
	MultiVectorMutableDense(
		const size_type                    rows = 0
		,const size_type                   cols = 0
		);
	///
	/** Calls <tt>this->initialize(gms,gms_trans,gms_release)</tt>.
	 */
	MultiVectorMutableDense(
		DMatrixSlice                     gms
		,BLAS_Cpp::Transp                  gms_trans
		,const release_resource_ptr_t&     gms_release
		);
	///
	/** Call <tt>this->initialize(v,v_release)</tt> with an allocated <tt>DenseLinAlgPack::DVector</tt>
	 * object.
	 */
	void initialize(
		const size_type                    rows
		,const size_type                   cols
		);
	///
	/** Initialize with a dense matrix slice.
	 *
	 * Note that solve of the method overrides have to allocate a temporary memory
	 * if <tt>gms_trans != no_trans</tt> (see \c get_gms_view() and \c output()).
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->rows() = ( gms_trans = no_trans ? gms.rows() : gms.cols())</tt>
	 * <li> <tt>this->cols() = ( gms_trans = no_trans ? gms.cols() : gms.rows())</tt>
	 * </ul>
	 */
	void initialize(
		DMatrixSlice                     gms
		,BLAS_Cpp::Transp                  gms_trans
		,const release_resource_ptr_t&     gms_release
		);

	//@}

	/** @name Access */
	//@{
	
	///
	/** Return the non-const dense matrix.
	 *
	 * Note that calling this method may result in the matrix implementation
	 * being modified.  Therefore, no other methods on \c this object should be
	 * called until the <tt>DMatrixSlice</tt> returned from this method is
	 * discarded.
	 */
	DMatrixSlice set_gms();
	///
	/** Return a const dense matrix.
	 */
	const DMatrixSlice get_gms() const;
	///
	/** Return if underlying matrix is being viewed as the transpose or non-transposed.
	 */
	BLAS_Cpp::Transp gms_trans() const;
	///
	/** Return a <tt>ref_count_ptr<></tt> pointer to the object that will
	 * release the associated resource.
	 */
	const release_resource_ptr_t& gms_release() const;

	//@}

	/** @name Overridden from MatrixOpGetGMS */
	//@{

	///
	const DMatrixSlice get_gms_view() const;
	///
	void free_gms_view(const DMatrixSlice* gms_view) const;

	//@}

	/** @name Overridden from MatrixOpGetGMSMutable */
	//@{

	///
	DMatrixSlice get_gms_view();
	///
	void commit_gms_view(DMatrixSlice* gms_view);

	//@}
	
	/** @name Overridden from MatrixBase */
	//@{

	///
	size_type rows() const;
	///
	size_type cols() const;

	//@}

	/** @name Overridden from MatrixOp */
	//@{

	///
	void zero_out();
	///
	void Mt_S( value_type alpha );
	///
	MatrixOp& operator=(const MatrixOp& mwo_rhs);
	///
	std::ostream& output(std::ostream& out) const;
	///
	bool syrk(
		 BLAS_Cpp::Transp M_trans, value_type alpha
		,value_type beta, MatrixSymOp* sym_lhs
		) const;
	///
	bool Mp_StMtM(
		MatrixOp* mwo_lhs, value_type alpha
		,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		,BLAS_Cpp::Transp trans_rhs2
		,value_type beta ) const;
	///
	bool Mp_StMtM(
		MatrixOp* mwo_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
		,value_type beta ) const;

	//@}

	/** @name Overridden from MultiVector */
	//@{

	///
	access_by_t access_by() const;

	//@}

	/** @name Overridden from MultiVectorMutable */
	//@{

	///
	vec_mut_ptr_t col(index_type j);
	///
	vec_mut_ptr_t row(index_type i);
	///
	vec_mut_ptr_t diag(int k);
	///
	multi_vec_mut_ptr_t mv_sub_view(const Range1D& row_rng, const Range1D& col_rng);

	//@}

	/** @name Overridden from MatrixOpSerial */
	//@{

	///
	void Vp_StMtV(
		DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const DVectorSlice& vs_rhs2, value_type beta) const;
	///
	void Vp_StMtV(
		DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;

	// ToDo: Add more overrides as they are needed!

	//@}

protected:

	/** @name Overridden from MultiVector */
	//@{

	///
	void apply_op(
		EApplyBy apply_by, const RTOpPack::RTOp& primary_op
		,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
		,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
		,RTOp_ReductTarget reduct_objs[]
		,const index_type primary_first_ele   , const index_type primary_sub_dim   , const index_type primary_global_offset
		,const index_type secondary_first_ele , const index_type secondary_sub_dim 
		) const;
	///
	void apply_op(
		EApplyBy apply_by, const RTOpPack::RTOp& primary_op, const RTOpPack::RTOp& secondary_op
		,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
		,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type primary_first_ele   , const index_type primary_sub_dim   , const index_type primary_global_offset
		,const index_type secondary_first_ele , const index_type secondary_sub_dim 
		) const;

	//@}

private:

	// ///////////////////////////////////////
	// Private data members
	
	DMatrixSlice            gms_;
	BLAS_Cpp::Transp          gms_trans_;
	release_resource_ptr_t    gms_release_;

}; // end class MultiVectorMutableDense			

// //////////////////////////////////////
// Inline members

inline
DMatrixSlice
MultiVectorMutableDense::set_gms()
{
	return gms_;
}

inline
const DMatrixSlice
MultiVectorMutableDense::get_gms() const
{
	return gms_;
}

inline
BLAS_Cpp::Transp
MultiVectorMutableDense::gms_trans() const
{
	return gms_trans_;
}

inline
const MultiVectorMutableDense::release_resource_ptr_t&
MultiVectorMutableDense::gms_release() const
{
	return gms_release_;
}

} // end namespace AbstractLinAlgPack

#endif // MULTI_VECTOR_MUTABLE_DENSE_H
