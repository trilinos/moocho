// /////////////////////////////////////////////////////////////////////////
// MultiVectorMutableDense.cpp
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

#include <assert.h>

#include "SparseLinAlgPack/include/MultiVectorMutableDense.h"
#include "SparseLinAlgPack/include/VectorWithOpMutableDense.h"
#include "SparseLinAlgPack/include/SpVectorOp.h"
#include "LinAlgPack/include/GenMatrixOp.h"
#include "LinAlgPack/include/GenMatrixOut.h"
#include "ReleaseResource_ref_count_ptr.h"
#include "WorkspacePack.h"
#include "ThrowException.h"

namespace SparseLinAlgPack {

// Constructors / initializers

MultiVectorMutableDense::MultiVectorMutableDense(
	const size_type                    rows
	,const size_type                   cols
	)
{
	this->initialize(rows,cols);
}

MultiVectorMutableDense::MultiVectorMutableDense(
	GenMatrixSlice                     gms
	,BLAS_Cpp::Transp                  gms_trans
	,const release_resource_ptr_t&     gms_release
	)
{
	this->initialize(gms,gms_trans,gms_release);
}

void MultiVectorMutableDense::initialize(
	const size_type                    rows
	,const size_type                   cols
	)
{
	namespace rcp = ReferenceCountingPack;
	namespace rmp = ResourceManagementPack;
	typedef rcp::ref_count_ptr<GenMatrix> vec_ptr_t;
	vec_ptr_t gms_ptr = rcp::rcp(new GenMatrix(rows,cols));
	this->initialize(
		(*gms_ptr)()
		,BLAS_Cpp::no_trans
		,rcp::rcp(
			new rmp::ReleaseResource_ref_count_ptr<GenMatrix>(
				gms_ptr
				)
			)
		);
}

void MultiVectorMutableDense::initialize(
	GenMatrixSlice                     gms
	,BLAS_Cpp::Transp                  gms_trans
	,const release_resource_ptr_t&     gms_release
	)
{
	gms_.bind(gms);
	gms_trans_   = gms_trans;
	gms_release_ = gms_release;
}

// Overridden from MatrixWithOpGetGMS

const GenMatrixSlice MultiVectorMutableDense::get_gms_view() const
{
	if(gms_trans_ == BLAS_Cpp::trans) {
		assert(0); // ToDo: We need to create a copy and transpose it!
	}
	return get_gms(); // No memory to allocate!
}

void MultiVectorMutableDense::free_gms_view(const GenMatrixSlice* gms_view) const
{
	if(gms_trans_ == BLAS_Cpp::trans) {
		assert(0); // ToDo: We need to free the copy that we created in get_gms_view()
	}
	else {
		// Nothing to free!
	}
}

// Overridden from MatrixWithOpGetGMSMutable

GenMatrixSlice MultiVectorMutableDense::get_gms_view()
{
	if(gms_trans_ == BLAS_Cpp::trans) {
		assert(0); // ToDo: We need to create a copy and transpose it!
	}
	return set_gms(); // No memory to allocate!
}

void MultiVectorMutableDense::commit_gms_view(GenMatrixSlice* gms_view)
{
	if(gms_trans_ == BLAS_Cpp::trans) {
		assert(0); // ToDo: We need to free the copy that we created in get_gms_view()
	}
	else {
		// Nothing to free!
	}
}

// Overridden from MultiVector

MultiVectorMutableDense::access_by_t
MultiVectorMutableDense::access_by() const
{
	return ROW_ACCESS & COL_ACCESS & DIAG_ACCESS; // row, column and diagonal access is available!
}

// Overridden from MultiVectorMutable

MultiVectorMutableDense::vec_mut_ptr_t
MultiVectorMutableDense::col(index_type j)
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp(
		new VectorWithOpMutableDense(
			LinAlgPack::col( set_gms(), gms_trans(), j )
			,NULL ) );
}

MultiVectorMutableDense::vec_mut_ptr_t
MultiVectorMutableDense::row(index_type i)
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp(
		new VectorWithOpMutableDense(
			LinAlgPack::row( set_gms(), gms_trans(), i )
			,NULL ) );
}

MultiVectorMutableDense::vec_mut_ptr_t
MultiVectorMutableDense::diag(int k)
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp(
		new VectorWithOpMutableDense(
			gms_.diag( gms_trans() == BLAS_Cpp::no_trans ? k : -k )
			,NULL ) );
}

MultiVectorMutableDense::multi_vec_mut_ptr_t
MultiVectorMutableDense::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng)
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp(
		new MultiVectorMutableDense(
			gms_(
				gms_trans() == BLAS_Cpp::no_trans   ? row_rng : col_rng
				,gms_trans() == BLAS_Cpp::no_trans  ? col_rng : row_rng )
			,gms_trans()
			,NULL ) );
}

// Overridden from MatrixBase

size_type MultiVectorMutableDense::rows() const
{
	return BLAS_Cpp::rows( get_gms().rows(), get_gms().cols(), gms_trans() );
}

size_type MultiVectorMutableDense::cols() const
{
	return BLAS_Cpp::cols( get_gms().rows(), get_gms().cols(), gms_trans() );
}

// Overridden from MatrixWithOp

void MultiVectorMutableDense::zero_out()
{
	gms_ = 0.0;
}

void MultiVectorMutableDense::Mt_S( value_type alpha )
{
	LinAlgPack::Mt_S(&gms_,alpha);
}

MatrixWithOp& MultiVectorMutableDense::operator=(const MatrixWithOp& mwo_rhs)
{
	LinAlgPack::assign( &set_gms(), MatrixDenseEncap(mwo_rhs)(), gms_trans() );
	return *this;
}

std::ostream& MultiVectorMutableDense::output(std::ostream& out) const
{
	if(gms_trans() == BLAS_Cpp::no_trans)
		return out << gms_;
	return MatrixWithOpSerial::output(out);
}

// Overridden from MatrixWithOpSerial

void MultiVectorMutableDense::Vp_StMtV(
	VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2, value_type beta) const
{
	LinAlgPack::Vp_StMtV(
		vs_lhs,alpha,gms_,BLAS_Cpp::trans_trans(gms_trans(),trans_rhs1),vs_rhs2,beta);
}

void MultiVectorMutableDense::Vp_StMtV(
	VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	SparseLinAlgPack::Vp_StMtV(
		vs_lhs,alpha,gms_,BLAS_Cpp::trans_trans(gms_trans(),trans_rhs1),sv_rhs2,beta);
}

// protected

// Overridden from MultiVector

void MultiVectorMutableDense::apply_op(
	EApplyBy apply_by, const RTOpPack::RTOp& primary_op
	,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
	,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
	,RTOp_ReductTarget reduct_objs[]
	,const index_type primary_first_ele   , const index_type primary_sub_dim   , const index_type primary_global_offset
	,const index_type secondary_first_ele , const index_type secondary_sub_dim 
	) const
{
	MultiVector::apply_op(
		apply_by, primary_op, num_multi_vecs, multi_vecs, num_targ_multi_vecs, targ_multi_vecs
		,reduct_objs
		,primary_first_ele, primary_sub_dim, primary_global_offset, secondary_first_ele, secondary_sub_dim
		); // ToDo: Provide Specialized implementation if needed?
}

void MultiVectorMutableDense::apply_op(
	EApplyBy apply_by, const RTOpPack::RTOp& primary_op, const RTOpPack::RTOp& secondary_op
	,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
	,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
	,RTOp_ReductTarget reduct_obj
	,const index_type primary_first_ele   , const index_type primary_sub_dim   , const index_type primary_global_offset
	,const index_type secondary_first_ele , const index_type secondary_sub_dim 
	) const
{
	MultiVector::apply_op(
		apply_by, primary_op, secondary_op, num_multi_vecs, multi_vecs, num_targ_multi_vecs, targ_multi_vecs
		,reduct_obj
		,primary_first_ele, primary_sub_dim, primary_global_offset, secondary_first_ele, secondary_sub_dim
		); // ToDo: Provide Specialized implementation if needed?
}

} // end namespace SparseLinAlgPack