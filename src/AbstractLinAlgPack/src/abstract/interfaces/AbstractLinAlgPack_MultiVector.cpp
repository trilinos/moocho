// /////////////////////////////////////////////////////////////////////////
// MultiVector.cpp

#include <assert.h>

#include "AbstractLinAlgPack/include/MultiVectorMutable.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "WorkspacePack.h"
#include "ThrowException.h"

namespace {

// Map from EApplyBy to Transp
inline
BLAS_Cpp::Transp
to_trans(AbstractLinAlgPack::MultiVector::EApplyBy apply_by)
{
	return ( apply_by == AbstractLinAlgPack::MultiVector::APPLY_BY_ROW
			? BLAS_Cpp::no_trans
			: BLAS_Cpp::trans
			);
}

// Return a row or a column vector from a multi-vector

inline 
AbstractLinAlgPack::MultiVector::vec_ptr_t
vec(
	const AbstractLinAlgPack::MultiVector&      multi_vec
	,const AbstractLinAlgPack::size_type        k
	,AbstractLinAlgPack::MultiVector::EApplyBy  apply_by
	)
{
	return ( apply_by == AbstractLinAlgPack::MultiVector::APPLY_BY_ROW
			? multi_vec.row(k)
			: multi_vec.col(k)
			);
}

inline 
AbstractLinAlgPack::MultiVectorMutable::vec_mut_ptr_t
vec(
	AbstractLinAlgPack::MultiVectorMutable*         multi_vec
	,const AbstractLinAlgPack::size_type            k
	,AbstractLinAlgPack::MultiVector::EApplyBy      apply_by
	)
{
	return ( apply_by == AbstractLinAlgPack::MultiVector::APPLY_BY_ROW
			? multi_vec->row(k)
			: multi_vec->col(k)
			);
}

} // end namespace

namespace AbstractLinAlgPack {

MultiVector::multi_vec_ptr_t
MultiVector::mv_sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
	assert(0); // ToDo: return a MultiVectorSubView object.
	// Note that the MultiVectorSubView class should derive from MatrixWithOpSubView
	// so that a client can rely on the MatrixWithOpSubView interface.
	return NULL;
}

void MultiVector::apply_reduction(
	EApplyBy apply_by, const RTOpPack::RTOp& prim_op
	,const size_t num_multi_vecs_in,   const MultiVector**   multi_vecs_in
	,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
	,RTOp_ReductTarget reduct_objs[]
	,const index_type prim_first_ele_in, const index_type prim_sub_dim_in, const index_type prim_global_offset_in
	,const index_type sec_first_ele_in, const index_type sec_sub_dim_in
	) const
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	// ToDo: Validate the input!

	wsp::Workspace<const MultiVector*> multi_vecs(wss,num_multi_vecs_in+1);
	multi_vecs[0] = this; // I am the first!
	for(size_type k = 1; k <= num_multi_vecs_in; ++k)
		multi_vecs[k] = multi_vecs_in[k-1];

	this->apply_op(
		apply_by, prim_op
		,num_multi_vecs_in+1,  &multi_vecs[0]
		,num_targ_multi_vecs,  targ_multi_vecs
		,reduct_objs
		,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
		,sec_first_ele_in, sec_sub_dim_in
		);
}

void MultiVector::apply_reduction(
	EApplyBy apply_by, const RTOpPack::RTOp& prim_op, const RTOpPack::RTOp& sec_op
	,const size_t num_multi_vecs_in,   const MultiVector**   multi_vecs_in
	,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
	,RTOp_ReductTarget reduct_obj
	,const index_type prim_first_ele_in, const index_type prim_sub_dim_in, const index_type prim_global_offset_in
	,const index_type sec_first_ele_in, const index_type sec_sub_dim_in
	) const
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	// ToDo: Validate the input!

	wsp::Workspace<const MultiVector*> multi_vecs(wss,num_multi_vecs_in+1);
	multi_vecs[0] = this; // I am the first!
	for(size_type k = 1; k <= num_multi_vecs_in; ++k)
		multi_vecs[k] = multi_vecs_in[k-1];

	this->apply_op(
		apply_by, prim_op, sec_op
		,num_multi_vecs_in+1,  &multi_vecs[0]
		,num_targ_multi_vecs,  targ_multi_vecs
		,reduct_obj
		,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
		,sec_first_ele_in, sec_sub_dim_in
		);
}

// Overridden form MatrixWithOp

MatrixWithOp::mat_ptr_t
MultiVector::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
	return mv_sub_view(row_rng,col_rng);
}

// Combined implementations for apply_reduction() and apply_transformation() methods

void MultiVector::apply_op(
	EApplyBy apply_by, const RTOpPack::RTOp& prim_op
	,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
	,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
	,RTOp_ReductTarget reduct_objs[]
	,const index_type prim_first_ele_in, const index_type prim_sub_dim_in, const index_type prim_global_offset_in
	,const index_type sec_first_ele_in, const index_type sec_sub_dim_in
	) const
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	// ToDo: Validate the input!

	// Get the primary and secondary dimmensions.
	const index_type
		prim_dim     = ( apply_by == APPLY_BY_ROW ? rows()          : cols()   ),
		sec_dim      = ( apply_by == APPLY_BY_ROW ? cols()          : rows()   ),
		prim_sub_dim = ( prim_sub_dim_in != 0     ? prim_sub_dim_in : prim_dim ),
		sec_sub_dim  = ( sec_sub_dim_in != 0      ? sec_sub_dim_in  : sec_dim  );
	assert(0 < prim_sub_dim && prim_sub_dim < prim_dim );
	assert(0 < sec_sub_dim && sec_sub_dim < sec_dim );

	//
	// Apply the reduction/transformation operator and trnasform the target
	// vectors and reduce each of the reduction objects.
	//

	wsp::Workspace<MultiVector::vec_ptr_t>             vecs_s(wss,num_multi_vecs);
	wsp::Workspace<const VectorWithOp*>                vecs(wss,num_multi_vecs);
	wsp::Workspace<MultiVectorMutable::vec_mut_ptr_t>  targ_vecs_s(wss,num_targ_multi_vecs);
	wsp::Workspace<VectorWithOpMutable*>               targ_vecs(wss,num_multi_vecs);

	{for(size_type j = sec_first_ele_in; j <= sec_first_ele_in - 1 + sec_sub_dim; ++j) {
		// Fill the arrays of vector arguments 
		{for(size_type k = 0; k < num_multi_vecs; ++k) {
			vecs_s[k] = vec( *multi_vecs[k], j, apply_by );
			vecs[k] = vecs_s[k].get();
		}}
		{for(size_type k = 0; k < num_targ_multi_vecs; ++k) {
			targ_vecs_s[k] = vec( targ_multi_vecs[k], j, apply_by );
			targ_vecs[k] = targ_vecs_s[k].get();
		}}
		// Apply the reduction/transformation operator
		if( num_multi_vecs > 0 )
			vecs[0]->apply_reduction(
				prim_op
				,num_multi_vecs-1,    (num_multi_vecs-1    ? &vecs[1]      : NULL)
				,num_targ_multi_vecs, (num_targ_multi_vecs ? &targ_vecs[0] : NULL)
				,reduct_objs[j - sec_first_ele_in]
				,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
			);
		else
			targ_vecs[0]->apply_transformation(
				prim_op
				,num_multi_vecs,        (num_multi_vecs        ? &vecs[0]      : NULL)
				,num_targ_multi_vecs-1, (num_targ_multi_vecs-1 ? &targ_vecs[1] : NULL)
				,reduct_objs[j - sec_first_ele_in]
				,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
			);
	}}

	// At this point all of the designated targ vectors in the target multi-vectors have
	// been transformed and all the reduction objects in reduct_obj[] have accumulated
	// the reductions.

}

void MultiVector::apply_op(
	EApplyBy apply_by, const RTOpPack::RTOp& prim_op, const RTOpPack::RTOp& sec_op
	,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
	,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
	,RTOp_ReductTarget reduct_obj
	,const index_type prim_first_ele_in, const index_type prim_sub_dim_in, const index_type prim_global_offset_in
	,const index_type sec_first_ele_in, const index_type sec_sub_dim_in
	) const
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	// ToDo: Validate the input!

	// Get the primary and secondary dimmensions.
	const index_type
		prim_dim    = ( apply_by == APPLY_BY_ROW ? rows()         : cols()  ),
		sec_dim     = ( apply_by == APPLY_BY_ROW ? cols()         : rows()  ),
		sec_sub_dim = ( sec_sub_dim_in != 0      ? sec_sub_dim_in : sec_dim );
	assert(0 < sec_sub_dim && sec_sub_dim < sec_dim );

	// Create a temporary buffer for the reduction objects of the primary reduction
	// so that we can call the companion version of this method.
	wsp::Workspace<RTOp_ReductTarget>   reduct_objs(wss,sec_sub_dim);
	{for(index_type k = 0; k < sec_sub_dim; ++k) {
		prim_op.reduct_obj_create_raw( &(reduct_objs[k]=RTOp_REDUCT_OBJ_NULL) );
	}}
	
	// Call the campanion version that accepts an array of reduction objects
	this->apply_op(
		apply_by, prim_op
		,num_multi_vecs,       multi_vecs
		,num_targ_multi_vecs,  targ_multi_vecs
		,&reduct_objs[0]
		,prim_first_ele_in, prim_sub_dim_in, prim_global_offset_in
		,sec_first_ele_in,  sec_sub_dim_in
		);

	// Reduce all the reduction objects using the secondary reduction operator
	// into one reduction object and free the intermedate reduction objects.
	{for(index_type k = 0; k < sec_sub_dim; ++k) {
		sec_op.reduce_reduct_objs( reduct_objs[k] ,reduct_obj );
		prim_op.reduct_obj_free( &reduct_objs[k] );
	}}
}

} // end namespace