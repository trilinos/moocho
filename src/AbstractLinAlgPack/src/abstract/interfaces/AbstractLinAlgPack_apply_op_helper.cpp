// /////////////////////////////////////////////////////////////////////////////
// apply_op_helper.cpp

#include <assert.h>

#include "apply_op_helper.hpp"
#include "VectorSpace.hpp"
#include "VectorMutable.hpp"
#include "WorkspacePack.hpp"
#include "Teuchos_TestForException.hpp"

void AbstractLinAlgPack::apply_op_validate_input(
	const char func_name[]
	,const RTOpPack::RTOp& op
	,const size_t num_vecs, const Vector* vecs[]
	,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
	,RTOp_ReductTarget reduct_obj
	,const index_type first_ele_in, const index_type sub_dim_in, const index_type global_offset_in
	)
{
	const VectorSpace
		&space = (num_vecs ? vecs[0]->space() : targ_vecs[0]->space() );
	const index_type
		dim = space.dim();
	TEST_FOR_EXCEPTION(
		global_offset_in < 0, std::logic_error
		,func_name << " : Error!  global_offset_in = "
		<<global_offset_in<<" is not valid" );
	TEST_FOR_EXCEPTION(
		first_ele_in > dim, std::logic_error
		,func_name << " : Error!  first_ele_in = "
		<<first_ele_in<<" is not compatible with space.dim() = " << dim );
	TEST_FOR_EXCEPTION(
		sub_dim_in < 0 || (sub_dim_in > 0 && sub_dim_in > dim-(first_ele_in-1)), std::logic_error
		,func_name << " : Error!  first_ele_in = "
		<<first_ele_in<<" and sub_dim_in = "<<sub_dim_in
		<<" is not compatible with space.dim() = " << dim );
	{for(int k = 0; k < num_vecs; ++k) {
		const bool is_compatible = space.is_compatible(vecs[k]->space());
		TEST_FOR_EXCEPTION(
			!is_compatible, VectorSpace::IncompatibleVectorSpaces
			,func_name << " : Error!  vecs["<<k<<"]->space() of type \'"
			<< typeid(vecs[k]->space()).name() << "\' "
			<< " with dimension vecs["<<k<<"].dim() = " << vecs[k]->dim()
			<< " is not compatible with space of type \'"
			<< typeid(space).name() << " with dimmension space.dim() = " << dim );
	}}
	{for(int k = 0; k < num_targ_vecs; ++k) {
		const bool is_compatible = space.is_compatible(targ_vecs[k]->space());
		TEST_FOR_EXCEPTION(
			!is_compatible, VectorSpace::IncompatibleVectorSpaces
			,func_name << " : Error!  targ_vecs["<<k<<"]->space() of type \'"
			<< typeid(targ_vecs[k]->space()).name() << "\' "
			<< " with dimension targ_vecs["<<k<<"].dim() = " << targ_vecs[k]->dim()
			<< " is not compatible with space of type \'"
			<< typeid(space).name() << " with dimmension space.dim() = " << dim );
	}}
}

void AbstractLinAlgPack::apply_op_serial(
	const RTOpPack::RTOp& op
	,const size_t num_vecs, const Vector* vecs[]
	,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
	,RTOp_ReductTarget reduct_obj
	,const index_type first_ele_in, const index_type sub_dim_in, const index_type global_offset_in
	)
{
 	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	
	// Dimension of global sub-vector
	const VectorSpace
		&space = ( num_vecs ? vecs[0]->space() : targ_vecs[0]->space() );
	const index_type
		full_dim       = space.dim(),
		global_sub_dim = sub_dim_in ? sub_dim_in     : full_dim - (first_ele_in-1);
	const Range1D
		global_sub_rng = Range1D(first_ele_in,(first_ele_in-1)+global_sub_dim);

	//
	// Get explicit views of the vector elements
	//

	wsp::Workspace<RTOpPack::SubVector >         local_vecs(wss,num_vecs);
	wsp::Workspace<RTOpPack::MutableSubVector>  local_targ_vecs(wss,num_targ_vecs);
	int k;
	for(k = 0; k < num_vecs; ++k) {
		RTOpPack::SubVector &v = local_vecs[k];
		vecs[k]->get_sub_vector( global_sub_rng, &v );
		v.setGlobalOffset( v.globalOffset() + global_offset_in );
	}
	for(k = 0; k < num_targ_vecs; ++k) {
		RTOpPack::MutableSubVector &v = local_targ_vecs[k];
		targ_vecs[k]->get_sub_vector( global_sub_rng, &v );
		v.setGlobalOffset( v.globalOffset() + global_offset_in );
	}

	//
	// Apply the reduction/transformation operator on all elements all at once!
	//

	op.apply_op(
		num_vecs,       num_vecs      ? &local_vecs[0]      : NULL
		,num_targ_vecs, num_targ_vecs ? &local_targ_vecs[0] : NULL
		,reduct_obj
		);

	//
	// Free (and commit) the explicit views of the vector elements
	// which should also inform the vectors that they have
	// changed.
	//

	for(k = 0; k < num_vecs; ++k) {
		RTOpPack::SubVector &v = local_vecs[k];
		v.setGlobalOffset( v.globalOffset() - global_offset_in );
		vecs[k]->free_sub_vector(&v);
	}
	for(k = 0; k < num_targ_vecs; ++k) {
		RTOpPack::MutableSubVector &v = local_targ_vecs[k];
		v.setGlobalOffset( v.globalOffset() - global_offset_in );
		targ_vecs[k]->commit_sub_vector(&v);
	}

}
