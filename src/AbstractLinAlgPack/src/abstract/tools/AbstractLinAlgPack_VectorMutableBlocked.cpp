// ///////////////////////////////////////////////////////////////////
// VectorMutableBlocked.cpp
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

#include <typeinfo>
#include <algorithm>

#include "AbstractLinAlgPack/src/abstract/tools/VectorMutableBlocked.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorMutableSubView.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/VectorSpaceBlocked.hpp"
#include "Range1D.hpp"
#include "WorkspacePack.hpp"
#include "ThrowException.hpp"

namespace {
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
template< class T >
inline
T my_min( const T& v1, const T& v2 ) { return v1 < v2 ? v1 : v2; }
} // end namespace

namespace AbstractLinAlgPack {

VectorMutableBlocked::VectorMutableBlocked(
	VectorMutable::vec_mut_ptr_t*  vecs
	,const vec_space_comp_ptr_t&         vec_space
	)
	// The members will be initialized in this->has_changed()
{
	initialize(vecs,vec_space);
}

void VectorMutableBlocked::initialize(
	VectorMutable::vec_mut_ptr_t*  vecs
	,const vec_space_comp_ptr_t&         vec_space
	)
{
	THROW_EXCEPTION(
		vec_space.get() == NULL, std::logic_error
		,"VectorMutableBlocked::initialize(...): Error, Must be constructed from "
		"a non-null block vector space!" );
	const int num_vec_spaces = vec_space->num_vector_spaces();
	vecs_.resize(num_vec_spaces);
	std::copy(vecs,vecs+num_vec_spaces,vecs_.begin());
	vec_space_ = vec_space;
	this->has_changed();
}
	
// overriddend form Vector

index_type VectorMutableBlocked::dim() const
{
	return vec_space_->dim();
}

const VectorSpace& VectorMutableBlocked::space() const
{
	return *vec_space_;
}

void VectorMutableBlocked::apply_op(
	const RTOpPack::RTOp& op
	,const size_t num_vecs, const Vector* vecs[]
	,const size_t num_targ_vecs, VectorMutable* targ_vecs[]
	,RTOp_ReductTarget reduct_obj
	,const index_type first_ele_in, const index_type sub_dim_in, const index_type global_offset_in
	) const
{
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	
	const index_type
		n = this->dim();

	// Validate the compatibility of the vectors!
#ifdef _DEBUG
	THROW_EXCEPTION(
		!(1 <= first_ele_in && first_ele_in <= n), std::out_of_range
		,"VectorMutableBlocked::apply_op(...): Error, "
		"first_ele_in = " << first_ele_in << " is not in range [1,"<<n<<"]" );
	THROW_EXCEPTION(
		sub_dim_in < 0, std::invalid_argument
		,"VectorMutableBlocked::apply_op(...): Error, "
		"sub_dim_in = " << sub_dim_in << " is not acceptable" );
	THROW_EXCEPTION(
		global_offset_in < 0, std::invalid_argument
		,"VectorMutableBlocked::apply_op(...): Error, "
		"global_offset_in = " << global_offset_in << " is not acceptable" );
	THROW_EXCEPTION(
		sub_dim_in > 0 && sub_dim_in - (first_ele_in - 1) > n, std::length_error
		,"VectorMutableBlocked::apply_op(...): Error, "
		"global_offset_in = " << global_offset_in << ", sub_dim_in = " << sub_dim_in
		<< "first_ele_in = " << first_ele_in << " and n = " << n
		<< " are not compatible" );
	bool test_failed;
	{for(int k = 0; k < num_vecs; ++k) {
		test_failed = !this->space().is_compatible(vecs[k]->space());
		THROW_EXCEPTION(
			test_failed, VectorSpace::IncompatibleVectorSpaces
			,"VectorMutableBlocked::apply_op(...): Error vecs["<<k<<"]->space() "
			<<"of type \'"<<typeid(vecs[k]->space()).name()<<"\' is not compatible with this "
			<<"\'VectorSpaceBlocked\' vector space!"
			);
	}}
	{for(int k = 0; k < num_targ_vecs; ++k) {
		test_failed = !this->space().is_compatible(targ_vecs[k]->space());
		THROW_EXCEPTION(
			test_failed, VectorSpace::IncompatibleVectorSpaces
			,"VectorMutableBlocked::apply_op(...): Error targ_vecs["<<k<<"]->space() "
			<<"of type \'"<<typeid(vecs[k]->space()).name()<<"\' is not compatible with this "
			<<"\'VectorSpaceBlocked\' vector space!"
			);
	}}
#endif
	
	// Dynamic cast the pointers for the vector arguments
	wsp::Workspace<const VectorMutableBlocked*>
		vecs_args(wss,num_vecs);
	{for(int k = 0; k < num_vecs; ++k) {
		vecs_args[k] = dynamic_cast<const VectorMutableBlocked*>(vecs[k]);
#ifdef _DEBUG
		THROW_EXCEPTION(
			vecs_args[k] == NULL, VectorSpace::IncompatibleVectorSpaces
			,"VectorMutableBlocked::apply_op(...): Error vecs["<<k<<"] "
			<<"of type \'"<<typeid(*vecs[k]).name()<<"\' does not support the "
			<<"\'VectorMutableBlocked\' interface!"
			);
#endif
	}}
	wsp::Workspace<VectorMutableBlocked*>
		targ_vecs_args(wss,num_targ_vecs);
	{for(int k = 0; k < num_targ_vecs; ++k) {
		targ_vecs_args[k] = dynamic_cast<VectorMutableBlocked*>(targ_vecs[k]);
#ifdef _DEBUG
		THROW_EXCEPTION(
			targ_vecs_args[k] == NULL, VectorSpace::IncompatibleVectorSpaces
			,"VectorMutableBlocked::apply_op(...): Error targ_vecs["<<k<<"] "
			<<"of type \'"<<typeid(*targ_vecs[k]).name()<<"\' does not support the "
			<<"\'VectorMutableBlocked\' interface!"
			);
#endif
	}}

	// Perform the reduction on each vector segment at a time.
	// ToDo: There could be a parallel version of this method!
	const index_type this_dim = this->dim();
	const index_type sub_dim  = ( sub_dim_in == 0
								  ? this_dim - (first_ele_in - 1)
								  : sub_dim_in );
	index_type num_elements_remaining = sub_dim;
	const int  num_vec_spaces = vec_space_->num_vector_spaces();
	wsp::Workspace<const Vector*>
		sub_vecs(wss,num_vecs);
	wsp::Workspace<VectorMutable*>
		sub_targ_vecs(wss,num_targ_vecs);
	index_type g_off = -first_ele_in + 1;
	for(int k = 0; k < num_vec_spaces; ++k) {
		const index_type local_dim = vecs_[k]->dim();
		if( g_off < 0 && -g_off > local_dim ) {
			g_off += local_dim;
			continue;
		}
		const index_type
			local_sub_dim = ( g_off >= 0
							  ? my_min( local_dim, num_elements_remaining )
							  : my_min( local_dim + g_off, num_elements_remaining ) );
		if( local_sub_dim <= 0 )
			break;
		for( int i = 0; i < num_vecs; ++i ) // Fill constituent vectors for segment k
			sub_vecs[i] = vecs_args[i]->vecs_[k].get();
		for( int j = 0; j < num_targ_vecs; ++j ) // Fill constituent target vectors for segment k
			sub_targ_vecs[j] = targ_vecs_args[j]->vecs_[k].get();
		AbstractLinAlgPack::apply_op( 
			op
			,num_vecs,num_vecs?&sub_vecs[0]:NULL
			,num_targ_vecs,num_targ_vecs?&sub_targ_vecs[0]:NULL
			,reduct_obj
			,g_off < 0 ? -g_off + 1 : 1                                // first_ele
			,local_sub_dim                                             // sub_dim
			,g_off < 0 ? global_offset_in : global_offset_in + g_off   // global_offset
			);
		g_off += local_dim;
		num_elements_remaining -= local_sub_dim;
	}
	assert(num_elements_remaining == 0);
	// Must allert all of the block vectors that they may have changed!
	{for(int k = 0; k < num_targ_vecs; ++k) {
		targ_vecs[k]->has_changed();
	}}
}

index_type VectorMutableBlocked::nz() const
{
	if( nz_ < 0.0 ) {
		const int num_vec_spaces = vec_space_->num_vector_spaces();
		nz_ = 0;
		for( int k = 0; k < num_vec_spaces; ++k )
			nz_ += vecs_[k]->nz();
	}
	return nz_;
}

std::ostream& VectorMutableBlocked::output(
	std::ostream& out, bool print_dim, bool newline
	,index_type global_offset
	) const
{
	if(print_dim)
		out << this->dim() << std::endl;
	size_type off = global_offset;
	const int num_vec_spaces = vec_space_->num_vector_spaces();
	for( int k = 0; k < num_vec_spaces; ++k ) {
		vecs_[k]->output(out,false,false,off);
		off += vecs_[k]->dim();
	}
	if(newline)
		out << std::endl;
	return out;
}

value_type VectorMutableBlocked::get_ele(index_type i) const
{
	int         kth_vector_space  = -1;
	index_type  kth_global_offset = 0;
	vec_space_->get_vector_space_position(i,&kth_vector_space,&kth_global_offset);
#ifdef _DEBUG
	assert( 0 <= kth_vector_space && kth_vector_space <= vecs_.size() );
#endif
	return vecs_[kth_vector_space]->get_ele( i - kth_global_offset );
}

value_type VectorMutableBlocked::norm_1() const
{
	if( norm_1_ < 0.0 ) {
		const int num_vec_spaces = vec_space_->num_vector_spaces();
		norm_1_ = 0.0;
		for( int k = 0; k < num_vec_spaces; ++k )
			norm_1_ += vecs_[k]->norm_1();
	}
	return norm_1_;
}

value_type VectorMutableBlocked::norm_inf() const
{
	if( norm_inf_ < 0.0 ) {
		const int num_vec_spaces = vec_space_->num_vector_spaces();
		norm_inf_ = 0.0;
		for( int k = 0; k < num_vec_spaces; ++k )
			norm_inf_  = my_max( norm_inf_, vecs_[k]->norm_inf() );
	}
	return norm_inf_;
}

value_type VectorMutableBlocked::inner_product( const Vector& v ) const
{
	return Vector::inner_product(v); // ToDo: Specialize? (why bother?)
}

void VectorMutableBlocked::get_sub_vector( const Range1D& rng_in, RTOpPack::SubVector* sub_vec ) const
{
	const Range1D
		rng = rng_in.full_range() ? Range1D( 1, this->dim()) : rng_in;
	int         kth_vector_space  = -1;
	index_type  kth_global_offset = 0;
	vec_space_->get_vector_space_position(rng.lbound(),&kth_vector_space,&kth_global_offset);
#ifdef _DEBUG
	assert( 0 <= kth_vector_space && kth_vector_space <= vecs_.size() );
#endif
	if( rng.lbound() + rng.size() <= kth_global_offset + 1 + vecs_[kth_vector_space]->dim() ) {
		// This involves only one sub-vector so just return it.
		static_cast<const Vector*>(vecs_[kth_vector_space].get())->get_sub_vector(
			rng - kth_global_offset, sub_vec );
		sub_vec->setGlobalOffset( sub_vec->globalOffset() + kth_global_offset );
	}
	else {
		// Just let the default implementation handle this.  ToDo: In the futrue
		// we could manually construct an explicit sub-vector that spanned
		// two or more consitituent vectors but this would be a lot of work.
		// However, this would require the use of temporary memory but
		// so what.
		Vector::get_sub_vector(rng_in,sub_vec);
	}
}

void VectorMutableBlocked::free_sub_vector( RTOpPack::SubVector* sub_vec ) const
{
	if( sub_vec->values() == NULL ) return;
	int         kth_vector_space  = -1;
	index_type  kth_global_offset = 0;
	vec_space_->get_vector_space_position(
		sub_vec->globalOffset()+1,&kth_vector_space,&kth_global_offset);
#ifdef _DEBUG
	assert( 0 <= kth_vector_space && kth_vector_space <= vecs_.size() );
#endif
	if( sub_vec->globalOffset() + sub_vec->subDim() <= kth_global_offset +  vecs_[kth_vector_space]->dim() ) {
		// This sub_vec was extracted from a single constituent vector
		sub_vec->setGlobalOffset( sub_vec->globalOffset() - kth_global_offset );
		vecs_[kth_vector_space].get()->free_sub_vector(sub_vec);
	}
	else {
		// This sub_vec was created by the default implementation!
		Vector::free_sub_vector(sub_vec);
	}
}

void VectorMutableBlocked::has_changed() const
{
	nz_ = -1;                   // set to uninitalized!
	norm_1_ = norm_inf_ = -1.0;
	Vector::has_changed();
}

// overridden from VectorMutable

VectorMutable::vec_mut_ptr_t
VectorMutableBlocked::sub_view( const Range1D& rng_in )
{
	namespace rcp = MemMngPack;
	const index_type dim = this->dim();
	const Range1D    rng = rng_in.full_range() ? Range1D(1,dim) : rng_in;
	// Validate the preconditions
#ifdef _DEBUG
	THROW_EXCEPTION(
		dim < rng.ubound(), std::out_of_range
		,"VectorMutableBlocked::sub_view(...): Error, rng = "
		<< "["<<rng.lbound()<<","<<rng.ubound()<<"] is not in range [1,"<<dim<<"]" );
#endif
	vecs_t &vecs = this->vecs_; // Need to examine in debugger!
	if( rng.lbound() == 1 && rng.ubound() == dim ) {
		if( vecs.size() == 1 )   return vecs[0]->sub_view(Range1D());
		else                     return rcp::rcp(this,false);
	}
	// Get the position of the vector space object of interest
	int           kth_vector_space  = -1;
	index_type    kth_global_offset = 0;
	vec_space_->get_vector_space_position(rng.lbound(),&kth_vector_space,&kth_global_offset);
	const VectorSpace::space_ptr_t*  vector_spaces      = vec_space_->vector_spaces();
	const index_type*                vec_spaces_offsets = vec_space_->vector_spaces_offsets();
#ifdef _DEBUG
	assert( 0 <= kth_vector_space && kth_vector_space <= vecs.size() );
#endif
	if( rng.lbound() == kth_global_offset + 1
		&& rng.size() == vec_spaces_offsets[kth_vector_space+1] - vec_spaces_offsets[kth_vector_space] )
		// The client selected a whole single constituent vector.
		return vecs[kth_vector_space]->sub_view(Range1D());
	if( rng.ubound() <= vec_spaces_offsets[kth_vector_space+1] )
		// The client selected a sub-vector of a single consituent vector
		return vecs[kth_vector_space]->sub_view( rng - vec_spaces_offsets[kth_vector_space] );
	// The client selected a sub-vector that spans two or more constituent vectors.
	// Get the position of the vector space object with the last element of interest
	int           end_kth_vector_space  = -1;
	index_type    end_kth_global_offset = 0;
	vec_space_->get_vector_space_position(rng.ubound(),&end_kth_vector_space,&end_kth_global_offset);
#ifdef _DEBUG
	assert( 0 <= end_kth_vector_space && end_kth_vector_space <= vecs.size() );
	assert( end_kth_vector_space > kth_vector_space );
#endif
	// Create a VectorWithOpMutableCompsiteStd object containing the relavant constituent vectors
	rcp::ref_count_ptr<VectorMutableBlocked>
		vec_comp = rcp::rcp(new VectorMutableBlocked(
			&vecs[kth_vector_space]
			,rcp::rcp(new VectorSpaceBlocked(
				&vector_spaces[kth_vector_space]
				,end_kth_vector_space - kth_vector_space + 1 ))
			));
	if( rng.lbound() == kth_global_offset + 1
		&& rng.size() == vec_spaces_offsets[end_kth_vector_space+1] - vec_spaces_offsets[kth_vector_space] )
		// The client selected exactly a contigous set of vectors
		return vec_comp;
	// The client selected some sub-set of elements in the contigous set of vectors
	return rcp::rcp(
		new VectorMutableSubView(
			vec_comp
			,Range1D( 
				rng.lbound()-vec_spaces_offsets[kth_vector_space]
				,rng.ubound()-vec_spaces_offsets[kth_vector_space] )
			) );
}

void VectorMutableBlocked::axpy( value_type alpha, const Vector& x )
{
	VectorMutable::axpy(alpha,x); // ToDo: Specialize? (why bother?)
}

VectorMutable& VectorMutableBlocked::operator=(value_type alpha)
{
	const int num_vec_spaces = vec_space_->num_vector_spaces();
	for( int k = 0; k < num_vec_spaces; ++k )
		*vecs_[k] = alpha;
	return *this;
}

VectorMutable& VectorMutableBlocked::operator=(const Vector& v)
{
	VectorMutable::operator=(v); // ToDo: Specialize (why bother?)
	return *this;
}

void VectorMutableBlocked::set_ele( index_type i, value_type val )
{
	int         kth_vector_space  = -1;
	index_type  kth_global_offset = 0;
	vec_space_->get_vector_space_position(i,&kth_vector_space,&kth_global_offset);
#ifdef _DEBUG
	assert( 0 <= kth_vector_space && kth_vector_space <= vecs_.size() );
#endif
	vecs_[kth_vector_space]->set_ele( i - kth_global_offset, val );
}

void VectorMutableBlocked::set_sub_vector( const RTOpPack::SparseSubVector& sub_vec )
{
	int         kth_vector_space  = -1;
	index_type  kth_global_offset = 0;
	vec_space_->get_vector_space_position(
		sub_vec.globalOffset()+1,&kth_vector_space,&kth_global_offset);
#ifdef _DEBUG
	assert( 0 <= kth_vector_space && kth_vector_space <= vecs_.size() );
#endif
	if( sub_vec.globalOffset() + sub_vec.subDim() <= kth_global_offset +  vecs_[kth_vector_space]->dim() ) {
		// This sub-vector fits into a single constituent vector
		RTOpPack::SparseSubVector sub_vec_g = sub_vec;
		sub_vec_g.setGlobalOffset( sub_vec_g.globalOffset() - kth_global_offset );
		vecs_[kth_vector_space]->set_sub_vector(sub_vec_g);
	}
	else {
		// Let the default implementation take care of this.  ToDo: In the futrue
		// it would be possible to manualy set the relavent constituent
		// vectors with no temp memory allocations.
		VectorMutable::set_sub_vector(sub_vec);
	}
}

void VectorMutableBlocked::assert_in_range(int k) const
{
	assert_initialized();
	THROW_EXCEPTION(
		k >= vec_space_->num_vector_spaces(), std::logic_error
		,"VectorMutableBlocked::assert_in_range(k) : Error, k = " << k << " is not "
		"in range [0," << vec_space_->num_vector_spaces() << "]" );
}

void VectorMutableBlocked::assert_initialized() const
{
	THROW_EXCEPTION(!vec_space_.get(),std::logic_error,"Error, this is not initalized!");
}

} // end namespace AbstractLinAlgPack