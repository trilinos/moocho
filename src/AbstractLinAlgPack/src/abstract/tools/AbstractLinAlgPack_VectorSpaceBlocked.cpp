// //////////////////////////////////////////////////////////////////////
// VectorSpaceBlocked.cpp
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

#include <algorithm>

#include "AbstractLinAlgPack/src/abstract/tools/VectorSpaceBlocked.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/VectorMutableBlocked.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorSpaceSubSpace.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorSpaceFactory.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/GenPermMatrixSlice.hpp"
#include "Range1D.hpp"
#include "WorkspacePack.hpp"
#include "ThrowException.hpp"

namespace AbstractLinAlgPack {

VectorSpaceBlocked::VectorSpaceBlocked(
	const VectorSpace::space_ptr_t         vec_spaces[]
	,int                                   num_vec_spaces
	,const VectorSpace::space_fcty_ptr_t   &small_vec_spc_fcty
	)
{
	initialize(vec_spaces,num_vec_spaces,small_vec_spc_fcty);
}
	
void VectorSpaceBlocked::initialize(
	const VectorSpace::space_ptr_t         vec_spaces[]
	,int                                   num_vec_spaces
	,const VectorSpace::space_fcty_ptr_t   &small_vec_spc_fcty
	)
{
	vector_spaces_.resize(num_vec_spaces);
	std::copy(vec_spaces,vec_spaces+num_vec_spaces,vector_spaces_.begin());
	vec_spaces_offsets_.resize(num_vec_spaces+1);
	vec_spaces_offsets_[0] = 0;
	for( int k = 1; k <= num_vec_spaces; ++k )
		vec_spaces_offsets_[k] = vec_spaces_offsets_[k-1] + vec_spaces[k-1]->dim();
	small_vec_spc_fcty_ = small_vec_spc_fcty;
}

void VectorSpaceBlocked::get_vector_space_position(
	index_type i, int* kth_vector_space, index_type* kth_global_offset ) const
{
	// Validate the preconditions
#ifdef _DEBUG
	THROW_EXCEPTION(
		i < 1 || this->dim() < i, std::out_of_range
		,"VectorSpaceBlocked::get_vector_space_position(...): Error, i = "
		<< i << " is not in range [1,"<<this->dim()<<"]"
		);
#endif
	*kth_vector_space  = 0;
	*kth_global_offset = 0;
	while( *kth_vector_space < vector_spaces_.size() ) {
		const RTOp_index_type off_kp1 = vec_spaces_offsets_[*kth_vector_space+1];
		if( off_kp1 + 1 > i ) {
			*kth_global_offset = vec_spaces_offsets_[*kth_vector_space];
			break;
		}
		++(*kth_vector_space);
	}
#ifdef _DEBUG
	assert(*kth_vector_space < vector_spaces_.size());
#endif
}

// overridden from VectorSpace

bool VectorSpaceBlocked::is_compatible(const VectorSpace& vec_space) const
{
	const VectorSpaceBlocked
		*vec_space_comp = dynamic_cast<const VectorSpaceBlocked*>(&vec_space);
	if( !vec_space_comp || vec_space_comp->vector_spaces_.size() != this->vector_spaces_.size() )
		return false;
	// There is some hope that these are compatible vector spaces.
	for( int k = 0; k < vector_spaces_.size(); ++k ) {
		if( !vec_space_comp->vector_spaces_[k]->is_compatible(*this->vector_spaces_[k]) )
			return false;
	}
	return true; // If we get here we are compatible!
}
	
index_type VectorSpaceBlocked::dim() const
{
	return vec_spaces_offsets_[vector_spaces_.size()];
}

VectorSpace::space_fcty_ptr_t
VectorSpaceBlocked::small_vec_spc_fcty() const
{
	return small_vec_spc_fcty_;
}

VectorSpace::vec_mut_ptr_t
VectorSpaceBlocked::create_member() const
{
	namespace rcp = MemMngPack;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	const int num_vec_spaces = this->num_vector_spaces();

	// Create the vector objects array.
	wsp::Workspace<VectorMutable::vec_mut_ptr_t>
		vecs(wss,num_vec_spaces);
	for( int k = 0; k < num_vec_spaces; ++k )
		vecs[k] = vector_spaces_[k]->create_member();

	return rcp::rcp(
		new VectorMutableBlocked(
			&vecs[0]
			,rcp::rcp(new VectorSpaceBlocked(*this))
			) );
}

VectorSpace::multi_vec_mut_ptr_t
VectorSpaceBlocked::create_members(size_type num_vecs) const
{
	assert(0); // ToDo: Implement using MultiVectorMutableCompositeStd!
	return MemMngPack::null;
}

VectorSpace::space_ptr_t
VectorSpaceBlocked::clone() const
{
	return MemMngPack::rcp(new VectorSpaceBlocked(*this)); // ToDo: Fix the behavior when needed!
}

VectorSpace::space_ptr_t
VectorSpaceBlocked::sub_space(const Range1D& rng_in) const
{
	namespace rcp = MemMngPack;
	const index_type dim = this->dim();
	const Range1D    rng = rng_in.full_range() ? Range1D(1,dim) : rng_in;
	// Validate the preconditions
#ifdef _DEBUG
	THROW_EXCEPTION(
		dim < rng.ubound(), std::out_of_range
		,"VectorSpaceBlocked::sub_space(...): Error, rng = "
		<< "["<<rng.lbound()<<","<<rng.ubound()<<"] is not in range [1,"<<dim<<"]" );
#endif
	if( rng.lbound() == 1 && rng.ubound() == dim )
		return space_ptr_t( this, false ); // Client selected the whole composite vector space.
	// Get the position of the vector space object of interest
	int           kth_vector_space  = -1;
	index_type    kth_global_offset = 0;
	this->get_vector_space_position(rng.lbound(),&kth_vector_space,&kth_global_offset);
	const vector_spaces_t      &vector_spaces      = vector_spaces_;      // Need to examine in debugger!
	const vec_spaces_offsets_t &vec_spaces_offsets = vec_spaces_offsets_;
#ifdef _DEBUG
	assert( 0 <= kth_vector_space && kth_vector_space <= vector_spaces.size() );
#endif
	if( rng.lbound() == kth_global_offset + 1
		&& rng.size() == vec_spaces_offsets[kth_vector_space+1] - vec_spaces_offsets[kth_vector_space] )
		// The client selected a whole single constituent vector space.
		return vector_spaces[kth_vector_space];
	if( rng.ubound() <= vec_spaces_offsets[kth_vector_space+1] )
		// The client selected a sub-space of a single consituent vector space
		return vector_spaces[kth_vector_space]->sub_space(rng-vec_spaces_offsets[kth_vector_space]);
	// The client selected a sub-space that spans two or more constituent vector spaces
	// Get the position of the vector space object with the last element of interest
	int           end_kth_vector_space  = -1;
	index_type    end_kth_global_offset = 0;
	this->get_vector_space_position(rng.ubound(),&end_kth_vector_space,&end_kth_global_offset);
#ifdef _DEBUG
	assert( 0 <= end_kth_vector_space && end_kth_vector_space <= vector_spaces.size() );
	assert( end_kth_vector_space > kth_vector_space );
#endif
	// Create a VectorSpaceComposite object containing the relavant constituent vector spaces
	rcp::ref_count_ptr<VectorSpaceBlocked>
		vec_space_comp = rcp::rcp(
			new VectorSpaceBlocked(
				&vector_spaces[kth_vector_space]
				,end_kth_vector_space - kth_vector_space + 1 )
			);
	if( rng.lbound() == kth_global_offset + 1
		&& rng.size() == vec_spaces_offsets[end_kth_vector_space+1] - vec_spaces_offsets[kth_vector_space] )
		// The client selected exactly a contigous set of vector spaces
		return vec_space_comp;
	// The client selected some sub-set of elements in the contigous set of vector spaces
	return rcp::rcp(
		new VectorSpaceSubSpace(
			vec_space_comp
			,Range1D( 
				rng.lbound()-vec_spaces_offsets[kth_vector_space]
				,rng.ubound()-vec_spaces_offsets[kth_vector_space] )
			) );
}

} // end namespace AbstractLinAlgPack