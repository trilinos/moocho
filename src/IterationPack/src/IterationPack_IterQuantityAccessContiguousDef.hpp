// /////////////////////////////////////////////////////////////////////////////////////
// IterQuantityAccessContiguousDef.h
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
//
// Definitions to template functions

#ifndef ITER_QUANITY_ACCESS_CONTINUOUS_DEF_H
#define ITER_QUANITY_ACCESS_CONTINUOUS_DEF_H

#include <typeinfo>
#include <algorithm>
#include <iterator>

#include "IterQuantityAccessContiguousDecl.h"
#include "ThrowException.h"

namespace GeneralIterationPack {

// Constructors/initializers

template<class T_info>
IterQuantityAccessContiguous<T_info>::IterQuantityAccessContiguous(
	int                              num_quantities
	,const std::string&              name
	,const abstract_factory_ptr_t&   abstract_factory
	)
	:num_quantities_(0)
	,name_(name)
	,abstract_factory_(abstract_factory)
	,max_offset_(0)
{
	resize( num_quantities );
}

template<class T_info>
IterQuantityAccessContiguous<T_info>::~IterQuantityAccessContiguous() {
	release_mem();
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::set_factory(
	const abstract_factory_ptr_t& abstract_factory
	)
{
	release_mem();
	max_offset_ = std::numeric_limits<int>::min() + num_quantities_;  // uninitialized
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::resize( int num_quantities ) {
	THROW_EXCEPTION(
		num_quantities < 1, std::length_error
		,"IterQuantityAccessContiguous::resize(num_quantities): Error, "
		"name = "<<name_<<", num_quantities = "<<num_quantities<<" must be greater than zero" );
	if( num_quantities_ != num_quantities )
		release_mem();
	num_quantities_ = num_quantities;
	max_offset_ = std::numeric_limits<int>::min() + num_quantities_; // uninitialized
}

// Overridden from IterQuantity

template<class T_info>
IterQuantity* IterQuantityAccessContiguous<T_info>::clone() const {
	return 0;
	// ToDo: replace above with the following when the copy
	// constructor is implemented.
	// return new IterQuantityAccessContiguous<T_info>(*this);
}

template<class T_info>
const char* IterQuantityAccessContiguous<T_info>::name() const {
	return name_.c_str();
}

template<class T_info>
bool IterQuantityAccessContiguous<T_info>::has_storage_k(int offset) const {
	return is_initialized()
		? offset >= max_offset_ - num_quantities_ + 1
		: true;
}

template<class T_info>
bool IterQuantityAccessContiguous<T_info>::updated_k(int offset) const {
	if( !is_initialized() )
		return false;
	if( ( offset > max_offset_ ) || ( offset < max_offset_ - num_quantities_ + 1 ) )
		return false;
	return updated_[max_offset_ - offset];
}

template<class T_info>
int IterQuantityAccessContiguous<T_info>::last_updated() const {
	if( !is_initialized() )
		return NONE_UPDATED;
	// Find the last still set as updated.
	for(	int offset = max_offset_;
			offset >= max_offset_ - num_quantities_ + 1;
			--offset										)
	{
		if( updated_[max_offset_ - offset] )
			return offset;
	}
	return NONE_UPDATED;
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::set_not_updated(int offset) {
	assert_updated_k(offset);
	updated_[max_offset_ - offset] = false;
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::set_all_not_updated() {
	if(!is_initialized()) return;
	std::fill( updated_.begin(), updated_.end(), false );
}

template<class T_info>
bool IterQuantityAccessContiguous<T_info>::will_loose_mem(int offset, int set_offset) const {
	assert_updated_k(offset);
	return set_offset - max_offset_ > num_quantities_ - (max_offset_ - offset) - 1;
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::next_iteration()
{
	if( !is_initialized() ) return;
	--max_offset_;
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::print_concrete_type( std::ostream& out ) const
{
	const int last_updated = this->last_updated();
	if(last_updated != NONE_UPDATED)
		out << typeid(get_k(last_updated)).name();
	else if( abstract_factory_.get() == NULL )
		out << "NULL";
	else
		out << typeid(*abstract_factory_->create()).name();
}

// Overridden from IterQuantityAccess

template<class T_info>
T_info& IterQuantityAccessContiguous<T_info>::get_k(int offset) {
	assert_updated_k(offset);
	return *quantities_[max_offset_ - offset];
}
template<class T_info>
const T_info& IterQuantityAccessContiguous<T_info>::get_k(int offset) const {
	assert_updated_k(offset);
	return *quantities_[max_offset_ - offset];
}

template<class T_info>
T_info& IterQuantityAccessContiguous<T_info>::set_k(int offset) {

	lazy_initialization();

	assert_has_storage_k(offset);	// assert that we are not trying to iterate backwards

	if(offset > max_offset_ + num_quantities_ - 1) {
		// There will be no back memory so you don't need to adjust the pointers
		max_offset_ = offset;
		std::fill(updated_.begin(), updated_.end(), false);
	}
	else {
		// Pointers may have to be rearranged
		if(offset > max_offset_) {
			// We need to rearrange quantities_ and updated_.
			int shifted = offset - max_offset_;

			// ///////////////////////////////////
			// Set the updated flags

			// Example: We are shifting from:
			//		[1, 0, -1, -2] to [ 3, 2, 1, 0 ]

			// Shift the flags for the memory that may be saved
			updated_t::iterator
				itr_updated_from        = updated_.begin(),
				itr_updated_from_end    = itr_updated_from + num_quantities_ - shifted,
				itr_updated_to          = itr_updated_from + shifted;

			std::copy(itr_updated_from, itr_updated_from_end, itr_updated_to);

			// make updated[] for the new quantities false
			std::fill_n( updated_.begin(), shifted, false );

			// /////////////////////////////////////
			// rotate the quantitiy pointer vector

			// example: [1, 0, -1, -2] => [3, 2, 1, 0] so reverse rotate 0 to the back.

			// Rotate the (num_quantities_ - shifted) pointer to the back and the
			// remainder back arround again.

#if defined(_INTEL_CXX)
			typedef std::reverse_iterator<T_info**, T_info*, T_info*&
				, T_info**, ptrdiff_t>									rev_t;
#else
			typedef std::reverse_iterator<T_info**>						rev_t;
#endif

			std::rotate(
				  rev_t(&quantities_[0] + num_quantities_)
				, rev_t(&quantities_[0] + num_quantities_ - shifted)
				, rev_t(&quantities_[0])                              );

			max_offset_ = offset;

		}
		// else, no pointers need to be rearranged since we are not advancing the range
	}

	updated_[max_offset_ - offset] = true;
	return *quantities_[max_offset_ - offset];
}

// Private member functions

template<class T_info>
bool IterQuantityAccessContiguous<T_info>::is_initialized() const {
	return store_.size() > 0;
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::lazy_initialization() {
	if( !is_initialized() ) {
		THROW_EXCEPTION(
			abstract_factory_.get() == NULL, std::logic_error
			,"IterQuantityAccessContiguous::lazy_initialization(): Error, "
			"iq_name = "<<name_<<" the abstract factory can not be NULL" );
		// Allocate storage
		updated_.resize(num_quantities_,false);
		store_.resize(num_quantities_);
		quantities_.resize(num_quantities_,NULL);
		// Set initial points to locations
		updated_t::iterator       itr_updated         = updated_.begin();
		store_t::iterator         itr_store           = store_.begin();
		quantities_t::iterator    itr_quantities      = quantities_.begin();
		for( ; itr_store != store_.end(); ++itr_updated, ++itr_store, ++itr_quantities )
		{
			*itr_updated     = false;
			*itr_store       = abstract_factory_->create();
			*itr_quantities  = itr_store->get();
		}
		max_offset_ = std::numeric_limits<int>::min() + num_quantities_ + 1;
	}
}

template<class T_info>
void IterQuantityAccessContiguous<T_info>::release_mem() {
	updated_.resize(0);
	store_.resize(0);
	quantities_.resize(0);
}

}	// end namespace GeneralIterationPack

#endif	// ITER_QUANITY_ACCESS_CONTINUOUS_DEF_H
