// //////////////////////////////////////////////////////////////////////
// SparseVectorClassDef.h
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

#ifndef SPARSE_VECTOR_CLASS_DEF_H
#define SPARSE_VECTOR_CLASS_DEF_H

#include <algorithm>
#include <sstream>

#include "SparseVectorClassDecl.h"
#include "compare_element_indexes.h"

namespace AbstractLinAlgPack {

// Non-member non-public utility functions

// /////////////////////////////////////////////////////////////////////////////////////
// Definitions of non-member functions

template<class T_Element>
SparseVectorSlice<T_Element>
create_slice(const SparseVectorUtilityPack::SpVecIndexLookup<T_Element>& index_lookup
	, size_type size, Range1D rng)
{
	// Check preconditions
	if(rng.full_range()) {
		rng = Range1D(1,size);
	}
	else {
		if(rng.ubound() > size)
			throw std::out_of_range("create_slice(...): Error, rng.ubound() > size");
	}

	// If there are no elements then any subregion will also not have any elements.
	if(!index_lookup.nz())
		return SparseVectorSlice<T_Element>(0, 0, 0, rng.ubound() - rng.lbound() + 1, true);

	// Create the slice (assumed sorted oviously).
	typedef SparseVectorUtilityPack::SpVecIndexLookup<T_Element> SpVecIndexLookup;
	index_lookup.validate_state();
	SpVecIndexLookup::poss_type
		lower_poss = index_lookup.find_poss(rng.lbound(), SpVecIndexLookup::LOWER_ELE),
		upper_poss = index_lookup.find_poss(rng.ubound(), SpVecIndexLookup::UPPER_ELE);
	if( lower_poss.poss == upper_poss.poss
		&& (	lower_poss.rel == SpVecIndexLookup::AFTER_ELE
			||	upper_poss.rel == SpVecIndexLookup::BEFORE_ELE )
	  )
	{
		// The requested subvector does not contain any elements.
		return SparseVectorSlice<T_Element>(0, 0, 0, rng.ubound() - rng.lbound() + 1, true);
	}
	else {
		// There are nonzero elements
		return SparseVectorSlice<T_Element>(index_lookup.ele() + lower_poss.poss
			, upper_poss.poss - lower_poss.poss + 1, index_lookup.offset() - rng.lbound() + 1
			, rng.ubound() - rng.lbound() + 1, true);
	}
}

// /////////////////////////////////////////////////////////////////////////////////////
// Definitions of members for SparseVector<>

template <class T_Element, class T_Alloc>
SparseVector<T_Element,T_Alloc>::SparseVector(
		const SparseVector<T_Element,T_Alloc>& sp_vec )
	:
	    alloc_(sp_vec.alloc_), size_(sp_vec.size_), max_nz_(sp_vec.max_nz_)
	  , assume_sorted_(sp_vec.assume_sorted_)
	  , know_is_sorted_(sp_vec.know_is_sorted_)
{
	// Allocate the memory for the elements and set the memory of the sparse vector.
	index_lookup_.set_sp_vec(
#ifdef _PG_CXX
		new element_type[max_nz_]
#else
		alloc_.allocate(max_nz_,NULL)
#endif
		,sp_vec.nz(),sp_vec.offset());
	// Perform an uninitialized copy of the elements
	iterator		ele_to_itr		= index_lookup_.ele();
	const_iterator	ele_from_itr	= sp_vec.begin();
	while(ele_from_itr != sp_vec.end()) {
#ifdef _PG_CXX
		new (ele_to_itr++) element_type(*ele_from_itr++);
#else
		alloc_.construct(ele_to_itr++,*ele_from_itr++);
#endif
	}
}

template <class T_Element, class T_Alloc>
SparseVector<T_Element,T_Alloc>::SparseVector(
		  SparseVectorSlice<T_Element> sp_vec_slc
		, const allocator_type& alloc )
	:
		alloc_(alloc), size_(sp_vec_slc.dim()), max_nz_(sp_vec_slc.nz())
	  , assume_sorted_(sp_vec_slc.is_sorted())
	  , know_is_sorted_(false)
{
	// Allocate the memory for the elements and set the memory of the sparse vector.
	index_lookup_.set_sp_vec(
#ifdef _PG_CXX
		new element_type[max_nz_]
#else
		alloc_.allocate(max_nz_,NULL)
#endif
		, sp_vec_slc.nz()
		,sp_vec_slc.offset() );
	// Perform an uninitialized copy of the elements
	iterator		ele_to_itr		= index_lookup_.ele();
	const_iterator	ele_from_itr	= sp_vec_slc.begin();
	while(ele_from_itr != sp_vec_slc.end()) {
#ifdef _PG_CXX
		new (ele_to_itr++) element_type(*ele_from_itr++);
#else
		alloc_.construct(ele_to_itr++,*ele_from_itr++);
#endif
	}
}

template <class T_Element, class T_Alloc>
SparseVector<T_Element,T_Alloc>&
SparseVector<T_Element,T_Alloc>::operator=(
	const SparseVector<T_Element,T_Alloc>& sp_vec)
{
	if(this == &sp_vec) return *this;	// assignment to self

	know_is_sorted_ = sp_vec.know_is_sorted_;
	assume_sorted_ = sp_vec.assume_sorted_;

	if( max_nz() < sp_vec.nz() ) {
		// need to allocate more storage
		resize(0,0);	// free current storage
		resize(sp_vec.dim(),sp_vec.nz(),sp_vec.offset());
	}
	else if( nz() ) {
		// Don't allocate new memory, just call distructors on current elements
		// and reset to uninitialized.
		for(iterator ele_itr = begin(); ele_itr != end();) {
#ifdef _PG_CXX
			(ele_itr++)->~element_type();
#else
			alloc_.destroy(ele_itr++);
#endif
		}
		// Set the other data
		size_ = sp_vec.size_;
	}
	
	// set nz and offset
	index_lookup_.set_sp_vec(index_lookup_.ele(),sp_vec.nz(),sp_vec.offset()); 

	if( sp_vec.nz() ) {
		// Perform an uninitialized copy of the elements
		iterator		ele_to_itr		= index_lookup_.ele();
		const_iterator	ele_from_itr	= sp_vec.begin();
		while(ele_from_itr != sp_vec.end()) {
#ifdef _PG_CXX
			new (ele_to_itr++) element_type(*ele_from_itr++);
#else
			alloc_.construct(ele_to_itr++,*ele_from_itr++);
#endif
		}
	}

	return *this;
}

template <class T_Element, class T_Alloc>
SparseVector<T_Element,T_Alloc>&
SparseVector<T_Element,T_Alloc>::operator=(
	const SparseVectorSlice<T_Element>& sp_vec_slc )
{
	know_is_sorted_ = false;
	assume_sorted_ = sp_vec_slc.is_sorted();

	if(max_nz() < sp_vec_slc.nz()) {
		// need to allocate more storage
		resize(0,0);	// free current storage
		resize(sp_vec_slc.dim(),sp_vec_slc.nz(),sp_vec_slc.offset());
	}
	else if( nz() ) {
		// Don't allocate new memory, just call distructors on current elements
		// and reset to uninitialized.
		for(iterator ele_itr = begin(); ele_itr != end();) {
#ifdef _PG_CXX
			(ele_itr++)->~element_type();
#else
			alloc_.destroy(ele_itr++);
#endif
		}
		// Set the other data
		size_ = sp_vec_slc.dim();
	}
	
	// set nz and offset
	index_lookup_.set_sp_vec(index_lookup_.ele(),sp_vec_slc.nz()
		,sp_vec_slc.offset()); 

	if( sp_vec_slc.nz() ) {
		// Perform an uninitialized copy of the elements
		iterator		ele_to_itr		= index_lookup_.ele();
		const_iterator	ele_from_itr	= sp_vec_slc.begin();
		while(ele_from_itr != sp_vec_slc.end()) {
#ifdef _PG_CXX
			new (ele_to_itr++) element_type(*ele_from_itr++);
#else
			alloc_.construct(ele_to_itr++,*ele_from_itr++);
#endif
		}
	}

	return *this;
}

template <class T_Element, class T_Alloc>
EOverLap SparseVector<T_Element,T_Alloc>::overlap(const SparseVectorSlice<T_Element>& sv) const
{
	if(!sv.dim()) return AbstractLinAlgPack::NO_OVERLAP;

	const_iterator										this_begin	= begin();
	SparseVectorSlice<T_Element>::const_iterator		sv_begin	= sv.begin();

	if( this_begin == sv_begin && end() == sv.end() )
	{
		return AbstractLinAlgPack::SAME_MEM;
	}

	if(		( this_begin < sv_begin && end() < sv_begin )
		||	( sv_begin < this_begin && sv.end() < this_begin )	)
	{
		return AbstractLinAlgPack::NO_OVERLAP;
	}

	return AbstractLinAlgPack::SOME_OVERLAP;
}

template <class T_Element, class T_Alloc>
void SparseVector<T_Element,T_Alloc>::resize(size_type size, size_type max_nz
	, difference_type offset)
{
	// free existing storage
	if(index_lookup_.ele()) {
		for(element_type* p = index_lookup_.ele(); p < index_lookup_.ele() + index_lookup_.nz(); ++p) {
#ifdef _PG_CXX
			p->~element_type();
#else
			alloc_.destroy(p);
#endif
		}
#ifdef _PG_CXX
		delete [] index_lookup_.ele();
#else
		alloc_.deallocate(index_lookup_.ele(), max_nz_);
#endif
	}
	
	// reinitialize
	max_nz_ = 0;
	know_is_sorted_ = false;

	if(max_nz) {
		// reallocate storage
#ifdef _PG_CXX
		index_lookup_.set_sp_vec( new element_type[max_nz_ = max_nz], 0, offset );
#else
		index_lookup_.set_sp_vec( alloc_.allocate( max_nz_ = max_nz, 0 ), 0, offset );
#endif
		size_ = size;
	}
	else {
		// reinitialize to no storage
		index_lookup_.set_sp_vec( 0, 0, offset );
		size_ = size;	// allow size to be nonzero with nz = 0
	}
}

template <class T_Element, class T_Alloc>
void SparseVector<T_Element,T_Alloc>::uninitialized_resize(size_type size, size_type nz, size_type max_nz
	, difference_type offset)
{
	if(nz > max_nz)
		throw std::length_error(
			"SparseVector<...>::uninitialized_resize(...) : nz can not be greater"
			" than max_nz" );
	resize(size,max_nz,offset);
	index_lookup_.set_sp_vec(index_lookup_.ele(), nz, index_lookup_.offset());
}

template <class T_Element, class T_Alloc>
void SparseVector<T_Element,T_Alloc>::insert_element(element_type ele)
{
	assert_space(1);
	assert_is_sorted();
	// Find the insertion point
	if( nz() ) {
		typedef SparseVectorUtilityPack::SpVecIndexLookup<T_Element> SpVecIndexLookup;
		index_lookup_.validate_state();
		SpVecIndexLookup::poss_type
			poss = ( nz() ? index_lookup_.find_poss(ele.index(), SpVecIndexLookup::LOWER_ELE)
					 : SpVecIndexLookup::poss_type(0,SpVecIndexLookup::BEFORE_ELE) );
		// Make sure this element does not already exist!
		if( nz() && poss.rel == SpVecIndexLookup::EQUAL_TO_ELE )
			throw std::length_error(
				"SparseVector<...>::insert_element(...) : Error, this index"
				" all ready exists!" );
		const size_type
			insert_poss = (poss.rel == SpVecIndexLookup::BEFORE_ELE ? poss.poss : poss.poss+1);
		// Copy elements out of the way to make room for inserted element
		std::copy_backward( // This assumes element_type supports assignment!
			index_lookup_.ele() + insert_poss, index_lookup_.ele() + index_lookup_.nz()
			, index_lookup_.ele() + index_lookup_.nz() + 1 );
		index_lookup_.ele()[insert_poss] = ele;
		index_lookup_.incr_nz();
	}
	else { // The first element we are adding!
		index_lookup_.ele()[0] = ele;
		index_lookup_.incr_nz();
	}
}

template <class T_Element, class T_Alloc>
void SparseVector<T_Element,T_Alloc>::sort() {
	if( index_lookup_.nz() > 0 )
		std::stable_sort(begin(),end(),compare_element_indexes_less<element_type>());
	know_is_sorted_ = true;
}

template <class T_Element, class T_Alloc>
void SparseVector<T_Element,T_Alloc>::assert_valid_and_sorted() const
{

	if(!index_lookup_.nz()) return;	// An empty sparse vector is certainly valid

	// Loop through the elements.  If they are sorted then duplicate
	// elements will be adjacent to each other so they will be easy to 
	// find.
	typename T_Element::index_type last_index;
	for(T_Element* p = index_lookup_.ele();
		p < index_lookup_.ele() + index_lookup_.nz(); ++p)
	{
		typename T_Element::index_type curr_index = p->index() + offset();
		if((1 > curr_index) || (curr_index > dim())) {
			// element index + offset is not in range
			std::ostringstream omsg;
			omsg	<< "SparseVector<...>::assert_valid_and_sorted():"
					<< " Error, not in range:  element (0-based) " << p - index_lookup_.ele() - 1
					<< " with index + offset = " << curr_index
					<< " is not in range";
			throw	std::out_of_range(omsg.str());
		}
		if(p == index_lookup_.ele()) { // skip these tests for the first element
			last_index = curr_index;
			continue;
		}
		if(curr_index < last_index) {
			// the elements are not sorted in accending order
			std::ostringstream omsg;
			omsg	<< "SparseVector<...>::assert_valid_and_sorted():"
					<< " Error, not sorted:  element (0-based) " << p - index_lookup_.ele() - 1
					<< " and " << p - index_lookup_.ele() << " are not in assending order";
			throw	NotSortedException(omsg.str());
		}
		if(curr_index == last_index) {
			// These are duplicate indexes
			std::ostringstream omsg;
			omsg	<< "SparseVector<...>::assert_valid_and_sorted():"
					<< " Error, duplicate indexes:  element (0-based) " << p - index_lookup_.ele() - 1
					<< " and " << p - index_lookup_.ele() << " have duplicate indexes";
			throw	DuplicateIndexesException(omsg.str());
		}
		last_index = curr_index;
	}
}


// /////////////////////////////////////////////////////////////////////////////////////
// Definitions of members for SparseVectorSlice<>

template <class T_Element>
EOverLap SparseVectorSlice<T_Element>::overlap(const SparseVectorSlice<T_Element>& sv) const
{
	if(!sv.dim()) return AbstractLinAlgPack::NO_OVERLAP;

	const_iterator					this_begin	= begin(),
									sv_begin	= sv.begin();

	if( this_begin == sv_begin && end() == sv.end() )
	{
		return AbstractLinAlgPack::SAME_MEM;
	}

	if(		( this_begin < sv_begin && end() < sv_begin )
		||	( sv_begin < this_begin && sv.end() < this_begin )	)
	{
		return AbstractLinAlgPack::NO_OVERLAP;
	}

	return AbstractLinAlgPack::SOME_OVERLAP;
}

} // end namespace AbstractLinAlgPack 

#endif // SPARSE_VECTOR_CLASS_DEF_H