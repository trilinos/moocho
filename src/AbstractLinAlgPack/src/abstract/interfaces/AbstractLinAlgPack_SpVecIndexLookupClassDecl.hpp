// //////////////////////////////////////////////////////////////////////
// SpVecIndexLookupClassDecl.h
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

#ifndef SPVEC_INDEX_LOOKUP_CLASS_DECL_H
#define SPVEC_INDEX_LOOKUP_CLASS_DECL_H

#include <stdexcept>

#include "AbstractLinAlgPackTypes.h"

namespace AbstractLinAlgPack {

namespace SparseVectorUtilityPack {

// ///////////////////////////////////////////////////////////////////////
///
/** Sparse Vector Index Lookup and Caching class.
  *
  * This class is used to perform a lookup for elements in a sparse vector
  * stored as an array of nonzero elements of a templated type T_Element.
  * The type T_Element must conform to the SparseElementTemplateInterface
  * specification.  These elements must be sorted in accending order.
  *
  * The default C++ assignment operator and copy constructor are allowed.
  */
template <class T_Element>
class SpVecIndexLookup {
public:

	/** @name Public types */
	//@{

	///
	typedef T_Element							element_type;
	///
	typedef typename element_type::index_type	index_type;
	///
	typedef ptrdiff_t							difference_type;
	///
	class NoSpVecSetException : public std::logic_error
	{public: NoSpVecSetException(const std::string& what_arg) : std::logic_error(what_arg) {}};
	///
	class InvalidInternalStateException : public std::logic_error
	{public: InvalidInternalStateException(const std::string& what_arg) : std::logic_error(what_arg) {}};
	///
	enum UpperLower { UPPER_ELE, LOWER_ELE };
	///
	enum ElementRelation { BEFORE_ELE, AFTER_ELE, EQUAL_TO_ELE };
	/// Struct with members: size_type poss; ElementRelation rel;
	struct poss_type {
		poss_type() : poss(0), rel(EQUAL_TO_ELE) {} 
		poss_type(size_type _poss, ElementRelation _rel) : poss(_poss), rel(_rel) {} 
		size_type			poss;
		ElementRelation		rel;
	};

	//@}

	/** @name Constructors */
	//@{

	///
	/** Construct uninitialized with not sparse vector set (#ele() == 0#) */
	SpVecIndexLookup()
		: ele_(0), nz_(0), offset_(0), index_cached_(0)
	{}

	///
	/** Construct initialize with a sparse vector */
	SpVecIndexLookup(element_type* ele, size_type nz, difference_type offset)
		: ele_(ele), nz_(nz), offset_(offset), index_cached_(0)
	{}

	//@}

	/** @name Sparse vector representation setup */
	//@{

	///
	/** Set a new sparse vector.
	  *
	  * This will wipe out any cache that may be stored.
	  */
	void set_sp_vec(element_type* ele, size_type nz, difference_type offset) {
		ele_ = ele;		nz_ = nz;		offset_ = offset;
		sp_vec_was_modified();
	}

	/// Increment nz only
	void incr_nz() {
		nz_++;
	}

	//@}

	/** @name Sparse vector representation access */
	//@{

	///
	element_type*		ele() const {
		return ele_;
	}
	///
	size_type			nz() const {
		return nz_;
	}
	///
	difference_type		offset() const {
		return offset_;
	}

	//@}

	/** @name Element lookup and caching */
	//@{

	///
	/** Lookup an element and cache the result if a binary search was performed.
	  *
	  * This function should only be used if it can be assumed that the elements
	  * are sorted in assending order.
	  *
	  * If #index# is the same as the previously cached lookup then this function
	  * will execute in O(1) time, otherwise a O(log(nz)) binary search will be
	  * performed to find the element and the result of the lookup will be cached.
	  *
	  * To be able to utilize a previously cached search this function must know
	  * if an upper element or lower element is to be found.\\
	  *
	  * Preconditions:\begin{itemize}
	  *	\item #ele() > 0# (throw #NoSpVecSetException#)
	  * \end{itemize}
	  *
	  * Postconditions:\begin{itemize}
	  *	\item [uplow == lower_ele] #index <= ele()[return.poss].index() + offset()# 
	  *	\item [uplow == upper_ele] #ele()[return.poss].index() + offset() <= index# 
	  * \end{itemize}
	  *
	  * @return		#poss_type# object where #return.poss# gives a position in tye
	  *				underlying array and #return.rel# gives where the element with
	  *				#index# is in relation to this possition.
	  *									[ #BEFORE_ELE#		if #ele()[return.poss].index() + offset() < index#	]
	  *					#return.rel# =	[ #EQUAL_TO_ELE#	if #ele()[return.poss].index() + offset() == index#	]
	  *									[ #AFTER_ELE#		if #ele()[return.poss].index() + offset() > index#	]
	  */
	poss_type find_poss(index_type index, UpperLower uplow) const;

	///
	/** Lookup an element.
	  *
	  * Lookup an exact match for an element.  If the element is not found, the
	  * end of the sequence will be returned.
	  * 
	  * If is_sorted == true then a binary search will be used (O(log(nz)).  If is_sorted==false
	  * a sequential search will be used (O(nz)).  No result is cached here.
	  * 
	  * Preconditions:\begin{itemize}
	  *	\item #ele() > 0# (throw #NoSpVecSetException#)
	  * \end{itemize}
	  *
	  * Postconditions:\begin{itemize}
	  *	\item [return == nz()] No element exits with this index 
	  *	\item [return < nz()] #index == ele()[return].index() + offset()# 
	  * \end{itemize}
	  */
	size_type find_element( index_type index, bool is_sorted ) const;
	
	//@}

	/** @name State management */
	//@{

	///
	/** Called by client to inform the object that the sparse vector was modified
	  * so that the cache may be invalidated.
	  */
	void sp_vec_was_modified() {
		index_cached_ = 0;
	}

	///
	/** Called by client to ensure that the internal state is valid.
	  *
	  * If #ele() != 0# but #nz == 0# or #ele()[0].index() + offset() < 1# then this
	  * is an invalid sparse vector and the function will throw a #NoSpVecSetException#
	  * exception.  It is up to the client to ensure that a valid sparse vector is set.
	  *
	  * If there in some other problem with internal state then an exception
	  * #InvalidInternalStateException# will be thrown.  The error message will be
	  * meant for a developer to inspect.
	  */
	void validate_state() const;

	//@}

private:
	// ///////////////////////////////////////////////////////////////////////////////
	// Private types

	// //////////////////////////////////////////////////////////////////////////////
	// Private data members

	element_type*		ele_;		// pointer to array of elements
	size_type			nz_;		// number of nonzero elements in ele_
	difference_type		offset_;	// offset for each element in ele_.  The actuall index for
									// each element i is ele_[i].index() + offset().
	mutable index_type		index_cached_;	// index of last binary search
	mutable size_type		poss_cached_;	// possition last looked up
	mutable ElementRelation	ele_rel_cached_;// Specifies where the element looked up is in relation
											// to the element at poss_cached:
											//		before_ele:	zero element before poss_cached_
											//		after_ele:	zero element after poss_cached_
											//		equal_to_ele: nonzero element at poss_cached_

	// ////////////////////////
	// Private member functions

	/// Assert that a sparse vector has been set up
	void assert_sp_vec_setup() const {
		if(!ele_)
			throw NoSpVecSetException("The sparse vector (ele, nz, offset) has not been set yet");
	}

	/// Adjust the cached possition
	size_type adjust_cached_poss(UpperLower uplow) const;

	///
	/** Perform a binary search for an element in the sparse vector.
	  *
	  * @param	index	[I]	index begin looked for
	  * @param	uplow	[I]	whether it is an upper (#UPPER_ELE#) or lower (#LOWER_ELE#) element needed
	  * @param	poss	[O]	possition where the element was found.  If #uplow == LOWER_ELE# then
	  *						#poss# will be the largest possible integer that satisfies:
	  *						#index <= ele_[poss].index()#.  If #uplow == UPPER_ELE# then #poss#
	  *						will be the smallest possible integer that satisfies:
	  *						#ele_[poss].index() <= index#
	  *	@param	ele_rel	[O]	Where the element with #index# is in relation to the element at #poss#
	  *						returned.  There are three possible values.
	  *						#BEFORE_ELE# :	The element saught is a zero element that is before
	  *										the element at #poss#
	  *						#AFTER_ELE# :	The element saught is a zero element that is after the 
	  *										element at #poss#
	  *						#EQUAL_TO_POSS#: This is a nonzero elment that exists at possition #poss#
	  */
	poss_type binary_ele_search(index_type index, UpperLower uplow) const;


};	// end class SpVecIndexLookup

}	// namespace SparseVectorUtilityPack

} // end namespace AbstractLinAlgPack

#endif // SPVEC_INDEX_LOOKUP_CLASS_DECL_H
