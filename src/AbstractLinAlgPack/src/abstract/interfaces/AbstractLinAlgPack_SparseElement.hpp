// //////////////////////////////////////////////////////////////////////
// SparseElement.h
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

#ifndef SPARSE_ELEMENT_H
#define SPARSE_ELEMENT_H

#include "AbstractLinAlgPackTypes.h"

namespace AbstractLinAlgPack {

///
/** Sparse storage element type.
  *
  * This class abstracts a sparse element of a templated
  * type.  It is ment to be used in a sparse vector.  Objects of
  * this type are designed so that the size of the object is
  * the same at least two value_type objects.
  *
  * The default assignment operator and copy constructor
  * are allowed.
  */
template <class T_Index, class T_Value>
class SparseElement {
public:
	/** @name Public Typedefs. */
	//@{

	///
	typedef T_Value						value_type;
	///
	typedef T_Index						index_type;

	//@}

	/** @name Constructors */
	//@{

	/// Construct uninitialized (#value() == 0.0#, #index() == 0#).
	SparseElement()
	{
		idx_pad_.index_  = 0;
		value_           = 0.0;
	}

	/// Construct with a value and index set
	SparseElement(index_type index, value_type value)
	{
		idx_pad_.index_ = index;
		value_          = value;
	}
	
	//@}

	/** @name Value and index access */
	//@{ 

	///
	value_type& value()
	{
		return value_;
	}
	///
	const value_type& value() const
	{
		return value_;
	}
	///
	const index_type& index() const
	{
		return idx_pad_.index_;
	}
	/// Initialize
	void initialize(index_type index, value_type value) {
		idx_pad_.index_ = index;
		value_          = value;
	}	
	/// Change the index
	void change_index(index_type index)
	{
		idx_pad_.index_ = index;
	}

	//@}

private:
	union index_and_padding {
		value_type  dummy;     // This is just included for alignment
		index_type index_;   // so that sizeof(this) == 2*sizeof(value_type)
	};
	index_and_padding		idx_pad_;
	value_type				value_;

};	// end class SparseElement

} // end namespace AbstractLinAlgPack 

#endif // SPARSE_ELEMENT_H
