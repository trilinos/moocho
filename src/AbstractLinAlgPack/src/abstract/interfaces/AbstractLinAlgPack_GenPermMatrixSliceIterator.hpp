// /////////////////////////////////////////////////////////
// GenPermMatrixSliceIterator.h
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

#ifndef GEN_PERM_MATRIX_SLICE_ITERATOR_H
#define GEN_PERM_MATRIX_SLICE_ITERATOR_H

#include <assert.h>

#include <iterator>

#include "AbstractLinAlgPackTypes.h"

namespace AbstractLinAlgPack {

namespace GenPermMatrixSliceIteratorPack {

///
enum EOrderedBy { BY_ROW, BY_COL, BY_ROW_AND_COL, UNORDERED };

///
/** External storage of a row and column indice.
  * This is required for creating a temporary in an assignment operation
  * in a sorting algorithm (like std::sort(...)).
  */
template< class T >
class external_row_col_value_type {
public:
	///
	typedef T			index_type;
	///
	typedef ptrdiff_t	difference_type;
	///
	external_row_col_value_type(
		  difference_type	row_off
		, difference_type	col_off
		, index_type		row_i
		, index_type		col_j
		)
	:
		row_off_(row_off), col_off_(col_off), row_i_(row_i), col_j_(col_j)
	{}
	difference_type	row_off_;
	difference_type	col_off_;
	index_type		row_i_;
	index_type		col_j_;
};

///
/** Internal storage for the iterator of the
  * row and column indices.
  */
template< class T >
class row_col_value_type {
public:
	///
	typedef T			index_type;
	///
	typedef ptrdiff_t	difference_type;
	///
	row_col_value_type( 
		  difference_type	row_off
		, difference_type	col_off
		, index_type		row_i[]
		, index_type		col_j[]
		, size_type			nz
		);
	///
	void bind_view( const row_col_value_type<T>& val );
	///
	void increment(difference_type);
	///
	index_type	row_i() const;
	///
	index_type	col_j() const;
	/// May be NULL
	index_type* row_i_ptr() const;
	///
	row_col_value_type<T>& operator=( const row_col_value_type<T>& val );
	///
	operator const external_row_col_value_type<T>() const
	{
		return external_row_col_value_type<T>(row_off_,col_off_,*row_i_,*col_j_);
	}
	///
	row_col_value_type<T>& operator=( const external_row_col_value_type<T>& val )
	{
		assert( row_off_ == val.row_off_ );
		assert( col_off_ == val.col_off_ );
		*row_i_ = val.row_i_;
		*col_j_ = val.col_j_;
		return *this;
	}
	
private:
	difference_type	row_off_;
	difference_type	col_off_;
	index_type		*row_i_;
	index_type		*col_j_;
	size_type		nz_;
	int				k_;	// zero based
	///
	void assert_in_range() const;
	/// Not defined and not to be called
	row_col_value_type();

};	// end class row_col_value_type

/// Swap row_col_value_type<T> objects
template< class T >
inline
void swap( row_col_value_type<T>& v1, row_col_value_type<T>& v2 )
{
	row_col_value_type<T> tmp = v1;
	v1 = v2;
	v2 = tmp;
}

///
/** This is a full random access iterator for accessing row and colunmn
  * indices.
  */
template< class T >
class row_col_iterator
#if defined(_WINDOWS) || defined(_INTEL_CXX) || defined(_PG_CXX) 
	: public std::iterator< std::random_access_iterator_tag, external_row_col_value_type<T>, ptrdiff_t >
#endif
{
public:
	///
	typedef T								index_type;
	///
	typedef	std::random_access_iterator_tag	iterator_category;
	///
	typedef	external_row_col_value_type<T>	value_type;
	///
	typedef row_col_value_type<T>&			reference;
	///
	typedef row_col_value_type<T>*			pointer;
	///
	typedef	ptrdiff_t						difference_type;
	/// Null pointer!
	row_col_iterator();
	///
	row_col_iterator(
		 difference_type	row_off
		,difference_type	col_off
		,index_type		row_i[]
		,index_type		col_j[]
		,size_type			nz			// Number of elements in row_i[] and col_j[]
		);
	///
	row_col_iterator<T>& operator=( const row_col_iterator<T>& itr );
	///
	reference operator*();
	///
	reference operator*() const;
	///
	pointer operator->() const;
	/// itr + a
	row_col_iterator<T>	operator+(difference_type);
	/// itr - a
	row_col_iterator<T>	operator-(difference_type);
	/// itr += a
	row_col_iterator<T>&	operator+=(difference_type);
	/// itr -= a
	row_col_iterator<T>&	operator-=(difference_type);
	/// ++itr
	row_col_iterator<T>&	operator++();
	/// itr++
	const row_col_iterator<T> operator++(int);
	/// --itr
	row_col_iterator<T>&	operator--();
	/// itr--
	const row_col_iterator<T> operator--(int);
	/// Difference
	difference_type operator-(const row_col_iterator<T>& itr) const;
	/// itr1 < itr2
	bool operator<( const row_col_iterator<T>& itr);
	/// itr1 <= itr2
	bool operator<=( const row_col_iterator<T>& itr);
	/// itr1 > itr 2
	bool operator>( const row_col_iterator<T>& itr);
	/// itr1 >= itr2
	bool operator>=( const row_col_iterator<T>& itr);
	/// itr1 == itr2
	bool operator==( const row_col_iterator<T>& itr);
	/// itr1 != itr2
	bool operator!=( const row_col_iterator<T>& itr);
	/// !itr (check for null)
	bool operator!() const;
	
private:
	mutable row_col_value_type<T>	value_;
	
};	// end class row_col_iterator<T>

// //////////////////////////////////////////////////////////
// Inline members for row_col_value_type<T>

template<class T>
inline
row_col_value_type<T>::row_col_value_type( 
		  difference_type	row_off
		, difference_type	col_off
		, index_type		row_i[]
		, index_type		col_j[]
		, size_type			nz
		)
	:
		row_off_(row_off)
		,col_off_(col_off)
		,row_i_(row_i)
		,col_j_(col_j)
		,nz_(nz)
		,k_(0)
{}

template<class T>
inline
void row_col_value_type<T>::bind_view( const row_col_value_type<T>& val )
{
		row_off_	= val.row_off_;
		col_off_	= val.col_off_;
		row_i_		= val.row_i_;
		col_j_		= val.col_j_;
		nz_			= val.nz_;
		k_			= val.k_;
}

template< class T >
inline
void row_col_value_type<T>::increment(difference_type d)
{
	row_i_	+= d;
	col_j_	+= d;
	k_ 		+= d;
}

template< class T >
inline
row_col_value_type<T>::index_type row_col_value_type<T>::row_i() const
{
	assert_in_range();
	return *row_i_ + row_off_;
}

template< class T >
inline
row_col_value_type<T>::index_type row_col_value_type<T>::col_j() const
{
	assert_in_range();
	return *col_j_ + col_off_;
}

template< class T >
inline
row_col_value_type<T>::index_type* row_col_value_type<T>::row_i_ptr() const
{
	return row_i_;
}

template< class T >
inline
row_col_value_type<T>& row_col_value_type<T>::operator=(
	const row_col_value_type<T>& val )
{
	*row_i_ = *val.row_i_;
	*col_j_ = *val.col_j_;
	return *this;
}

template< class T >
inline
void row_col_value_type<T>::assert_in_range() const
{
	// ToDo: Finish this!
	assert( 0 <= k_ && k_ < nz_ );
}

/// Assert not null
void GPMS_row_col_iterator_assert_not_null(const void* p);

// //////////////////////////////////////////////////////////
// Inline members for row_col_iterator<T>

template< class T >
inline
row_col_iterator<T>::row_col_iterator()
	:
		value_(0,0,NULL,NULL,0)
{}

template< class T >
inline
row_col_iterator<T>::row_col_iterator(
		 difference_type	row_off
		,difference_type	col_off
		,index_type         row_i[]
		,index_type         col_j[]
		,size_type			nz
		)
	:
		value_(row_off,col_off,row_i,col_j,nz)
{}

template< class T >
inline
row_col_iterator<T>& row_col_iterator<T>::operator=( const row_col_iterator<T>& itr )
{
	value_.bind_view( itr.value_ );
	return *this;
}

template< class T >
inline
row_col_iterator<T>::reference
row_col_iterator<T>::operator*()
{
	GPMS_row_col_iterator_assert_not_null(value_.row_i_ptr());
	return value_;
}


template< class T >
inline
row_col_iterator<T>::reference
row_col_iterator<T>::operator*() const
{
	GPMS_row_col_iterator_assert_not_null(value_.row_i_ptr());
	return value_;
}

template< class T >
inline
row_col_iterator<T>::pointer
row_col_iterator<T>::operator->() const
{
	GPMS_row_col_iterator_assert_not_null(value_.row_i_ptr());
	return &value_;
}

template< class T >
inline
row_col_iterator<T>
row_col_iterator<T>::operator+(difference_type d)
{
	row_col_iterator<T> itr = *this;
	itr.value_.increment(d);
	return itr;
}

template< class T >
inline
row_col_iterator<T>
row_col_iterator<T>::operator-(difference_type d)
{
	row_col_iterator<T> itr = *this;
	itr.value_.increment(-d);
	return itr;
}

template< class T >
inline
row_col_iterator<T>&
row_col_iterator<T>::operator+=(difference_type d)
{
	value_.increment(d);
	return *this;
}

template< class T >
inline
row_col_iterator<T>&
row_col_iterator<T>::operator-=(difference_type d)
{
	value_.increment(-d);
	return *this;
}

template< class T >
inline
row_col_iterator<T>&
row_col_iterator<T>::operator++()
{
	value_.increment(1);
	return *this;
}

template< class T >
inline
const row_col_iterator<T>
row_col_iterator<T>::operator++(int)
{
	row_col_iterator<T> itr = *this;
	value_.increment(1);
	return itr;
}

template< class T >
inline
row_col_iterator<T>&
row_col_iterator<T>::operator--()
{
	value_.increment(-1);
	return *this;
}

template< class T >
inline
const row_col_iterator<T>
row_col_iterator<T>::operator--(int)
{
	row_col_iterator<T> itr = *this;
	value_.increment(-1);
	return itr;
}

template< class T >
inline
row_col_iterator<T>::difference_type
row_col_iterator<T>::operator-(const row_col_iterator<T>& itr) const
{
	return value_.row_i_ptr() - itr.value_.row_i_ptr();
}

template< class T >
inline
bool row_col_iterator<T>::operator<( const row_col_iterator<T>& itr)
{
	return ( value_.row_i_ptr() && itr.value_.row_i_ptr() )
			&&  ( value_.row_i_ptr() < itr.value_.row_i_ptr() );
}

template< class T >
inline
bool row_col_iterator<T>::operator<=( const row_col_iterator<T>& itr)
{
	return ( value_.row_i_ptr() && itr.value_.row_i_ptr() )
			&&  ( value_.row_i_ptr() <= itr.value_.row_i_ptr() );
}

template< class T >
inline
bool row_col_iterator<T>::operator>( const row_col_iterator<T>& itr)
{
	return ( value_.row_i_ptr() && itr.value_.row_i_ptr() )
			&&  ( value_.row_i_ptr() > itr.value_.row_i_ptr() );
}

template< class T >
inline
bool row_col_iterator<T>::operator>=( const row_col_iterator<T>& itr)
{
	return ( value_.row_i_ptr() && itr.value_.row_i_ptr() )
			&&  ( value_.row_i_ptr() >= itr.value_.row_i_ptr() );
}

template< class T >
inline
bool row_col_iterator<T>::operator==( const row_col_iterator<T>& itr)
{
	return ( value_.row_i_ptr() && itr.value_.row_i_ptr() )
			&&  ( value_.row_i_ptr() == itr.value_.row_i_ptr() );
}

template< class T >
inline
bool row_col_iterator<T>::operator!=( const row_col_iterator<T>& itr)
{
	return ( value_.row_i_ptr() && itr.value_.row_i_ptr() )
			&&  ( value_.row_i_ptr() != itr.value_.row_i_ptr() );
}

template< class T >
inline
bool row_col_iterator<T>::operator!() const
{
	return  value_.row_i_ptr() == NULL;
}

}	// end namespace GenPermMatrixSliceIteratorPack

}	// end namespace AbstractLinAlgPack

#endif   // GEN_PERM_MATRIX_SLICE_ITERATOR_H
