// //////////////////////////////////////////////////////////////////////
// TransSparseCOOElementViewIter.hpp
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

#ifndef TRANS_SPARSE_COO_ELEMENT_VIEW_ITER_H
#define TRANS_SPARSE_COO_ELEMENT_VIEW_ITER_H

#include <iterator>

#include "SparseLinAlgPackTypes.hpp"

namespace SparseLinAlgPack {

///
/** Templateded iterator for iterating through a set of COO Matix 
  * elements but viewing them in transpose.
  *
  * The iterator type T_Iter must yield a type that conforms to
  * the SparseCOOElementTemplateInterface specification.  In particular
  * it must support itr->value(), itr->row_i(), and itr->col_i().
  *
  * Any category of iterator can be used with this class.  The relavant iterator traits
  * have to be added to the template list since MS VC++ 5.0 does not support
  * partial template specialization and therefore std::iterator_traits<> does
  * not exist.
  *
  * There is a lot of code here but for an iterator object, only the underlying
  * iterator is stored.  Thus, because of inlining accessing a COO matrix through
  * this transposed iterator is just as efficient as assessing it through
  * its non-transposed iterator.  This is important for efficiency reasons.
  *
  * The default assignment operator and copy constructor
  * are allowed.
  *
  * The template arguments are:\begin{itemize}
  *	\item	T_Iter		The type of the iterator that is used to yield COO sparse elements
  *	\item	T_IterCat	The category of iterator (forward, random access etc.)
  *	\item	T_Indice	The type returned by itr->indice()
  *	\item	T_ValRef	The type for the value returned by value().  For example, if T_Ele
  *						is a const type then T_ValRef should be a constant reference and
  *						if it is a non-const type it should be a non-const reference.
  */
template <class T_Iter, class T_IterCat, class T_Indice, class T_ValRef, class T_Diff>
class TransSparseCOOElementViewIter {
public:
	/** @name Public types. */
	//@{

	///
	/** Type for the object that is returned for the transpose sparse
	  * element.
	  *
	  * The default copy constructor is allowed but not the default constructor
	  * or assignment operator.
	  *
	  * Here the type of element may be const or nonconst in accordance with the iterator's
	  * type.  The element view (const or nonconst) is just a view to this type.
	  */
	template <class TT_Iter, class TT_IterCat, class TT_Indice, class TT_ValRef, class TT_Diff>
	class ElementView {
	public:

		// friends

		friend
		class TransSparseCOOElementViewIter<TT_Iter,TT_IterCat,TT_Indice,TT_ValRef,TT_Diff>;

		// typedefs

		///
		typedef TT_Indice						indice_type;
		///
		typedef TT_ValRef						value_ref_type;

		///
		/** Construct with an iterator to the first element.
		  */
		ElementView(const TT_Iter& iter) : encap_iter_(iter)
		{}

		// access functions

		///
		value_ref_type value() const {
			return encap_iter_->value();
		}

		/// returns col_j() of the underlying COO element
		indice_type row_i() const {
			return encap_iter_->col_j();
		}

		/// returns row_i() of the underlying COO element
		indice_type col_j() const {
			return encap_iter_->row_i();
		}

	private:
		TT_Iter	encap_iter_;// iterator that is manipulated by
							// the host iterator

		// not defined and not to be called
		ElementView();
		ElementView& operator=(const ElementView&);
	};

	// Local typedefs

	///
	typedef ElementView<T_Iter,T_IterCat,T_Indice
							,T_ValRef,T_Diff>				element_view_type;
	///
	typedef T_Iter											encap_iter_type;
	///
	typedef TransSparseCOOElementViewIter<T_Iter,T_IterCat
								,T_Indice,T_ValRef,T_Diff>	iterator_type;

	// Standard C+ library types for iterators

	///
	typedef	T_IterCat										iterator_category;
	///
	typedef	element_view_type								value_type;
	///
	typedef element_view_type&								reference_type;
	///
	typedef element_view_type*								pointer_type;
	///
	typedef	T_Diff											difference_type;
	///
	typedef	size_t											distance_type;

	//	end Public Types
	//@}

	/** @name Constructors. */
	//@{

	///
	/** Construct with the iterator of COO elements to transpose.
	  *
	  * Note that this constructor also allows an implicit conversion
	  * from the nontransposed iterator to the transposed iterator.
	  * This may not be desireable.
	  */
	TransSparseCOOElementViewIter(T_Iter itr) : element_view_(itr)
	{}

	//@}

	/** @name Iterator control functions. */
	//@{

	// Access
	reference_type		operator*()	const {
		return const_cast<element_view_type&>(element_view_);
	}
	pointer_type		operator->() const {
		return const_cast<element_view_type*>(&element_view_);
	}
	value_type			operator[](distance_type n)	const {
		return element_view_type(element_view_.encap_iter_ + n);
	}
	// Incrementation
	// ++itr
	iterator_type&		operator++() {
		element_view_.encap_iter_++;
		return *this;
	}
	// itr++
	const iterator_type	operator++(int) {
		iterator_type tmp = *this;
		++*this;
		return tmp;
	}
	// --itr
	iterator_type&		operator--() {
		element_view_.encap_iter_--;
		return *this;
	}
	// itr--
	const iterator_type	operator--(int) {
		iterator_type tmp = *this;
		--*this;
		return tmp;
	}
	// itr + n 
	iterator_type			operator+(distance_type n) {
		return iterator_type(element_view_.encap_iter_ + n);
	}
	const iterator_type		operator+(distance_type n) const {
		return iterator_type(element_view_.encap_iter_ + n);
	}
	// itr += n
	iterator_type&		operator+=(distance_type n) {
		element_view_.encap_iter_ += n;
		return *this;
	}
	// itr - n
	iterator_type			operator-(distance_type n) {
		return iterator_type(element_view_.encap_iter_ - n);
	}
	const iterator_type	operator-(distance_type n) const {
		return iterator_type(element_view_.encap_iter_ - n);
	}
	// itr -= n
	iterator_type&		operator-=(distance_type n) {
		element_view_.encap_iter_ -= n;
		return *this;
	}
	// distance = itr1 - itr2 (distance between elements)
	distance_type		operator-(const iterator_type& itr) const	{
		return element_view_.encap_iter_ - itr.element_view_.encap_iter_;
	}

	//@}

	/** @name Comparison Operators.
	 */
	//@{

	// <
	bool operator<(const iterator_type& itr)
	{	
		return element_view_.encap_iter_ < itr.element_view_.encap_iter_;
	}
	// <=
	bool operator<=(const iterator_type& itr)
	{	
		return element_view_.encap_iter_ <= itr.element_view_.encap_iter_;
	}
	// >
	bool operator>(const iterator_type& itr)
	{	
		return element_view_.encap_iter_ > itr.element_view_.encap_iter_;
	}
	// >=
	bool operator>=(const iterator_type& itr)
	{	
		return element_view_.encap_iter_ >= itr.element_view_.encap_iter_;
	}
	// ==
	bool operator==(const iterator_type& itr)
	{	
		return element_view_.encap_iter_ == itr.element_view_.encap_iter_;
	}
	// !=
	bool operator!=(const iterator_type& itr)
	{	
		return element_view_.encap_iter_ != itr.element_view_.encap_iter_;
	}

	//@}

private:
	element_view_type	element_view_;
	
	// not defined and not to be called
	TransSparseCOOElementViewIter();

};	// end class TransSparseCOOElementViewIter

// ///////////////////////////////////////////////////////////////////////////////////////
// Nonmember functions

// Allow distance_type as lhs argument in n + itr

//template <class Iter, class Cat, class Indice, class ValRef, class Diff>
//inline TransSparseCOOElementViewIter<Iter,Cat,Indice,ValRef,Diff>
//operator+(Diff n, TransSparseCOOElementViewIter<Iter,Cat,Indice,ValRef,Diff> itr)
//{
//	return itr + n;
//}

template <class Iter, class Cat, class Indice, class ValRef, class Diff>
inline TransSparseCOOElementViewIter<Iter,Cat,Indice,ValRef,Diff>
operator+(Diff n, const TransSparseCOOElementViewIter<Iter,Cat,Indice,ValRef,Diff> itr)
{
	return itr + n;
}

} // end namespace SparseLinAlgPack 

#endif // TRANS_SPARSE_COO_ELEMENT_VIEW_ITER_H
