// //////////////////////////////////////////////////////////////////////////////////////
// COOMatrixClass.h
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

#ifndef COO_MATRIX_CLASS_H
#define COO_MATRIX_CLASS_H

#include <valarray>
#include <vector>

#include "SparseLinAlgPackTypes.h"
#include "MiRefCount.h"

namespace SparseLinAlgPack {

///
/** Sparse Coordinate Matrix abstraction storage class.
  *
  * This class abstracts a fortran style sparse coordinate matrix which is stored
  * is three vectors: (val, ivect, jvect).  This class allows direct access to
  * these arrays for integration with fortran function calls.
  *
  * The row and column indices represented by this class must be 1-based.
  */
class COOMatrix {
public:
	// ///////////////////////////////////////////////////////////////
	// Friends
	
	// ///////////////////////////////////////////////////////////////
	// Public interface

	/** @name {\bf Public Typedefs} */
	//@{
	///
	typedef	SparseLinAlgPack::size_type					size_type;
	///
	typedef	SparseLinAlgPack::indice_type				indice_type;
	///
	typedef SparseLinAlgPack::value_type				value_type;

	//@}

	/** @name Constructors.
	  *
	  * The default copy constructor is used since it has the correct
	  * implementation (memberwise assignment and not binary copy).
	  * In the default copy  constructor, initially these objects will
	  * share memory for #ivect#, #jvect#.  Initially only #val# is
	  * allocated and made unique.  If the sparsity information of the
	  * matrix does not change then the value of the nonzero elements
	  * can be changed without further allocations.  If however, the
	  * row and/or column access in changed, or ivect or jvect is changed
	  * then new allocations will be performed.
	  */
	//@{
	
	/// Consturct with no storage allocated
	COOMatrix();

	//@}

	///
	/** Assignment operator.
	  *
	  * This function has similar behavior w.r.t sharing that the copy constructor
	  * has accept any current storage will be deallocated.
	  */
	COOMatrix& operator=(const COOMatrix& coom);

	///
	/** Resize for a #rows# by #cols# sparse matrix with #nz# elements.
	  *
	  * Any sharing if row or column indices is lost.
	  */
	void resize(size_type rows, size_type cols, size_type nz);

	/// Return the number of rows in the row access view.
	size_type rows() const;
	/// Returns the number of columns in the column access view.
	size_type cols() const;
	//// Return the number of non-zero elements in the underlying COO matrix
	size_type nz() const;

	/** @name COO matrix representation access.
	  *
	  * These member functions give access the the #val#, #ivect#
	  * and #jvect# arrays used to store the underlying nonzero
	  * elements of the full matrix.
	  * 
	  * Since sharing can be used between COOMatrix objects
	  * the client should call the 'const_' functions if
	  * they are only going to be reading these data since if
	  * the others are called on a nonconst object then if
	  * #ivect# or #jvect# is being shared then a freash
	  * copy would be made unnecesarily.
	  */

	/// Return pointer to raw storage array (length #nz()#) for the values of the non-zero elements
	value_type*							val();
	///
	const value_type*					val() const;
	///
	const value_type*					const_val() const;
	///
	/// Return pointer to raw storage array (length #nz()#) for the row indices of the non-zero elements
	indice_type*						ivect();
	///
	const indice_type*					ivect() const;
	///
	const indice_type*					const_ivect() const;
	/// Return pointer to raw storage array (length #nz()#) for the column indices of the non-zero elements
	indice_type*						jvect();
	///
	const indice_type*					jvect() const;
	///
	const indice_type*					const_jvect() const;

	///
	/** Initialize from an input stream.
	  *
	  * The format for the imput is:
	  *
	  * #m  n  nz#\\
	  * #a1:i1:j1  a2:i2:j2 .... anz:inz:jnz#\\
	  *
	  * In the above format, each non-zero element is given as a three item pair:
	  * value of the non-zero element, row indice (1-based) of the non-zero element,
	  * and the column indice (1-based) of the non-zero element.  There must
	  * be no spaces between the numbers and the \':\' charachter and there must
	  * be at least one whitespace character between elements.
	  */
	void initialize(std::istream& istrm);

private:
	// ///////////////////////////////////////////////////////////////////////////
	// Private types

	typedef MemMngPack::RefCount<
		std::valarray<indice_type> >					va_indice_ref_type;

	typedef std::valarray<value_type>					va_value_type;

	typedef std::valarray<indice_type>					va_indice_type;

	// ///////////////////////////////////////////////////////////////////////////
	// Private data members

	size_type					rows_,	// The number of rows this sparse matrix represents
								cols_,	// The numner of columns this sparse matrix represents
								nz_;		// The number of non-zero elements this matrix holds.
	va_value_type				val_;		// The vector of non-zero elements
	va_indice_ref_type			ivect_ref_,	// The vector of row indices
								jvect_ref_;	// The vector of column indices

	// ///////////////////////////////////////////////////////////////////////////
	// Private member functions

};	// end class COOMatrix

// ///////////////////////////////////////////////////////////////////////////////////
// Inline member function definitions for COOMatrix

// constructors
inline COOMatrix::COOMatrix() : rows_(0), cols_(0), nz_(0)
{}
// dimensions
inline COOMatrix::size_type COOMatrix::rows() const {
	return rows_;
}
inline COOMatrix::size_type COOMatrix::cols() const {
	return cols_;
}
inline COOMatrix::size_type COOMatrix::nz() const {
	return nz_;
}
// representation vectors (val, ivect, jvect)
inline COOMatrix::value_type* COOMatrix::val() {
	return &val_[0];
}
inline const COOMatrix::value_type* COOMatrix::val() const {
	return const_val();
}
inline const COOMatrix::value_type* COOMatrix::const_val() const {
	return &const_cast<va_value_type&>(val_)[0];
}
inline COOMatrix::indice_type* COOMatrix::ivect() {
	return &(ivect_ref_.obj())[0];
}
inline const COOMatrix::indice_type* COOMatrix::ivect() const {
	return const_ivect();
}
inline const COOMatrix::indice_type* COOMatrix::const_ivect() const {
	return &const_cast<va_indice_type&>(ivect_ref_.const_obj())[0];
}
inline COOMatrix::indice_type* COOMatrix::jvect() {
	return &jvect_ref_.obj()[0];
}
inline const COOMatrix::indice_type* COOMatrix::jvect() const {
	return const_jvect();
}
inline const COOMatrix::indice_type* COOMatrix::const_jvect() const {
	return &const_cast<va_indice_type&>(jvect_ref_.const_obj())[0];
}

} // end namespace SparseLinAlgPack

#endif // COO_MATRIX_CLASS_H
