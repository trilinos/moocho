// /////////////////////////////////////////////////////////
// GenPermMatrixSlice.h
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

#ifndef GEN_PERM_MATRIX_SLICE_H
#define GEN_PERM_MATRIX_SLICE_H

#include "GenPermMatrixSliceIterator.h"

namespace AbstractLinAlgPack {

///
/** Concrete matrix type to represent general permutation (mapping) matrices.
  *
  * These are matrices who's rows or columns represent eta vectors
  * (i.e. only one nonzero element with the value 1).  These matrices
  * can be rectangular and have one or more zero rows & columns.  Therefore, these
  * matrices can be used to represent gathering and scattering operations
  * on other vectors and matrices.
  *
  * This is only a view type.  The client specifies the mapping arrays and then
  * this class provides a clean encapsulation for the mapping.  Objects of this
  * type can also represent the identity matrix which is constructed with
  * the initialize_identity(...) function.
  *
  * The default copy constructor is allowd but the default 
  * assignment operator function is not.
  */
class GenPermMatrixSlice {
public:

	// ////////////////////////////////////////
	/** @name Public types.
	  */
	//@{

	///
	enum EIdentityOrZero { IDENTITY_MATRIX, ZERO_MATRIX };

	///
	typedef GenPermMatrixSliceIteratorPack::EOrderedBy		EOrderedBy;

	///
	typedef GenPermMatrixSliceIteratorPack::row_col_iterator<const index_type>
											const_iterator;
	///
	typedef	ptrdiff_t						difference_type;

	//@}

	/// Construct to an uninitialzied, unsized matrix.
	GenPermMatrixSlice();

	/// Construct to a matrix intialized to identity or zero (see initialize(,,,)).
	GenPermMatrixSlice( index_type rows, index_type cols, EIdentityOrZero type );

	///
	/** Initialize an identity or zero permutation.
	 *
	 * If #type == IDENTITY_MATRIX# then after this function is called #this# will
	 * represent #Q = [ I; 0 ]# if #rows > cols# or #Q = [ I, 0 ]#
	 * if #rows < cols# or #Q = I# if #rows == cols#.  If #type == ZERO_MATRIX]#
	 * then #this# will represent a #rows# x #cols# zero matrix.
	 *
	 * Postconditions:\begin{itemize}
	 * \item #this->rows() == rows#
	 * \item #this->cols() == cols#
	 * \item [#type == IDENTITY_MATRIX#] #this->nz() == min(rows,cols)#
	 * \item [#type == ZERO_MATRIX#]     #this->nz() == 0#
	 * \item [#type == IDENTITY_MATRIX#] #this->is_identity() == true#
	 * \item [#type == ZERO_MATRIX#]     #this->is_identity() == false#
	 * \item #this->ordered_by() == BY_ROW_AND_COL#
	 * \end{itemize}
	 */
	void initialize( index_type rows, index_type cols, EIdentityOrZero type );

	///
	/** Initialize.
	  *
	  * This function sets up general permutation view.
	  * The client must setup the arrays row_i[] and col_j[] to define
	  * the mapping.  If nz == 0 then row_i and col_j can be NULL.
	  * Otherwise, row_i != NULL and col_j != NULL.  If nz > 0 and
	  * ordered_by == BY_ROW then row_i[k] must be sorted in assending order
	  * else if ordered_by == BY_COL or col_j[k] must be sorted in assnding order
	  * else if ordered_by == UNORDERED then no ordering for row_i[k] or
	  * col_j[k] is required.  If nz == 1 then the value of ordered_by
	  * is not really significant at it is ordered by row and by column.
	  * 
	  * It is required that if nz > 0 then:
	  * 1 <= row_i[k] + row_off <= rows, for k = 1...nz
	  * 1 <= col_j[k] + col_off <= cols, for k = 1...nz
	  * 
	  * All of these preconditions will be checked if test_setup == true.
	  *
	  * After setup, the memory pointed to by row_i[] and col_j[] must not
	  * be altered since this object does not make an independent copy
	  * of this data.
	  * 
	  * After construction the nonzero elements of this matrix are:
	  * M(row_i[k]+row_off,col_j[k]+col_off) = 1.0, for k = 1...nz.
	  *
	  *	@param	rows	[I]	Number of rows in matrix
	  *	@param	cols	[I]	Number of columns in matrix
	  *	@param	nz		[I]	Number of nonzero elements in the matrix
	  *	@param	row_off	[I]	Row offsets for row_i[]
	  *	@param	col_off	[I]	Column offsets for col_i[]
	  *	@param	ordered_by
	  *					[I]	The ordering of the nonzero elements
	  *	@param	row_i	[I]	Array (size nz):  If nz == 0 then
	  *						row_i can be NULL
	  *	@param	col_j	[I]	Array (size nz): If nz == 0 then
	  *						col_j can be NULL
	  *	@param	test_setup
	  *					[I] If true then all of the preconditions for
	  *						the input arguments will be checked.
	  */
	void initialize(
		index_type			rows
		,index_type			cols
		,index_type			nz
		,difference_type	row_off
		,difference_type	col_off
		,EOrderedBy			ordered_by
		,const index_type	row_i[]
		,const index_type	col_j[]
		,bool				test_setup = false
		);

	///
	/** Initialize and sort.
	  *
	  * This is the same as the initialize(...) function except that
	  * this function will actually sort the entries by row or
	  * by column or not at all..
	  *
	  * ToDo: Finish documentation.
	  *
	  *	@param	rows	[I]	Number of rows in matrix
	  *	@param	cols	[I]	Number of columns in matrix
	  *	@param	nz		[I]	Number of nonzero elements in the matrix
	  *	@param	row_off	[I]	Row offsets for row_i[]
	  *	@param	col_off	[I]	Column offsets for col_i[]
	  *	@param	ordered_by
	  *					[I]	The ordering of the nonzero elements
	  *	@param	row_i	[I/O]	Array (size nz):  If nz == 0 then
	  *						row_i can be NULL.  On output it will be
	  *						sorted according to ordered_by
	  *	@param	col_j	[I/O]	Array (size nz): If nz == 0 then
	  *						col_j can be NULL.  On output it will be
	  *						sorted according to ordered_by
	  *	@param	test_setup
	  *					[I] If true then all of the preconditions for
	  *						the input arguments will be checked.
	  */
	void initialize_and_sort(
		index_type			rows
		,index_type			cols
		,index_type			nz
		,difference_type	row_off
		,difference_type	col_off
		,EOrderedBy			ordered_by
		,index_type			row_i[]
		,index_type			col_j[]
		,bool				test_setup = false
		);
		
	///
	/** Bind the view of another GenPermMatrixSlice object.
	  *
	  * After construction these objects will share points to
	  * the same row_i[] and col_j[] arrays.
	  */
	void bind( const GenPermMatrixSlice& gpms );

	///
	index_type rows() const;
	///
	index_type cols() const;
	///
	index_type nz() const;
	///
	EOrderedBy ordered_by() const;
	///
	bool is_identity() const;

	///
	/** Lookup the ith row index for the nonzero entry in the jth column
	 * if it exists.
	 *
	 * This function will return 0 if the index is not found.  If 
	 * #this->ordered_by() == BY_COL || this->ordered_by() == BY_ROW_AND_COL#
	 * then this function will be executed in O(log(#this->nz()#)) time.
	 * Otherwise it will execute in O(#this->nz()#) time.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #(1 <= col_j && col_j <= this->cols())# (throw #std::out_of_range#)
	 * \end{itemize}
	 *
	 * Postconditions:\begin{itemize}
	 * \item #(1 <= return && return <= this->rows()) || return == 0#
	 * \end{itemize}
	 */
	index_type lookup_row_i(index_type col_j) const;

	///
	/** Lookup the jth column index for the nonzero entry in the ith row
	 * if it exists.
	 *
	 * This function will return 0 if the index is not found.  If 
	 * #this->ordered_by() == BY_ROW || this->ordered_by() == BY_ROW_AND_COL#
	 * then this function will be executed in O(log(#this->nz()#)) time.
	 * Otherwise it will execute in O(#this->nz()#) time.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #(1 <= row_i && row_i <= this->rows())# (throw #std::out_of_range#)
	 * \end{itemize}
	 *
	 * Postconditions:\begin{itemize}
	 * \item #(1 <= return && return <= this->cols()) || return == 0#
	 * \end{itemize}
	 */
	index_type lookup_col_j(index_type row_i) const;

	/** @name Iterator Access.
	  *
	  * These functions return forward iterators for accessing which
	  * rows and columns each nonzero 1 entry is located at.
	  * If this->is_identity() == true then these iterators are obviously
	  * unnecessary and will throw exceptions.
	  /begin{verbatim}
	for( GenPermMatrixSlice::const_iterator itr = gpms.begin(); itr != gpms.end(); ++itr )
	{
		std::cout << "row_i = " << itr->row_i();
		std::cout << "col_j = " << itr->col_j();
	}
	  /end{verbatim}
	  *
	  * You can also take a difference between iterators.
	  */
	//@{

	///
	const_iterator begin() const;

	///
	const_iterator end() const;

	//@}

	///
	/** Create a submatrix by row, by column.
	  *
	  * If nz > 1 and this->ordered_by() == BY_ROW then ordered_by
	  * must also equal BY_ROW or if this->ordered_by() == BY_COL
	  * then ordered_by must also equal BY_COL
	  * or an std::logic_error exception will be thrown.  If nz == 1, then
	  * obviously the nozeros are ordered by row and by column.
	  * This function should not be called if this->is_identity() == true.
	  * 
	  * The argument rng must be explicitly sized (rng.full_range() != true) or the
	  * exception std::logic_error will be thrown.  Also, the range argument
	  * must obey rng.lbound() >= 1, and rng.ubound() <= this->rows()
	  * if ordered_by == BY_ROW and rng.ubound() <= this->cols()
	  * if ordered_by == BY_COL.  The argument ordered_by == UNORDERED
	  * is not allowed and this function can not be called
	  * if this->ordered_by() == UNORDERED.  This operation just does not
	  * make any sense.
	  * 
	  * The returned submatrix will contain all the entries for the designated
	  * rows if ordered_by == BY_ROW or columns if ordered_by == BY_CO.
	  * 
	  * ToDo: Spell out the behavior of this operation more carefully.
	  */
	const GenPermMatrixSlice create_submatrix( const Range1D& rng
		, EOrderedBy ordered_by ) const;

private:

	// //////////////////////////////
	// Private data members

	index_type			rows_;
	index_type			cols_;
	index_type			nz_;
	difference_type		row_off_;
	difference_type		col_off_;
	EOrderedBy			ordered_by_;
	const index_type		*row_i_;
	const index_type		*col_j_;

	// //////////////////////////////
	// Private member functions

	// Validate the input data (not the ordering!)
	static void validate_input_data(
		index_type			rows
		,index_type			cols
		,index_type			nz
		,difference_type	row_off
		,difference_type	col_off
		,EOrderedBy			ordered_by
		,const index_type	row_i[]
		,const index_type	col_j[]
		,std::ostringstream &omsg
		);

	///
	void validate_not_identity() const;
	
	/// not defined and not to be called
	GenPermMatrixSlice& operator=( const GenPermMatrixSlice& );

};	// end class GenPermMatrixSlice

// //////////////////////////////////////////////////////////
// Inline members for GenPermMatrixSlice

inline
GenPermMatrixSlice::GenPermMatrixSlice( index_type rows, index_type cols, EIdentityOrZero type )
{
	initialize(rows,cols,type);
}

inline
index_type GenPermMatrixSlice::rows() const
{
	return rows_;
}

inline
index_type GenPermMatrixSlice::cols() const
{
	return cols_;
}

inline
index_type GenPermMatrixSlice::nz() const
{
	return nz_;
}

inline
bool GenPermMatrixSlice::is_identity() const
{
	return nz_ > 0 && row_i_ == NULL;
}

inline
GenPermMatrixSlice::EOrderedBy GenPermMatrixSlice::ordered_by() const
{
	return ordered_by_;
}

}	// end namespace AbstractLinAlgPack

#endif   // GEN_PERM_MATRIX_SLICE_H
