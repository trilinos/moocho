// //////////////////////////////////////////////////////////
// MatrixLoadSparseElements.h
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

#ifndef MATRIX_LOAD_SPARSE_FORTRAN_COMPATIBLE_ELEMENTS_H
#define MATRIX_LOAD_SPARSE_FORTRAN_COMPATIBLE_ELEMENTS_H

#include "AbstractLinAlgPack/include/MatrixBase.h"

namespace SparseLinAlgPack {

///
/** Mix-in interface for loading nonzero elements into a sparse matrix data structure.
 *
 * The formats supported are:
 *
 * Coordiante:
 \verbatim

 		Aval[k], Arow[k], Acol[k], k = 0..num_nonzeros(...)-1
 \endverbatim
 * Compressed Row (Column):
 \verbatim

		Aval[k], Acol[k], k = 0..num_nonzeros(...)-1
		Arow_start[j], j = 0..rows()-1
 \endverbatim
 * This is meant to be the do-all interface for clients to use to add nonzero elements
 * for sparse matrices.
 *
 * ToDo: Discuss element uniqueness!
 *
 * ToDo: Finish documentation!
 */
class MatrixLoadSparseElements
	: public virtual AbstractLinAlgPack::MatrixBase
{
public:

	///
	/** Resize the matrix and reserve space for nonzero elements to be added.
	 *
	 * All of the nonzeros in the current matrix are discarded and we start fresh.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void reinitialize(
		size_type  rows
		,size_type cols
		,size_type max_nz
		) = 0;

	///
	/** Reinitialize internal counter to load new nonzero values.
	 *
	 * The row and column index arrays are preserved from the last setup and here
	 * the client only wants to set the nonzero values for the same matrix
	 * structure.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void reset_to_load_values() = 0;

	///
	/** Get pointers to buffers to add nonzero elements.
	 *
	 * @param  max_nz_load
	 *                  [in] Maximum number of nonzero elements that will be set
	 *                  in the returned buffers.
	 * @param  val      [out] On output <tt>*val</tt> is set to a pointer to an contiguous array
	 *                  of memory of at least \c max_nz_load entries for which the values of the
	 *                  nonzero elements to add are to be set.
	 * @param  row_i    [out] On output <tt>*row_i</tt> is set to a pointer to an contiguous array
	 *                  of memory of at least \c max_nz_load entries for which the row indexes of the
	 *                  nonzero elements to add are to be set.  If <tt>row_i == NULL</tt> then
	 *                  no buffer is allocated for the row indexes.
	 * @param  col_j    [out] On output <tt>*col_J</tt> is set to a pointer to an contiguous array
	 *                  of memory of at least \c max_nz_load entries for which the column indexes of the
	 *                  nonzero elements to add are to be set.  If <tt>col_j == NULL</tt> then
	 *                  no buffer is allocated for the column indexes.
	 *
	 * Preconditions:<ul>
	 * <li> If <tt>reset_to_load_values()</tt> was called to setup nonzero elements, then
	 *      \c row_i and \c col_j must be \c NULL or a <tt>std::logic_error</tt> exception
	 *      is thrown.
	 * </ul>
	 *
	 * After entries in the arrays \c (*val)[], \c (*row_i)[] and \c (*col_j)[] are set, the client must
	 * call <tt>this->commit_load_nonzeros_buffers()</tt> to commit the nonzero entries that are set
	 * in these buffers.
	 */
	virtual void get_load_nonzeros_buffers(
		size_type      max_nz_load
		,value_type    **val
		,index_type    **row_i
		,index_type    **col_j
		) = 0;

	///
	/** Commit nonzeros in buffers obtained from \c get_load_nonzeros_buffers().
	 *
	 * @param  nz_commit
	 *                  [in] Number of nonzero elements to be loaded fron the buffers.  Note that
	 *                  <tt>nz_commit <= max_nz_load</tt> on the previous call to <tt>get_load_nonzeros_buffers()</tt>.
	 * @param  val      [in/out] On input <tt>(*val)[]</tt> contains an array of \c nz_commit value entries for
	 *                  nonzero elements to load.  This must point to the same buffer returned from the last call to
	 *                  \c get_load_nonzero_buffers(). On output <tt>*val</tt> is set to \c NULL>
	 * @param  row_i    [in/out] On input <tt>(*row_i)[]</tt> contains an array of \c nz_commit row index entries for
	 *                  nonzero elements to load.  This must point to the same buffer returned from the last call to
	 *                  \c get_load_nonzero_buffers(). If <tt>row_i == NULL</tt> then no row indexes are set.  Here
	 *                  it is assumed that the row indexes from a previous load have already been set.
	 *                  On output <tt>*row_i</tt> is set to \c NULL>
	 * @param  col_j    [in/out] On input <tt>(*col_j)[]</tt> contains an array of \c nz_commit column index entries for
	 *                  nonzero elements to load.  This must point to the same buffer returned from the last call to
	 *                  \c get_load_nonzero_buffers().  If <tt>col_J == NULL</tt> then no column indexes are set.  Here
	 *                  it is assumed that the column indexes from a previous load have already been set.
	 *                  On output <tt>*col_j</tt> is set to \c NULL>
	 *
	 * Preconditions:<ul>
	 * <li> If <tt>reset_to_load_values()</tt> was called to setup nonzero elements, then
	 *      \c row_i and \c col_j must be \c NULL or a <tt>std::logic_error</tt> exception
	 *      is thrown.
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void commit_load_nonzeros_buffers(
		size_type      nz_commit
		,value_type    **val
		,index_type    **row_i
		,index_type    **col_j
		) = 0;

	///
	/** To be called when the matrix construction is finally finished after all
	 * of the nonzero entries have been added.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->nz()</tt> returns the sum of all of <tt>nz_commit</tt> in all previous calls to
	 *      <tt>commit_load_nonzeros_buffers()</tt> since the last call to <tt>reinitialize()</tt>.
	 * </ul>
	 */
	virtual void finish_construction() = 0;

};	// end class MatrixLoadSparseElements

}	// end namespace SparseLinAlgPack 

#endif	// MATRIX_LOAD_SPARSE_FORTRAN_COMPATIBLE_ELEMENTS_H
