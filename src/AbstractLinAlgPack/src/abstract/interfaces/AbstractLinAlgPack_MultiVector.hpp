// ///////////////////////////////////////////////////////////
// MultiVector.hpp
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

#ifndef ALAP_MULTI_VECTOR_H
#define ALAP_MULTI_VECTOR_H

#include "MatrixWithOp.hpp"
#include "RTOpPack/src/RTOpCpp.hpp"
#include "ref_count_ptr.hpp"

namespace AbstractLinAlgPack {

///
/** Interface for a collection of non-mutable vectors (multi-vector, matrix).
 *
 * This interface is quite restrictive in that it allows a client to access a
 * matrix by accessing rows, columns and/or diagonals.
 * The vector objects returned from the provided access methods
 * \c row(), \c col() and \c diag() are abstract vectors so there is
 * still good implementation flexibility but many <tt>%MatrixWithOp</tt>
 * implementations will not be able to support this interface.
 *
 * The primary purpose for this interface is to allow for convienent aggregations
 * of column vectors.  Such an orderly arrangement allows for better optimized
 * linear algebra operations such as matrix-matrix multiplication and the solution
 * of linear systems for multiple right hand sides.  Every application area (serial
 * parallel, out-of-core etc.) should be able to define at least one reasonbly
 * efficient implementation of a <tt>%MultiVector</tt> (or a <tt>%MultiVectorMutable</tt>)
 * subclass.
 *
 * The <tt>%MultiVector</tt> interface is derived from the \c MatrixWithOp 
 * interface and therefore a <tt>%MultiVector</tt> can be considered as a matrix
 * which has some interesting implications.  As an extended matrix interface, this
 * is somewhat of a "last resort" interface that allows many matrix operations to
 * have default implementations based on vector operations. None of the linear
 * algebra methods in <tt>%MatrixWithOp</tt> or any of the other matrix interfaces
 * have methods that directly accept <tt>%MultiVector</tt> objects.  However,
 * since <tt>%MultiVector</tt> is derived from <tt>%MatrixWithOp</tt>, a
 * <tt>%MultiVector</tt> object can be used anywere a <tt>%MatrixWithOp</tt>
 * object is accepted.  In fact, many of the default implementations for the
 * linear algebra methods in <tt>%MatrixWithOp</tt> test to see if the matrix
 * arguments support <tt>%MultiVector</tt> (or \c MultiVectorMutable</tt>)
 * and will fail if these interfaces are not supported.
 *
 * Note that only certain kinds of access may be preferred and it is allowed
 * for subclasses to return \c NULL vector objects for some types of access.  For
 * example, a matrix may be naturally oriented by column (for the primary role of
 * as a multi-vector) but row or diagonal access may be very inefficient.  For this
 * reason, the client should call the \c access_by() method which returns a bit
 * field that the client can compare to the constants \c ROW_ACCESS, \c COL_ACCESS
 * and \c DIAG_ACCESS.  The method \c access_by() only returns the types of access
 * that are guarrentted to be efficient, but
 * does not necessarily imply that a type of access is not supported.  For example,
 * <tt>(this->access_by() & ROW_ACCESS) == false</tt> this does not mean that 
 * <tt>this->row(i)</tt> will return \c NULL, it only means that row access will
 * be inefficient.  To determine if a certain type of access is even possible, check
 * the return for \c row(), \c col() and/or \c diag().  For example, if <tt>this->rows(1)</tt>
 * returns \c NULL, then this is a flag that row access for every row is not supported.
 * Diagonal access may be a different story.  For some matrix subclasses, only the
 * center diagonal my be easily accessable in which case \c diag(0) may return
 * <tt>!= NULL</tt> but \c diag(k) for all <tt>k != 0</tt> may return \c NULL.
 *
 * Note that since, this interface is derived from \c MatrixWithOp that it must
 * support the methods \c space_rows() and \c space_cols().  This does not imply
 * however that either of the access methods \c row() or \c col() must return
 * non-<tt>NULL</tt>.
 *
 * Examples of matrix implementations that can support this interface are a dense
 * BLAS compatible matrix (\c ROW_ACCESS, \c COL_ACCESS and \c DIAG_ACCESS), a
 * compressed column sparse matrix (\c COL_ACCESS only), a compressed row sparse matrix
 * (\c ROW_ACCESS only), a unordered sparse matrix with the diagonal explicity stored
 * (<tt>diag(0) != NULL</tt>) etc.
 *
 * Another very powerfull feature of this interface is the ability to apply
 * reduction/transformation operators over a sub-set of rows and columns in a set
 * of multi-vector objects.  The behavior is identical as if the client extracted
 * the rows or columns in a set of multi-vectors and called
 * \c VectorWithOp::apply_reduction() or \c VectorWithOpMutable::apply_transformation()
 * itself.  However, the advantage of using the multi-vector methods is that there
 * may be greater opertunity to exploit parallelism.  Also, the intermediate reduction
 * objects over a set of rows or columns can be reduced by a secondary reduction object.
 *
 * ToDo: Finish documentation!
 */
class MultiVector : virtual public MatrixWithOp {
public:

	///
	typedef int  access_by_t;
	///
	enum {
		ROW_ACCESS    = 0x1 ///< 
		,COL_ACCESS   = 0x2 ///<
		,DIAG_ACCESS  = 0x4 ///<
	};
	///
	enum EApplyBy {
		APPLY_BY_ROW        ///<
		,APPLY_BY_COL       ///<
	};
	///
	typedef MemMngPack::ref_count_ptr<const VectorWithOp>   vec_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<const MultiVector>    multi_vec_ptr_t;

	/** @name Provide row, column and diagonal access as non-mutable vectors */
	//@{

	///
	/** Return a bit field for the types of access that are the most convenient.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return & COL_ACCESS || return & ROW_ACCESS || return & DIAG_ACCESS</tt>
	 * </ul>
	 */
	virtual access_by_t access_by() const = 0;

	///
	/** Get a non-mutable column vector.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->access_by() & COL_ACCESS</tt>] <tt>return.get() != NULL</tt>
	 * <li> [<tt>return.get() != NULL</tt>] <tt>space_cols().is_compatible(return->space()) == true</tt>
	 * </ul>
	 */
	virtual vec_ptr_t col(index_type j) const = 0;
	///
	/** Get a non-mutable row vector.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->access_by() & ROW_ACCESS</tt>] <tt>return.get() != NULL</tt>
	 * <li> [<tt>return.get() != NULL</tt>] <tt>space_rows().is_compatible(return->space()) == true</tt>
	 * </ul>
	 */
	virtual vec_ptr_t row(index_type i) const = 0;
	///
	/** Get a non-mutable diagonal vector.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->access_by() & DIAG_ACCESS</tt>] <tt>return.get() != NULL</tt>
	 * </ul>
	 */
	virtual vec_ptr_t diag(int k) const = 0;

	//@}

	/** @name Sub-view methods */
	//@{

	///
	/** Returns a sub-view of the multi vector.
	 *
	 * ToDo: Finish documentation!
	 *
	 * The default implementation returns a \c MultiVectorSubView object for
	 * any valid arbitary sub-view.
	 */
	virtual multi_vec_ptr_t mv_sub_view(const Range1D& row_rng, const Range1D& col_rng) const;
	
	///
	/** Inlined implementation calls <tt>this->mv_sub_view(Range1D(rl,ru),Range1D(cl,cu))</tt>.
	 */
	multi_vec_ptr_t mv_sub_view(
		const index_type& rl, const index_type& ru
		,const index_type& cl, const index_type& cu
		) const;

	//@}

	/** @name Collective apply_reduction() methods */
	//@{

	///
	/** Apply a reduction/transformation operator row by row, or column by column and return an array
	 * of the reduction objects.
	 *
	 * Preconditions:<ul>
	 * <li> [<tt>apply_by == APPLY_BY_COL</tt>] <tt>(this->access_by() & COL_ACCESS) == true)</tt> (throw <tt>???</tt>)
	 * <li> [<tt>apply_by == APPLY_BY_ROW</tt>] <tt>(this->access_by() & ROW_ACCESS) == true)</tt> (throw <tt>???</tt>)
	 * <li> ToDo: Finish!
	 * </ul>
	 *
	 * The default implementation calls \c this->apply_op().
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void apply_reduction(
		EApplyBy apply_by, const RTOpPack::RTOp& primary_op
		,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
		,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
		,RTOp_ReductTarget reduct_objs[]
		,const index_type primary_first_ele   = 1, const index_type primary_sub_dim   = 0, const index_type primary_global_offset = 0
		,const index_type secondary_first_ele = 1, const index_type secondary_sub_dim = 0
		) const;

	///
	/** Apply a reduction/transformation operator row by row, or column by column and reduce the intermediate
	 * reduction objects into one reduction object.
	 *
	 * Preconditions:<ul>
	 * <li> [<tt>apply_by == APPLY_BY_COL</tt>] <tt>(this->access_by() & COL_ACCESS) == true)</tt> (throw <tt>???</tt>)
	 * <li> [<tt>apply_by == APPLY_BY_ROW</tt>] <tt>(this->access_by() & ROW_ACCESS) == true)</tt> (throw <tt>???</tt>)
	 * <li> ToDo: Finish!
	 * </ul>
	 *
	 * The default implementation calls \c this->apply_op().
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void apply_reduction(
		EApplyBy apply_by, const RTOpPack::RTOp& primary_op, const RTOpPack::RTOp& secondary_op
		,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
		,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type primary_first_ele   = 1, const index_type primary_sub_dim   = 0, const index_type primary_global_offset = 0
		,const index_type secondary_first_ele = 1, const index_type secondary_sub_dim = 0
		) const;

	//@}

	/** @name Overridden from MatrixWithOp */
	//@{

	///
	/** Returns <tt>this->mv_sub_view(row_rng,col_rng)</tt> casted to a MatrixWithOp.
	 */
	mat_ptr_t sub_view(const Range1D& row_rng, const Range1D& col_rng) const;

	///
	/** Provides a specialized implementation for <tt>mwo_rhs1</tt> of type <tt>MatrixSymDiagonal</tt>.
	 *
	 * @return Returns <tt>true</tt> and implements the operation if
	 * <tt>dynamic_cast<MatrixSymDiagonal>(&mwo_rhs1) != NULL
	 * && op(*this).access_by() =& MultiVector::COL_ACCESS
	 * && (mvm_lhs = dynamic_cast<MultiVectorMutable*>(&mwo_lhs)) != NULL
	 * && mvm_lhs->access_by() & MultiVector::COL_ACCESS</tt>.
	 * Otherwise, this function returns <tt>false</tt> and does not implement the operation.
	 * or <tt>dynamic_cast<const MatrixSymDiagonal>(&mwo_rhs1) != NULL</tt>.
	 *
	 * The default implementation relies on column access of <tt>op(*this)</tt>
	 * and <tt>mwo_lhs</tt> to implement this method.
	 */
	bool Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		,const MatrixWithOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		,BLAS_Cpp::Transp trans_rhs2
		,value_type beta ) const;

	///
	/** Provides a specialized implementation for <tt>mwo_rhs2</tt> of type <tt>MatrixSymDiagonal</tt>.
	 *
	 * @return Returns <tt>true</tt> and implements the operation if
	 * <tt>dynamic_cast<MatrixSymDiagonal>(&mwo_rhs1) != NULL
	 * && op(*this).access_by() =& MultiVector::ROW_ACCESS
	 * && (mvm_lhs = dynamic_cast<MultiVectorMutable*>(&mwo_lhs)) != NULL
	 * && mvm_lhs->access_by() & MultiVector::ROW_ACCESS</tt>.
	 * Otherwise, this function returns <tt>false</tt> and does not implement the operation.
	 * or <tt>dynamic_cast<const MatrixSymDiagonal>(&mwo_rhs1) != NULL</tt>.
	 *
	 * The default implementation relies on row access of <tt>op(*this)</tt>
	 * and <tt>mwo_lhs</tt> to implement this method.
	 */
	bool Mp_StMtM(
		MatrixWithOp* mwo_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixWithOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
		,value_type beta ) const;

	//@}

protected:

	/** @name Combined implementations for apply_reduction() and apply_transformation() methods */
	//@{

	///
	/** Apply a reduction/transformation operator row by row, or column by column and return an array
	 * of the reduction objects.
	 *
	 * Note that \c this is already one of the vecs[] or targ_vecs[] arguments.
	 *
	 * Preconditions:<ul>
	 * <li> ToDo: Finish!
	 * </ul>
	 *
	 * The default implementation relys on \c row() or \c col() access methods to apply the operators.
	 */
	virtual void apply_op(
		EApplyBy apply_by, const RTOpPack::RTOp& primary_op
		,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
		,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
		,RTOp_ReductTarget reduct_objs[]
		,const index_type primary_first_ele   = 1, const index_type primary_sub_dim   = 0, const index_type primary_global_offset = 0
		,const index_type secondary_first_ele = 1, const index_type secondary_sub_dim = 0
		) const;

	///
	/** Apply a reduction/transformation operator row by row, or column by column and reduce the intermediate
	 * reduction objects into one reduction object.
	 *
	 * Note that \c this is already one of the vecs[] or targ_vecs[] arguments.
	 *
	 * Preconditions:<ul>
	 * <li> ToDo: Finish!
	 * </ul>
	 *
	 * The default implementation relys on \c row() or \c col() access methods to apply the operators.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void apply_op(
		EApplyBy apply_by, const RTOpPack::RTOp& primary_op, const RTOpPack::RTOp& secondary_op
		,const size_t num_multi_vecs,      const MultiVector**   multi_vecs
		,const size_t num_targ_multi_vecs, MultiVectorMutable**  targ_multi_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type primary_first_ele   = 1, const index_type primary_sub_dim   = 0, const index_type primary_global_offset = 0
		,const index_type secondary_first_ele = 1, const index_type secondary_sub_dim = 0
		) const;

	//@}

private:
	
#ifdef DOXYGEN_COMPILE
	VectorWithOp         *rows;
	VectorWithOp         *columns;
	VectorWithOp         *diagonals;
#endif	

}; // end class MultiVector

// //////////////////////////////////////////////////
// Inlined methods for MultiVector

inline
MultiVector::multi_vec_ptr_t
MultiVector::mv_sub_view(
	const index_type& rl, const index_type& ru
	,const index_type& cl, const index_type& cu
	) const
{
	return this->mv_sub_view(Range1D(rl,ru),Range1D(cl,cu));
}

} // end namespace AbstractLinAlgPack

#endif // ALAP_MULTI_VECTOR_H
