// /////////////////////////////////////////////////////////////////
// VectorWithOpMutableDense.hpp
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

#ifndef VECTOR_WITH_OP_MUTABLE_DENSE_H
#define VECTOR_WITH_OP_MUTABLE_DENSE_H

#include "VectorSpaceSerial.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorMutable.hpp"
#include "DenseLinAlgPack/src/DVectorClass.hpp"
#include "ref_count_ptr.hpp"
#include "ReleaseResource.hpp"

namespace SparseLinAlgPack {

///
/** DVector "Adaptor" subclass for <tt>DenseLinAlgPack::DVectorSlice</tt>
 * or <tt>DenseLinAlgPack::DVector</tt> objects.
 *
 * This class can be used either as a view of a <tt>DenseLinAlgPack::DVectorSlice</tt> object
 * or as a storage type for a <tt>DenseLinAlgPack::DVector</tt> object.
 *
 * To create a storage type with the dimension of \c dim just call the constructor
 * <tt>VectorWithOpMutableDense(dim)</tt> or after construction you can call
 * <tt>this->initialize(dim)</tt>.
 *
 * To simply create a view of a vector, say \c v, without ownership just call
 * <tt>VectorWithOpMutableDense(v(),NULL)</tt> or after construction call
 * <tt>this->initialize(v(),NULL)</tt>.
 *
 * Alternately, \c this can be given a vector with the responsibility to
 * delete any associated memory by calling <tt>this->initialize()</tt>
 * with a <tt>ReleaseResource</tt> object to perform the deallocation.
 *
 * If \c this has been initialized by <tt>this->initialize(dim)</tt> and if
 * the client really needs to get at the <tt>DenseLinAlgPack::DVector</tt> object
 * itself, then it can be obtained as:
 \code
 void f( VectorWithOpMutableDense* v )
     namespace rmp = MemMngPack;
     DVector &_v = *dynamic_cast<rmp::ReleaseResource_ref_count_ptr<DVector>&>(*v.vec_release()).ptr;

 \endcode
 * This is not pretty but it is not supposed to be.  Of course the above function will throw
 * an exception if the <tt>dynamic_cast<></tt> fails.
 */
class VectorWithOpMutableDense
	: virtual public AbstractLinAlgPack::VectorMutable
{
public:

	///
	typedef MemMngPack::ref_count_ptr<
		MemMngPack::ReleaseResource>  release_resource_ptr_t;

	/** @name Constructors/initializers */
	//@{

	///
	/** Calls <tt>this->initialize(dim)</tt>.
	 */
	VectorWithOpMutableDense(
		const size_type                    dim = 0
		);
	///
	/** Calls <tt>this->initialize(v,v_release)</tt>.
	 */
	VectorWithOpMutableDense(
		DVectorSlice                        v
		,const release_resource_ptr_t&     v_release
		);
	///
	/** Call <tt>this->initialize(v,v_release)</tt> with an allocated <tt>DenseLinAlgPack::DVector</tt>
	 * object.
	 */
	void initialize(
		const size_type                    dim
		);
	///
	/** Initialize with a dense vector slice.
	 */
	void initialize(
		DVectorSlice                        v
		,const release_resource_ptr_t&     v_release
		);

	//@}

	/** @name Access */
	//@{
	
	///
	/** Return the non-const dense vector.
	 *
	 * Note that calling this method will result in the vector implementation
	 * being modified.  Therefore, no other methods on \c this object should be
	 * called until the <tt>DVectorSlice</tt> returned from this method is
	 * discarded.
	 *
	 * Note that the underlying implementation calls <tt>this->has_changed()</tt>
	 * before this method returns.
	 */
	DVectorSlice set_vec();
	///
	/** Return a const dense vector.
	 */
	const DVectorSlice get_vec() const;
	///
	/** Return a <tt>ref_count_ptr<></tt> pointer to the object that will
	 * release the associated resource.
	 */
	const release_resource_ptr_t& vec_release() const;

	//@}

	/** @name Overriddenn from Vector */
	//@{

	///
	const VectorSpace& space() const;
	///
	void apply_reduction(
		const RTOpPack::RTOp& op
		,const size_t num_vecs, const Vector** vecs
		,const size_t num_targ_vecs, VectorMutable** targ_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type first_ele, const index_type sub_dim, const index_type global_offset
		) const;
	///
	index_type dim() const;
	///
	value_type get_ele(index_type i) const;
	///
	void get_sub_vector(
		const Range1D& rng, ESparseOrDense sparse_or_dense, RTOp_SubVector* sub_vec ) const;
	///
	void free_sub_vector( RTOp_SubVector* sub_vec ) const;

	//@}

	/** @name Overriddenn from VectorMutable */
	//@{

	///
	void apply_transformation(
		const RTOpPack::RTOp& op
		,const size_t num_vecs, const Vector** vecs
		,const size_t num_targ_vecs, VectorMutable** targ_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type first_ele, const index_type sub_dim, const index_type global_offset
		);
	///
	VectorMutable& operator=(value_type alpha);
	///
	VectorMutable& operator=(const Vector& v);
	///
	VectorMutable& operator=(const VectorMutable& v);
	///
	void set_ele( index_type i, value_type val );
	///
	vec_mut_ptr_t sub_view( const Range1D& rng );
	///
	void get_sub_vector(
		const Range1D& rng, RTOp_MutableSubVector* sub_vec );
	///
	void commit_sub_vector( RTOp_MutableSubVector* sub_vec );
	///
	void set_sub_vector( const RTOp_SubVector& sub_vec );
	///
	void Vp_StMtV(
		value_type                       alpha
		,const GenPermMatrixSlice        &P
		,BLAS_Cpp::Transp                P_trans
		,const Vector              &x
		,value_type                      beta
		);

	//@}

	/// Hack
	VectorWithOpMutableDense* operator&()
	{
		return this;
	}
private:

	// ///////////////////////////////////////
	// Private data members
	
	DVectorSlice               v_;
	release_resource_ptr_t    v_release_;
	VectorSpaceSerial         space_;

	// ///////////////////////////////////////
	// Private member functions

	///
	/** Performs a reduction/transformation operation
	 */
	void apply_op(
		const RTOpPack::RTOp& op
		,const size_t num_vecs,      const Vector**       vecs
		,const size_t num_targ_vecs, VectorMutable**      targ_vecs
		,RTOp_ReductTarget reduct_obj
		,const index_type first_ele, const index_type sub_dim, const index_type global_offset
		) const;

	// Not defined and not to be called
	VectorWithOpMutableDense(const VectorWithOpMutableDense&);
	VectorWithOpMutableDense& operator=(const VectorWithOpMutableDense&);

}; // end class VectorWithOpMutableDense

// //////////////////////////////////////
// Inline members

inline
DVectorSlice
VectorWithOpMutableDense::set_vec()
{
	this->has_changed();
	return v_;
}

inline
const DVectorSlice
VectorWithOpMutableDense::get_vec() const
{
	return v_;
}

inline
const VectorWithOpMutableDense::release_resource_ptr_t&
VectorWithOpMutableDense::vec_release() const
{
	return v_release_;
}

} // end namespace SparseLinAlgPack

#endif // VECTOR_WITH_OP_MUTABLE_DENSE_H
