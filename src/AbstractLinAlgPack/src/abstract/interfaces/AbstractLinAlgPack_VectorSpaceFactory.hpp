// ///////////////////////////////////////////////////////////////
// VectorSpaceFactory.hpp
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

#ifndef VECTOR_SPACE_FACTORY_H
#define VECTOR_SPACE_FACTORY_H

#include "AbstractLinAlgPack/src/AbstractLinAlgPackTypes.hpp"
#include "ref_count_ptr.hpp"

namespace AbstractLinAlgPack {

///
/** Abstract interface for objects that can create vector spaces of a specified dimension.
 *
 * ToDo: Finish documentation!
 */
class VectorSpaceFactory
{
public:

	///
	typedef MemMngPack::ref_count_ptr<const InnerProduct>   inner_prod_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<const VectorSpace>    space_ptr_t;

	/** @name Constructors / initializers */
	//@{

	/// Calls \c inner_prod()
	VectorSpaceFactory( const inner_prod_ptr_t& inner_prod = MemMngPack::null );

	///
	/** Initialize with an inner product object that will be given to vector.
	 *
	 * @param  inner_prod  [in] Smart pointer to inner product strategy object.
	 *                     If <tt>inner_prod.get()==NULL</tt> then an
	 *                     \c InnerProductDot object will be used instead.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>inner_prod.get() != NULL</tt>] <tt>this->inner_prod().get() == inner_prod.get()</tt>
	 * <li> [<tt>inner_prod.get() == NULL</tt>] <tt>dynamic_cast<InnerProductDot*>(this->inner_prod().get()) != NULL</tt>
	 * </ul>
	 */
	virtual void inner_prod( const inner_prod_ptr_t& inner_prod );

	///
	/** Return the smart pointer to the inner product strategy object.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get() != NULL</tt>
	 * </ul>
	 */
	virtual const inner_prod_ptr_t inner_prod() const;

	//@}

	/** @name Pure virtual functions that must be overridden */
	//@{

	///
	/** Create a vector space of the given dimension.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get() != NULL</tt>
	 * <li> <tt>return->dim() == dim</tt>
	 * <li> [<tt>this->inner_prod().get() != NULL</tt>] <tt>this->inner_prod().get() == return->inner_prod().get()</tt>
	 * <li> [<tt>this->inner_prod().get() == NULL</tt>] <tt>dynamic_cast<InnerProductDot*>(return->inner_prod().get()) != NULL</tt>
	 * </ul>
	 *
	 * @return  Returns a smart reference counted pointer to a dynamically
	 * allocated vector space object that can be used to create vector.
	 */
	virtual space_ptr_t create_vec_spc(index_type dim) const = 0;

	//@}
	
private:
#ifdef DOXYGEN_COMPILE
	InnerProduct       *inner_prod;
#else
	inner_prod_ptr_t   inner_prod_;
#endif

}; // end class VectorSpaceFactory

} // end namespace AbstractLinAlgPack

#endif  // VECTOR_SPACE_FACTORY_H
