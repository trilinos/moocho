// ///////////////////////////////////////////////////////////////////////
// NLPObjGradient.h
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

#ifndef NLP_OBJ_GRADIENT_H
#define NLP_OBJ_GRADIENT_H

#include "NLP.h"

namespace NLPInterfacePack {
///
/** %NLP interface class that adds gradient information for the objective function {abstract}.
 *
 * <b>Overview:</b>
 *
 * This class adds the ability to compute the gradient of the objective function
 * \c Gf(x) to the basic information given in the \c NLP interface class.  Note that
 * \c Gf is in the vector space \c space_x().
 *
 * <b>Client Usage:</b>
 *
 * As with the <tt>NLP</tt> base interface, the <tt>initialize()</tt> method must be called before
 * the %NLP object can be used.   The method <tt>set_Gf()</tt> is used to set a pointer to a vector
 * to update when the gradient of the objective \c Gf is computed when <tt>calc_Gf()</tt> is called.
 *
 * The number of evaluations of \c Gf using <tt>calc_Gf()</tt> is returned by <tt>num_Gf_evals()</tt>.
 * 
 * <b>Subclass developer's notes:</b>
 *
 * <A NAME="must_override"></A>
 * In addition to the methods that must be overridden by the <tt>NLP</tt> interface
 * (<A HREF="classNLPInterfacePack_1_1NLP.html#must_override">see</A>) the following methods
 * must also be overridden: <tt>imp_calc_Gf()</tt>.
 *
 * <A NAME="should_override"></A>
 * In addition to the methods that should be overridden from <tt>%NLP</tt> by most subclasses
 * (<A HREF="classNLPInterfacePack_1_1NLP.html#should_override">see</A>), the following
 * additional methods should be overridden: \c initialize().
 *
 * The following methods should never have to be overridden by most subclasses except in some very
 * strange situations: \c set_Gf(), \c get_Gf(), \c Gf(), \c num_Gf_evals().
 */
class NLPObjGradient : virtual public NLP {
public:

	/** @name Constructors */
	//@{

	/// Initialize to no reference set to calculation quanities
	NLPObjGradient();

	//@}

	/** @name NLP initialization */
	//@{

	///
	/** Initialize the NLP for its first use.
	  *
	  * This function implementation should be called by subclass implementations
	  * in order to reset counts for \c f(x), \c c(x), \c h(x) and \c Gf(x) evaluations.
	  * This implementation calls <tt>this->NLP::initialize()</tt>
	  *
	  * Postconditions:<ul>
	  * <li> See <tt>NLP::initialize()</tt>
	  * <li> <tt>this->num_Gf_evals() == 0</tt>
	  * </ul>
	  */
	void initialize();

	//@}

	/** @name <<std aggr>> members for the gradient of the objective function Gf(x) */
	//@{

	///
	/** Set a pointer to a vector to be updated when <tt>this->calc_Gf()</tt> is called.
	 *
	 * @param  Gf  [in] Pointer to gradient vector.  May be \c NULL.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> [<tt>Gf != NULL</tt>] <tt>Gf->space().is_compatible(*this->space_x()) == true</tt>
	 *      (throw <tt>VectorBase::IncompatibleVectors</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get_Gf() == Gf</tt>
	 * </ul>
	 */
	virtual void set_Gf(VectorWithOpMutable* Gf);
	///
	/** Return pointer passed to <tt>this->set_Gf()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual VectorWithOpMutable* get_Gf();
	///
	/** Returns non-<tt>const</tt> <tt>*this->get_Gf()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_Gf() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual VectorWithOpMutable& Gf();
	///
	/** Returns <tt>const</tt> <tt>*this->get_Gf()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_Gf() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual const VectorWithOp& Gf() const;

	//@}

	/** @name Calculation Members */
	//@{

	///
	/** Update the vector for \c Gf at the point \c x and put it in the stored reference.
	 *
	 * @param  x     [in] Point at which to calculate the gradient of the objective <tt>Gf(x)</tt>.
	 * @param  newx  [in] (default \c true) If \c true, the values in \c x are the same as
	 *               the last call to a <tt>this->calc_*(x,newx)</tt> member.
	 *               If \c false, the values in \c x are not the same as the last call to a
	 *               <tt>this->calc_*(x,newx)</tt> member.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorBase::IncompatibleVectors</tt>)
	 * <li> <tt>this->get_Gf() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->Gf()</tt> is updated to \c Gf(x)
	 * </ul>
	 *
	 * If <tt>set_mult_calc(true)</tt> was called then referenced storage for \c f and/or \c c and/or \c h
	 * may also be updated but are not guaranteed to be.  But no other quanities from possible subclasses are allowed
	 * to be updated as a side effect (i.e. no higher order derivatives).
	 */ 
	virtual void calc_Gf(const VectorWithOp& x, bool newx = true) const;

	//@}

	/** @name Function evaluation counts. */
	//@{

	///
	/** Objective gradient evaluations count.
	 *
	 * This function can be called to find out how many evaluations
	 * \c this->calc_Gf() the client requested since \c this->initialize() was called.
	 */
	virtual size_type num_Gf_evals() const;
	
	//@}
	
	///
	/** Struct for gradient (objective), objective and constriants (pointers)
	 */
	struct ObjGradInfo {
		///
		ObjGradInfo()
			: Gf(NULL), f(NULL), c(NULL)
		{}
		///
		ObjGradInfo( VectorWithOpMutable* Gf_in, const ZeroOrderInfo& first_order_info_in )
			: Gf(Gf_in), f(first_order_info_in.f), c(first_order_info_in.c), h(first_order_info_in.h)
		{}
		/// Pointer to gradient of objective function <tt>Gf</tt> (may be NULL if not set)
		VectorWithOpMutable*       Gf;
		/// Pointer to objective function <tt>f</tt> (may be NULL if not set)
		value_type*                f;
		/// Pointer to constraints residule <tt>c</tt> (may be NULL if not set)
		VectorWithOpMutable*       c;
		/// Pointer to constraints residule <tt>h</tt> (may be NULL if not set)
		VectorWithOpMutable*       h;
	}; // end struct ObjGradInfo

protected:

	/// Return objective gradient and zero order information.
	const ObjGradInfo obj_grad_info() const;

	/** @name Protected methods to be overridden by subclasses */
	//@{

	///
	/** Overridden to compute f(x) and perhaps c(x) (if multiple calculaiton = true).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
	 * <li> <tt>obj_grad_info.Gf != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*obj_grad_info.Gf</tt> is updated to \a Gf(x).
	 * </ul>
	 *
	 * @param x       [in]  Unknown vector (size n).
	 * @param newx    [in]  True if is a new point.
	 * @param obj_grad_info
	 *                [out] Pointers to \c f, \c c, \c h and \c Gf.
	 *                On output <tt>*obj_grad_info.Gf</tt> is updated to \a Gf(x).
	 *                Any of the other objects pointed to in
	 *                \c obj_grad_info may be set if <tt>this->mult_calc() == true</tt> but are
	 *                now guaranteed to be.
	 */
	virtual void imp_calc_Gf(const VectorWithOp& x, bool newx, const ObjGradInfo& obj_grad_info) const = 0;

	//@}

private:
	mutable VectorWithOpMutable     *Gf_;
	mutable size_type				num_Gf_evals_;

};	// end class NLPObjGradient

// //////////////////
// Inline members

inline
const NLPObjGradient::ObjGradInfo NLPObjGradient::obj_grad_info() const
{
	return ObjGradInfo(Gf_,zero_order_info());
}

}	// end namespace NLPInterfacePack 

#endif // NLP_OBJ_GRADIENT_H