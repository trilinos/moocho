// //////////////////////////////////////////
// ExampleNLPFirstOrderDirect.hpp
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

#ifndef EXAMPLE_NLP_FIRST_ORDER_DIRECT_H
#define EXAMPLE_NLP_FIRST_ORDER_DIRECT_H

#include "ExampleNLPObjGradient.hpp"
#include "NLPInterfacePack/src/NLPFirstOrderDirect.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpMutable.hpp"
#include "AbstractLinAlgPack/src/VectorSpace.hpp"
#include "AbstractLinAlgPack/src/VectorSpaceCompositeStd.hpp"

namespace NLPInterfacePack {

///
/** Simple example %NLP subclass to illustrate how to implement the
 * \c NLPFirstOrderDirect interface for a specialized \c NLP.
 *
 * For the NLP formulation see <tt>ExampleNLPObjGradient</tt>.
 *
 * In this %NLP we will select the first <tt>n/2</tt> variables as the dependent
 * or basis variables.  Note that the variable bounds are on the dependent
 * or independent variables variables if has_bounds = true.  Otherwise
 * there are no variable bounds.  Also the starting point can be
 * varied.
 *
 * The implementation of this %NLP subclass is actually independent from the vector
 * space for the dependent and independent variables (same vector space).
 * This implementation is defined entirely based on an arbitrary
 * <tt>VectorSpace</tt> object that is passed to the constructor
 * \c ExampleNLPFirstOrderDirect().  This %NLP subclass uses a
 * <tt>\ref AbstractLinAlgPack::VectorSpaceCompositeStd "VectorSpaceCompositeStd"</tt>
 * object to represent the space for <tt>[ x_dep; x_indep ]</tt>
 *
 * The quantities computed by this subclass include:
 \verbatim

  py = -inv(C)*c
   
  D = -inv(C)*N
  
  where:
  
    Gc' = [ C , N ]
  
  		[ x(m+1) - 1									]
  		[				x(m+2) - 1						]
  	C = [							.					]
  		[								.				]
  		[									x(m+m) - 1	]
  
  
  		[ x(1) - 10										]
  		[				x(2) - 10						]
  	N = [							.					]
  		[								.				]
  		[									x(m) - 10	]
  
 \endverbatim 
 * Here \c Gc is never computed explicitly.
 *
 * To make this possible this subclass relies on some specialized RTOp operators which
 * are implemented in C (for portability).  These operator classes are declared in the header
 * file <tt>ExampleNLPFirstOrderDirect.hpp</tt> and are documented \ref explnlp2_ops_grp "here".
 */
class ExampleNLPFirstOrderDirect
	: virtual public NLPFirstOrderDirect
	, virtual public ExampleNLPObjGradient
{
public:

	///
	/** Constructor.
	 *
	 * @param  vec_space  [in] Smart pointer to a vector space object that will
	 *                    be used to define the spaces of dependent and independent
	 *                    variables.
	 * @param  xo         [in] The initial starting guess for \a x.
	 * @param  has_bounds [in] If \c true, then the NLP will have bounds.  If \c false
	 *                    then it will not have bounds.
	 * @param  dep_bouned [in] If \c true, then the bounds will be on the dependent
	 *                    variables.  If \c false, then the bounds will be on the
	 *                    independent variable.  This argument is ignored if
	 *                    <tt>has_bounds == false</tt>.
	 */
	ExampleNLPFirstOrderDirect(
		const VectorSpace::space_ptr_t&  vec_space
		,value_type                      xo
		,bool                            has_bounds
		,bool                            dep_bounded
		);

	/** @name Overridden public members from NLP */
	//@{

	///
	void initialize(bool test_setup);
	///
	bool is_initialized() const;
	/// Returns <tt>return.get() == NULL</tt>.
	vec_space_ptr_t space_h() const;
	/// Throws exception.
	const VectorWithOp& hl() const;
	/// Throws exception.
	const VectorWithOp& hu() const;

	//@}

	/** @name Overridden public members from NLPFirstOrderDirect */
	//@{

	///
	Range1D var_dep() const;
	///
	Range1D var_indep() const;
	///
	const mat_fcty_ptr_t factory_D() const;
	///
	void calc_point(
		const VectorWithOp     &x
		,value_type            *f
		,VectorWithOpMutable   *c
		,bool                  recalc_c
		,VectorWithOpMutable   *h
		,VectorWithOpMutable   *Gf
		,VectorWithOpMutable   *py
		,VectorWithOpMutable   *rGf
		,MatrixWithOp          *GcU
		,MatrixWithOp          *Gh
		,MatrixWithOp          *D
		,MatrixWithOp          *V
		,MatrixWithOp          *P
		) const;
	///
	void calc_semi_newton_step(
		const VectorWithOp    &x
		,VectorWithOpMutable  *c
		,bool                 recalc_c
		,VectorWithOpMutable  *py
		) const;

	//@}

protected:

	/** @name Overridden protected members from NLP */
	//@{

	/// This implementation does nothing (should never be called though).
	void imp_calc_h(const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const;

	//@}

private:

	// /////////////////////////////////////////
	// Private data members

	mat_fcty_ptr_t              factory_D_;         // Matrix space object for D

	bool         initialized_;            // flag for if initialized has been called.

	// /////////////////////////////////////////
	// Private member functions

	///
	void assert_is_initialized() const;

};	// end class ExampleNLPFirstOrderDirect

// ///////////////////////////////////////////////
// Inline member functions

inline
void ExampleNLPFirstOrderDirect::assert_is_initialized() const
{
    using NLPInterfacePack::NLP;
	if( !is_initialized() )
		throw NLP::UnInitialized("ExampleNLPFirstOrderDirect::assert_is_initialized() : Error, "
			"ExampleNLPFirstOrderDirect::initialize() has not been called yet." );
}

}	// end namespace NLPInterfacePack

#endif	// EXAMPLE_NLP_FIRST_ORDER_DIRECT_H
