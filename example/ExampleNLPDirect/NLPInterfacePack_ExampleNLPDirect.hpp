// //////////////////////////////////////////
// ExampleNLPFirstOrderDirect.h
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

#include "NLPInterfacePack/include/NLPFirstOrderDirect.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/VectorSpaceCompositeStd.h"

namespace NLPInterfacePack {

///
/** Simple example %NLP subclass to illustrate how to implement the
 * \c NLPFirstOrderDirect interface for a specialized \c NLP.
 *
 * The example %NLP we will use is a scalable problem where
 * the basis of the jacobian of the constraints is a diagonal
 * matrix.
 \verbatim

    min    f(x) = (1/2) * sum( x(i)^2, for i = 1..n )
    s.t.   c(x)(j) = x(j) * (x(m+j) -1) - 10 * x(m+j) = 0, for j = 1..m
          0.01 < x(i) < 20, for i = p...p+m

    where:
        m = n/2
        p = 1 if dep_bounded == true or m+1 if dep_bounded = false
 \endverbatim
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
 * file <tt>ExampleNLPFirstOrderDirect.h</tt> and are documented \ref explnlp2_ops_grp "here".
 */
class ExampleNLPFirstOrderDirect
	: public NLPFirstOrderDirect
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
	void initialize();
	///
	bool is_initialized() const;
	///
	size_type n() const;
	///
	size_type m() const;
	///
	vec_space_ptr_t space_x() const;
	///
	vec_space_ptr_t space_c() const;
	///
    size_type num_bounded_x() const;
	///
	void force_xinit_in_bounds(bool force_xinit_in_bounds);
	///
	bool force_xinit_in_bounds() const;
	///
	const VectorWithOp& xinit() const;
	///
	const VectorWithOp& xl() const;
	///
	const VectorWithOp& xu() const;
	///
	void scale_f( value_type scale_f );
	///
	value_type scale_f() const;
	///
	void report_final_solution(
		const VectorWithOp&    x
		,const VectorWithOp*   lambda
		,const VectorWithOp*   lambdaI
		,const VectorWithOp*   nu
		,bool                  optimal
		) const;

	//@}

	/** @name Overridden public members from NLPFirstOrderDirect */
	//@{

	///
	const mat_space_ptr_t& space_D() const;
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

	///
	void imp_calc_f(
		const VectorWithOp& x, bool newx
		,const ZeroOrderInfo& zero_order_info) const;
	///
	void imp_calc_c(
		const VectorWithOp& x, bool newx
		,const ZeroOrderInfo& zero_order_info) const;

	//@}

	/** @name Overridden protected members from NLPObjGradient */
	//@{

	///
	void imp_calc_Gf(
		const VectorWithOp& x, bool newx
		,const ObjGradInfo& obj_grad_info) const;

	//@}

private:

	// /////////////////////////////////////////
	// Private types

	typedef ReferenceCountingPack::ref_count_ptr<VectorSpaceCompositeStd>  vec_space_comp_t;

	// /////////////////////////////////////////
	// Private data members

	VectorSpace::space_ptr_t    vec_space_;       // The vector space for dependent and indepenent variables and c(x).
	vec_space_comp_t            vec_space_comp_;  // Composite vector space for x = [ xD; xI ]
	mat_space_ptr_t             space_D_;         // Matrix space object for D

	bool         initialized_;            // flag for if initialized has been called.
	value_type   obj_scale_;              // default = 1.0;
	bool         has_bounds_;             // default = true
	bool         force_xinit_in_bounds_;  // default = true.

	size_type	n_;                       // Number of variables in the problem.
	VectorSpace::vec_mut_ptr_t  xinit_;   // Initial guess.
	VectorSpace::vec_mut_ptr_t  xl_;      // lower bounds.
	VectorSpace::vec_mut_ptr_t  xu_;      // upper bounds.

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