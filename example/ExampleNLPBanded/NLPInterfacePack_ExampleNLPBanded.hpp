// //////////////////////////////////////////
// ExampleNLPBanded.h
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

#ifndef EXAMPLE_NLP_BANDED_H
#define EXAMPLE_NLP_BANDED_H

#include "NLPInterfacePack/include/NLPSerialPreprocessExplJac.h"

namespace NLPInterfacePack {

///
/** Simple scalable serial %NLP subclass.
 *
 * This example %NLP is a scalable problem where the basis of the jacobian of the
 * equality constraints is a banded (band width = bw) symmetric positive definite
 * matrix.  Both the number of dependnet and independent variables can be varied.
 *
 * To setup this %NLP, the client specifies:<ul>
 * <li> \c nD : the number of dependent variables
 * <li> \c nI : the number of independent variables (<tt>nI <= nD</tt>)
 * <li> \c bw : the band width (defined as the number of diagonals, including the
 *              center diagonal i==j that have non-zero element)
 * <li> \c mU : the number of undecomposed dependent constraints
 * <li> \c mI : the number of general constraints
 * <li> \c diag_scal : Constant scaling factor for the diagonal elements
 * <li> \c diag_vary : The scaling factor for the diagonal elements (to produce
 *                     illconditioning) between the largets and the smallest.
 * <li> \c sym_basis : True if the basis (selected by the NLP) is symmetric or not.
 * </ul>
 *
 * This %NLP is defined as:
 \verbatim

    min    f(x) = (1/2) * sum( x(i)^2, for i = 1..n )
    s.t.
           c(j) = ( ds(j)*x(j)                                 \
                    - sum( 3/(k)*x(j-k), k=1...klu(j) )        |
                    - sum( fu/(k)*x(j+k), k=1...kuu(j) )       | for j  = 1...nD
                   ) * (x(nD+q(j)) + 1)^2                      |
                   + co(j) == 0                                /

            c(nD+jU) = c(jU) + co(nD+jU)  == 0                 } for jU = 1...mU

            hl(jI) <= x(jI) - x(nD+q(jI)) <= hu(jI)            } for jI = 1...mI

            xl(i) <= x(i) <= xu(i)                             } for i  = 1...n

    where:

        n = nD + nI

        m = nD + mU

        mI = mI

        ds(j) = diag_scal * ( (diag_vary - 1)/(nD -1) * (j - 1) + 1 )

             / 3  : if sym_basis = true
        fu = |
             \ 6  : if sym_basis = false

                                 / 2  : if floor((j-1)/nI) < nD % nI
        q(j) = floor((j-1)/nI) + |
                                 \ 1  : if floor((j-1)/nI) >= nD % nI

                                                                  
                  / bw-1 : if j - bw >= 0                         \
        klu(j) =  |                                               |
                  \ j-1  : if j - bw <= 1                         |
                                                                  | for j=1...nD
                  / bw-1 : if j + bw-1 <= nD                      |
        kuu(j) =  |                                               |
                  \ nD-j : if j - bw <= 1                         /
 \endverbatim
 * In the above formuation, the sums are not computed if the upper bounds on \c k
 * are zero.  The term <tt>co(j)</tt> is an adjustable term that can be used
 * to manipulate the solution.  Note that if <tt>co(nD+jI) != 0</tt> above, then
 * the undecomposed dependent equality constraints are inconsistent with the
 * decomposed equalities and therefore the NLP is infeasible.  An infeasible NLP
 * can also be created by manipulating \c xl(i), \c xu(i), \c hl(jI), \c hu(jI)
 * and \c co(j).
 *
 * For the above NLP, the Jacobian of the decomposed equalities has Jacobian elements:
 \verbatim
    
                      /  -3/(j-i) * (x(nD+q(j)) + 1)^2         : i - klu(i) <= j < i
                      |
                      |  ds(j) * (x(nD+q(j)) + 1)^2            : i == j
   d(c(j))/d(x(i)) =  |
                      |  -3/(i-j) * (x(nD+q(j)) + 1)^2         : i < j <= i + kuu(i)
                      |
                      |  2 * (c(j) - co(j)) / (x(nD+q(j)) + 1) : i == nD + q
                      | 
                      \  0                                     : otherwise
                      
                      , for j = 1...nD, i = 1...nD+nI
 \endverbatim
 * The above definition shows that for the independent variables, the Jacobian
 * elements are written in terms of the constraint \c c(j).  This fact
 * is exploited in the computational routines when <tt>this->multi_calc() == true</tt>.
 *
 * For <tt>nD == 7, nI == 2, bw = 2</tt> with <tt>floor(nD/nI) = 3</tt> and
 * <tt>nD \% nI = 1</tt>, the Jacobian <tt>Gc'</tt> looks like:
 \verbatim

   1 | x  x                 x    |
   2 | x  x  x              x    |
   3 |    x  x  x           x    |
   4 |       x  x  x        x    |
   5 |          x  x  x        x |
   6 |             x  x  x     x |
   7 |                x  x     x |
       -  -  -  -  -  -  -  -  -
       1  2  3  4  5  6  7  8  9
 \endverbatim
 * 
 * ToDo: Finish documentation!
 */
class ExampleNLPBanded
	: public NLPSerialPreprocessExplJac
{
public:

	/** @name Constructors / initializers */
	//@{

	///
	/** Constructor.
	 *
	 * ToDo: Finish documentation!
	 */
	ExampleNLPBanded(
		size_type     nD
		,size_type    nI
		,size_type    bw                = 1
		,size_type    mU                = 0
		,size_type    mI                = 0
		,value_type   xo                = 0.1
		,value_type   xDl               = -NLP::infinite_bound()
		,value_type   xDu               = +NLP::infinite_bound()
		,value_type   xIl               = -NLP::infinite_bound()
		,value_type   xIu               = +NLP::infinite_bound()
		,value_type   hl                = -NLP::infinite_bound()
		,value_type   hu                = +NLP::infinite_bound()
		,bool         nlp_selects_basis = false
		,value_type   diag_scal         = 10.0
		,value_type   diag_vary         = 1.0
		,bool         sym_basis         = false
		);

	//@}

	// Todo: Add methods to manipulate bounds and co ...

	/** @name Access */
	//@{
	
	// ToDo: Add these ...

	//@}

	/** @name Overridden public members from NLP */
	//@{

	///
	void initialize(bool test_setup);
	///
	bool is_initialized() const;
	///
	value_type max_var_bounds_viol() const;
	///
	void set_multi_calc(bool multi_calc) const;
	///
	bool multi_calc() const;

	//@}

	/** @name Overridden from NLPVarReductPerm */
	//@{
	
	///
	bool nlp_selects_basis() const;

	//@}

protected:

	/** @name Overridden protected methods from NLPSerialPreprocess */
	//@{

	///
	bool imp_nlp_has_changed() const;
	///
	size_type imp_n_orig() const;
	///
	size_type imp_m_orig() const;
	///
	size_type imp_mI_orig() const;
	///
	const VectorSlice imp_xinit_orig() const;
	///
	bool imp_has_var_bounds() const;
	///
	const VectorSlice imp_xl_orig() const;
	///
	const VectorSlice imp_xu_orig() const;
	///
	const VectorSlice imp_hl_orig() const;
	///
	const VectorSlice imp_hu_orig() const;
	///
	void imp_calc_f_orig(
		const VectorSlice            &x_full
		,bool                        newx
		,const ZeroOrderInfoSerial   &zero_order_info
		) const;
	///
	void imp_calc_c_orig(
		const VectorSlice            &x_full
		,bool                        newx
		,const ZeroOrderInfoSerial   &zero_order_info
		) const;
	///
	void imp_calc_h_orig(
		const VectorSlice            &x_full
		,bool                        newx
		,const ZeroOrderInfoSerial   &zero_order_info
		) const;
	///
	void imp_calc_Gf_orig(
		const VectorSlice            &x_full
		,bool                        newx
		,const ObjGradInfoSerial     &obj_grad_info
		) const;
	///
	bool imp_get_next_basis(
		IVector      *var_perm
		,IVector     *equ_perm
		,size_type   *rank
		);
	///
	void imp_report_orig_final_solution(
		const VectorSlice      &x_orig
		,const VectorSlice     *lambda_orig
		,const VectorSlice     *lambdaI_orig
		,const VectorSlice     *nu_orig
		,bool                  is_optimal
		) const;

	//@}
	
	/** @name Overridden protected methods from NLPSerialPreprocessExplJac */
	//@{

	///
	size_type imp_Gc_nz_orig() const;
	///
	size_type imp_Gh_nz_orig() const;
	///
	void imp_calc_Gc_orig(
		const VectorSlice& x_full, bool newx
		, const FirstOrderExplInfo& first_order_expl_info
		) const;
	///
	void imp_calc_Gh_orig(
		const VectorSlice& x_full, bool newx
		, const FirstOrderExplInfo& first_order_expl_info
		) const;

	//@}

private:

	// /////////////////////////////////////////
	// Private types

	// /////////////////////////////////////////
	// Private data members

	bool         is_initialized_;

	bool         nlp_selects_basis_;
	bool         basis_selection_was_given_;

	bool         has_var_bounds_;

	mutable bool multi_calc_;

	size_type    nD_;
	size_type    nI_;
	size_type    bw_;
	size_type    mU_;
	size_type    mI_;

	size_type    Gc_orig_nz_;
	size_type    Gh_orig_nz_;

	Vector       xinit_orig_;
	Vector       xl_orig_;
	Vector       xu_orig_;
	Vector       hl_orig_;
	Vector       hu_orig_;

	Vector       co_orig_;

	mutable bool  c_orig_updated_;

	value_type   diag_scal_;
	value_type   diag_vary_;
	value_type   fu_;

	// /////////////////////////////////////////
	// Private member functions

	///
	void assert_is_initialized() const;

	///
	void inform_new_point(bool newx) const;

	// Not defined and not to be called
	ExampleNLPBanded();
	ExampleNLPBanded(const ExampleNLPBanded&);
	ExampleNLPBanded& operator=(const ExampleNLPBanded&);

};	// end class ExampleNLPBanded

}	// end namespace NLPInterfacePack

#endif	// EXAMPLE_NLP_BANDED_H
