// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef NLPIP_NLP_DIRECT_THYRA_MODEL_EVALUATOR_HPP
#define NLPIP_NLP_DIRECT_THYRA_MODEL_EVALUATOR_HPP

#include "NLPInterfacePack_NLPThyraModelEvaluator.hpp"
#include "NLPInterfacePack_NLPDirect.hpp"

namespace NLPInterfacePack {

/** \brief Implement the %NLPFirstOrder interface using a
 * <tt>Thyra::ModelEvaluator</tt> object.
 *
 * ToDo: Finish documentation!
 */
class NLPDirectThyraModelEvaluator
  : virtual public NLPDirect
  , virtual public NLPThyraModelEvaluator
{
public:

	/** \brief Initialize to uninitialized */
	NLPDirectThyraModelEvaluator();

	/** \brief Calls <tt>initialize()</tt>. */
	NLPDirectThyraModelEvaluator(
		const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &model
    ,const int                                                      p_idx
    ,const int                                                      g_idx
		,const Thyra::VectorBase<value_type>                            *model_xL      = NULL
		,const Thyra::VectorBase<value_type>                            *model_xU      = NULL
		,const Thyra::VectorBase<value_type>                            *model_x0      = NULL
		,const Thyra::VectorBase<value_type>                            *model_pL      = NULL
		,const Thyra::VectorBase<value_type>                            *model_pU      = NULL
		,const Thyra::VectorBase<value_type>                            *model_p0      = NULL
		);

	/** \brief .Initialize given a <tt>Thyra::ModelEvaluator</tt> and
	 * a description of how to interpret it.
	 *
	 * ToDo: Finish documentation!
	 */
	void initialize(
		const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &model
    ,const int                                                      p_idx
    ,const int                                                      g_idx
		,const Thyra::VectorBase<value_type>                            *model_xL      = NULL
		,const Thyra::VectorBase<value_type>                            *model_xU      = NULL
		,const Thyra::VectorBase<value_type>                            *model_x0      = NULL
		,const Thyra::VectorBase<value_type>                            *model_pL      = NULL
		,const Thyra::VectorBase<value_type>                            *model_pU      = NULL
		,const Thyra::VectorBase<value_type>                            *model_p0      = NULL
		);

	/** @name Overridden public members from NLP */
	//@{

	/** \brief . */
	void initialize(bool test_setup);
	/** \brief . */
	void unset_quantities();

	//@}

	/** @name Overridden public members from NLPDirect */
	//@{

	/** \brief . */
	Range1D var_dep() const;
	/** \brief . */
	Range1D var_indep() const;
	/** \brief . */
	const mat_fcty_ptr_t factory_D() const;
  /** \brief . */
	const mat_sym_nonsing_fcty_ptr_t factory_S() const;
	/** \brief . */
	void calc_point(
		const Vector     &x
		,value_type      *f
		,VectorMutable   *c
		,bool            recalc_c
		,VectorMutable   *Gf
		,VectorMutable   *py
		,VectorMutable   *rGf
		,MatrixOp        *GcU
		,MatrixOp        *D
		,MatrixOp        *Uz
		) const;
	/** \brief . */
	void calc_semi_newton_step(
		const Vector    &x
		,VectorMutable  *c
		,bool           recalc_c
		,VectorMutable  *py
		) const;

	//@}

private:

  mutable Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<value_type> >  thyra_C_;
  mutable Teuchos::RefCountPtr<Thyra::MultiVectorBase<value_type> >        thyra_N_;
  
};	// end class NLPDirectThyraModelEvaluator

}	// end namespace NLPInterfacePack

#endif	// NLPIP_NLP_DIRECT_THYRA_MODEL_EVALUATOR_HPP