// //////////////////////////////////////////
// NLPBarrier.cpp
//
// Copyright (C) 2001
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

#include <math.h>
#include <iostream>
#include <limits>

#include "NLPInterfacePack/src/abstract/tools/NLPBarrier.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorSpace.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorOut.hpp"
#include "ThrowException.hpp"

namespace NLPInterfacePack {

NLPBarrier::NLPBarrier()
	:
	barrier_term_(0.0),
	objective_term_(0.0),
	nlp_(MemMngPack::null)
	{
	}


void NLPBarrier::InitializeFromNLP(
  MemMngPack::ref_count_ptr<NLP> original_nlp
  )
	{
	THROW_EXCEPTION(
	  !original_nlp.get(),
	  std::logic_error,
	  "null nlp passed to NLPBarrier decorator"
	  );

	nlp_ = MemMngPack::rcp_dynamic_cast<NLPObjGrad>(original_nlp);

	THROW_EXCEPTION(
	  !nlp_.get(),
	  std::logic_error,
	  "non NLPObjGrad NLP passed to NLPBarrier decorator"
	  );
	}

void NLPBarrier::mu(const value_type mu)
	{
	mu_ = mu;
	}

value_type NLPBarrier::barrier_term() const
	{
	return barrier_term_;
	}

value_type NLPBarrier::objective_term() const
	{
	return objective_term_;
	}

const MemMngPack::ref_count_ptr<Vector> NLPBarrier::grad_barrier_term() const
	{
	return grad_barrier_term_;
	}

const MemMngPack::ref_count_ptr<Vector>  NLPBarrier::grad_objective_term() const
	{
	return grad_objective_term_;
	}


void NLPBarrier::calc_f(const Vector& x, bool newx) const
	{
	nlp_->calc_f(x, newx);
	value_type* f_val = nlp_->get_f();

	objective_term_ = *f_val;
	barrier_term_   = CalculateBarrierTerm(x);

	(*f_val) += barrier_term_;
	}

void NLPBarrier::calc_Gf(const Vector& x, bool newx) const
	{
	using AbstractLinAlgPack::inv_of_difference;

   	nlp_->calc_Gf(x, newx);
	grad_objective_term_ = nlp_->get_Gf()->clone();

	//std::cout << "grad_objective_term=\n";
	//grad_objective_term_->output(std::cout);

	if (!grad_barrier_term_temp_.get())
		{ grad_barrier_term_temp_ = grad_objective_term_->space().create_member(); }

	if (!grad_barrier_term_.get())
		{ grad_barrier_term_ = grad_objective_term_->space().create_member(); }	

	*grad_barrier_term_temp_ = 0.0;
	*grad_barrier_term_ = 0.0;	

	inv_of_difference(mu_, nlp_->xu(), x, grad_barrier_term_.get());
 	//std::cout << "mu*invXU=\n";
	//grad_barrier_term_->output(std::cout);

	inv_of_difference(mu_, x, nlp_->xl(), grad_barrier_term_temp_.get());
 	//std::cout << "mu*invXL=\n";
	//grad_barrier_term_temp_->output(std::cout);

	grad_barrier_term_->axpy(-1.0, *grad_barrier_term_temp_);

	nlp_->get_Gf()->axpy(1.0, *grad_barrier_term_);

	//std::cout << "grad_objective_term with barrier=\n";
	//nlp_->get_Gf()->output(std::cout);
	}

void NLPBarrier::imp_calc_f(
  const Vector& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info
  ) const
	{
	assert(false && !"This should never get called.");
	}

void NLPBarrier::imp_calc_c(
  const Vector& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info
  ) const
	{
	assert(false && !"This should never get called.");
	}

void NLPBarrier::imp_calc_c_breve(
  const Vector& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info_breve
  ) const
	{	
	assert(false && !"This should never get called.");
	}

void NLPBarrier::imp_calc_h_breve(
  const Vector& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info_breve
  ) const
	{	
	assert(false && !"This should never get called.");
	}

void NLPBarrier::imp_calc_Gf(
  const Vector& x,
  bool newx, 
  const ObjGradInfo& obj_grad_info
  ) const
	{
	assert(false && !"This should never get called.");
	}


value_type NLPBarrier::CalculateBarrierTerm(const Vector& x) const
	{
	using AbstractLinAlgPack::log_bound_barrier;
	barrier_term_ = log_bound_barrier(x, xl(), xu());
//	std::cerr << "NLPBarrier::CalculateBarrierTerm(x) : (1) barrier_term_ = " << barrier_term_ << std::endl;
	barrier_term_ *= -mu_;
//	std::cerr << "NLPBarrier::CalculateBarrierTerm(x) : (2) barrier_term_ = " << barrier_term_ << std::endl;
	return barrier_term_;
	}

}	// end namespace NLPInterfacePack
