// ////////////////////////////////////////////////////////////////////////////
// PreEvalNewPointBarrier_Step.cpp
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

#include <ostream>
#include <typeinfo>
#include <iostream>

#include "AbstractLinAlgPack/src/abstract/tools/MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorStdOps.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorOut.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/assert_print_nan_inf.hpp"
#include "NLPInterfacePack/src/abstract/interfaces/NLPFirstOrder.hpp"
#include "IterationPack/src/print_algorithm_step.hpp"
#include "MoochoPack/src/IpState.hpp"
#include "MoochoPack/src/std/PreEvalNewPointBarrier_Step.hpp"
#include "MoochoPack/src/moocho_algo_conversion.hpp"

#include "StringToBool.hpp"

#include "dynamic_cast_verbose.hpp"

namespace MoochoPack {

PreEvalNewPointBarrier_Step::PreEvalNewPointBarrier_Step(
  const value_type relative_bound_push,
  const value_type absolute_bound_push
  )
	:
	relative_bound_push_(relative_bound_push),
	absolute_bound_push_(absolute_bound_push)
	{}

bool PreEvalNewPointBarrier_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
	{
	using DynamicCastHelperPack::dyn_cast;
	using IterationPack::print_algorithm_step;
	using AbstractLinAlgPack::force_in_bounds_buffer;

	NLPAlgo            &algo   = dyn_cast<NLPAlgo>(_algo);
	IpState             &s      = dyn_cast<IpState>(_algo.state());
	NLP                 &nlp    = algo.nlp();
	NLPFirstOrder   *nlp_foi = dynamic_cast<NLPFirstOrder*>(&nlp);
	
	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();
	
	if(!nlp.is_initialized())
		nlp.initialize(algo.algo_cntr().check_results());

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
		{
		using IterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
		}

	IterQuantityAccess<value_type>           &barrier_parameter_iq = s.barrier_parameter();
	IterQuantityAccess<VectorMutable>  &x_iq  = s.x();
	IterQuantityAccess<MatrixSymDiagStd> &Vl_iq = s.Vl();
	IterQuantityAccess<MatrixSymDiagStd> &Vu_iq = s.Vu();

	if( x_iq.last_updated() == IterQuantity::NONE_UPDATED ) 
		{
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
			{
			out << "\nInitialize x with x_k = nlp.xinit() ...\n"
				<< " and push x_k within bounds.\n";
			}
		VectorMutable& x_k = x_iq.set_k(0) = nlp.xinit();
  
		// apply transformation operator to push x sufficiently within bounds
		force_in_bounds_buffer(relative_bound_push_, 
							   absolute_bound_push_,
							   nlp.xl(),
							   nlp.xu(),
							   &x_k);

		// evaluate the func and constraints
		IterQuantityAccess<value_type>
			&f_iq    = s.f();
		IterQuantityAccess<VectorMutable>
			&Gf_iq   = s.Gf(),
			*c_iq    = nlp.m() > 0 ? &s.c() : NULL;
		IterQuantityAccess<MatrixOp>
			*Gc_iq   = nlp_foi ? &s.Gc() : NULL;

		using AbstractLinAlgPack::assert_print_nan_inf;
		assert_print_nan_inf(x_k, "x", true, NULL); // With throw exception if Inf or NaN!

		// Wipe out storage for computed iteration quantities (just in case?) : RAB: 7/29/2002
		if(f_iq.updated_k(0))
			f_iq.set_not_updated_k(0);
		if(Gf_iq.updated_k(0))
			Gf_iq.set_not_updated_k(0);
		if (c_iq)
			{
			if(c_iq->updated_k(0))
				c_iq->set_not_updated_k(0);
			}
		if (nlp_foi)
			{
			if(Gc_iq->updated_k(0))
				Gc_iq->set_not_updated_k(0);
			}
		}

	if (barrier_parameter_iq.last_updated() == IterQuantity::NONE_UPDATED)
		{
		barrier_parameter_iq.set_k(-1) = 0.1; // RAB: 7/29/2002: This should be parameterized! (allow user to set this!)
		}

	// Print vector information
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
		{
		out << "x_k =\n" << x_iq.get_k(0);
		}

	return true;
	}


void PreEvalNewPointBarrier_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
	{
	//const NLPAlgo   &algo = rsqp_algo(_algo);
	//const NLPAlgoState  &s    = algo.rsqp_state();
	out << L << "# Evaluate information specific to primal / dual barrier algorithms\n"
		<< L << "if (x never updated) then\n"
		<< L << "  x_k = nlp.xinit()\n"
		<< L << "  force_in_bounds(x_k)\n"
		<< L << "  set f_k not updated\n"
		<< L << "  set Gf_k not updated\n"
		<< L << "  if (m > 0) then\n"
		<< L << "    set c_k not updated\n"
		<< L << "    set Gc_k not updated\n"
		<< L << "end\n";
	}


namespace {

const int local_num_options = 2;

enum local_EOptions 
	{
	RELATIVE_BOUND_PUSH,
	ABSOLUTE_BOUND_PUSH
	};

const char* local_SOptions[local_num_options] = 
	{
	"relative_bound_push",
	"absolute_bound_push"
	};

}

 
PreEvalNewPointBarrier_StepSetOptions::PreEvalNewPointBarrier_StepSetOptions(
  PreEvalNewPointBarrier_Step* target
  , const char opt_grp_name[] )
	:
	OptionsFromStreamPack::SetOptionsFromStreamNode(
	  opt_grp_name, local_num_options, local_SOptions ),
	OptionsFromStreamPack::SetOptionsToTargetBase< PreEvalNewPointBarrier_Step >( target )
	{
	}

void PreEvalNewPointBarrier_StepSetOptions::set_option( 
  int option_num, const std::string& option_value )
	{
	using OptionsFromStreamPack::StringToBool;

	typedef PreEvalNewPointBarrier_Step target_t;
	switch( (local_EOptions)option_num ) 
		{
		case RELATIVE_BOUND_PUSH:
			target().relative_bound_push(::atof(option_value.c_str()));
			break;
		case ABSOLUTE_BOUND_PUSH:
			target().absolute_bound_push(::atof(option_value.c_str()));
			break;
		default:
			assert(0);	// Local error only?
		}
	}

}  // end namespace MoochoPack
