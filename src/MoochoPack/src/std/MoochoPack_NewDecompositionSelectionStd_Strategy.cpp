// ////////////////////////////////////////////////////////////////////////////
// NewDecompositionSelectionStd_Strategy.cpp
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

#include "ReducedSpaceSQPPack/include/std/NewDecompositionSelectionStd_Strategy.h"
#include "ReducedSpaceSQPPack/include/std/rSQPAlgorithmStepNames.h"
#include "ReducedSpaceSQPPack/include/rSQPAlgo.h"
#include "ReducedSpaceSQPPack/include/rSQPState.h"

namespace ReducedSpaceSQPPack {

NewDecompositionSelectionStd_Strategy::NewDecompositionSelectionStd_Strategy(
	const decomp_sys_handler_ptr_t   &decomp_sys_handler
	)
	:decomp_sys_handler_(decomp_sys_handler)
{}

bool NewDecompositionSelectionStd_Strategy::new_decomposition(
	rSQPAlgo& algo, Algorithm::poss_type step_poss
	,GeneralIterationPack::EDoStepType type, Algorithm::poss_type assoc_step_poss
	)
{
	rSQPState               &s     = algo.rsqp_state();
	EJournalOutputLevel     olevel = algo.algo_cntr().journal_output_level();
	std::ostream&           out    = algo.track().journal_out();

	// Check to see if we have a decomposition system set
	if( !get_decomp_sys_handler().get() ) {
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
			out	<< "\nWe are asked to select a new basis but there is no\n"
					"decomposition system set so we have no choice but to terminiate\n"
					"the algorithm"
				<< " (k = " << algo.state().k() << ")\n";
		}
		algo.terminate(false);
		return false;
	}

	// We may get an infinite loop here so make sure we are under the max
	// number of iterations.
	if( s.k() >= algo.algo_cntr().max_iter() ) {
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
			out	<< "\nThe maximum number of rSQP iterations\n"
				<< " have been exceeded so quit "
				<< " (k = " << algo.state().k() << ")\n";
		}
		algo.terminate(false);
		return false;
	}

	// Select a new decomposition
	decomp_sys_handler().select_new_decomposition(true);
	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		out	<< "x_kp1 = x_k\n"
			<< "k=k+1\n"
			<< "goto EvalNewPoint\n";
	}
	s.x().set_k(1) = s.x().get_k(0);
	s.alpha().set_k(0) = 0.0;	// Show that no step was taken.
	algo.track().output_iteration( algo );
	s.next_iteration();
	algo.do_step_next( EvalNewPoint_name );
	return false;
}

void NewDecompositionSelectionStd_Strategy::print_new_decomposition(
	const rSQPAlgo& algo, Algorithm::poss_type step_poss
	,GeneralIterationPack::EDoStepType type, Algorithm::poss_type assoc_step_poss
	,std::ostream& out, const std::string& L
	) const
{
	out
		<< L << "if k > max_iter then\n"
		<< L << "  terminate the algorithm\n"
		<< L << "end\n"		
		<< L << "Select a new basis at current point\n"
		<< L << "x_kp1 = x_k\n"
		<< L << "alpha_k = 0\n"
		<< L << "k=k+1\n"
		<< L << "goto EvalNewPoint\n";
}

}	// end namespace ReducedSpaceSQPPack 