// ////////////////////////////////////////////////////////////////////////////
// LineSearchFailureNewDecompositionSelection_Step.cpp
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

#include <ostream>
#include <typeinfo>

#include "ReducedSpaceSQPPack/src/std/LineSearchFailureNewDecompositionSelection_Step.hpp"
#include "ReducedSpaceSQPPack/src/std/rSQPAlgorithmStepNames.hpp"
#include "ReducedSpaceSQPPack/src/rsqp_algo_conversion.hpp"
#include "ReducedSpaceSQPPack/src/ReducedSpaceSQPPackExceptions.hpp"
#include "IterationPack/src/print_algorithm_step.hpp"

namespace ReducedSpaceSQPPack {

LineSearchFailureNewDecompositionSelection_Step::LineSearchFailureNewDecompositionSelection_Step(
	const line_search_step_ptr_t        &line_search_step
	,const new_decomp_strategy_ptr_t    &new_decomp_strategy
	)
	:line_search_step_(line_search_step)
	,new_decomp_strategy_(new_decomp_strategy)
	,last_ls_failure_k_(-100) // has not failed yet
{}

bool LineSearchFailureNewDecompositionSelection_Step::do_step(
	Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
	)
{
	try {
		return line_search_step().do_step(_algo,step_poss,type,assoc_step_poss);
	}
	catch(const LineSearchFailure& lsf_excpt) {
		rSQPAlgo	&algo	= rsqp_algo(_algo);
		rSQPState	&s		= algo.rsqp_state();
		NLP			&nlp	= algo.nlp();

		EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
		std::ostream& out = algo.track().journal_out();

		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
			out	<< "\nLine search failed "
				<< " (k = " << algo.state().k() << ")\n"
				<< "LineSearchFailure description: " << lsf_excpt.what() << "\n";
		}

		if( last_ls_failure_k_ == s.k() - 1 ) {
			if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
				out	<< "\nThe line search failed again even with a new decomposition!"
					<< " (k = " << algo.state().k() << ")\n"
					<< "We quit!\n";
			}
			throw;
		}

		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
			out	<< "\nSelecting a new decomposition..."
				<< " (k = " << algo.state().k() << ")\n";
		}

		last_ls_failure_k_ = s.k();
		return new_decomp_strategy().new_decomposition(algo,step_poss,type,assoc_step_poss);
	}
	return false;	// will never be executed.
}

void LineSearchFailureNewDecompositionSelection_Step::print_step(
	const Algorithm& algo
	,poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
	,std::ostream& out, const std::string& L
	) const
{
	out
		<< L << "do line search step : " << typeid(line_search_step()).name() << std::endl;
	line_search_step().print_step(algo,step_poss,type,assoc_step_poss,out,L+"    ");
	out
		<< L << "end line search step\n"
		<< L << "if thrown line_search_failure then\n"
		<< L << "  if line search failed at the last iteration also then\n"
		<< L << "    throw line_search_failure\n"
		<< L << "  end\n"
		<< L << "  new decomposition selection : " << typeid(new_decomp_strategy()).name() << std::endl
		;
	new_decomp_strategy().print_new_decomposition(
		rsqp_algo(algo),step_poss,type,assoc_step_poss,out, L + "    " );
	out
		<< L << "  end new decomposition selection\n"
		<< L << "end\n"
		;
}

} // end namespace ReducedSpaceSQPPack
