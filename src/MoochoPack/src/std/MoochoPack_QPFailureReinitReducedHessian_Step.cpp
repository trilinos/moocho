// ////////////////////////////////////////////////////////////////////////////
// QPFailureReinitReducedHessian_Step.cpp
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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>
#include <typeinfo>

#include "ReducedSpaceSQPPack/include/std/QPFailureReinitReducedHessian_Step.h"
#include "ReducedSpaceSQPPack/include/std/rSQPAlgorithmStepNames.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"

ReducedSpaceSQPPack::QPFailureReinitReducedHessian_Step::QPFailureReinitReducedHessian_Step(
	const null_space_step_ptr_t& null_space_step
	)
	:null_space_step_(null_space_step)
	,last_qp_failure_k_(-100) // has not failed yet.
{}

bool ReducedSpaceSQPPack::QPFailureReinitReducedHessian_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	try {
		return null_space_step().do_step(_algo,step_poss,type,assoc_step_poss);
	}
	catch(const QPFailure& qp_excpt) {
		rSQPAlgo	&algo	= rsqp_algo(_algo);
		rSQPState	&s		= algo.rsqp_state();
		NLP			&nlp	= algo.nlp();

		EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
		std::ostream& out = algo.track().journal_out();

		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
			out	<< "\nQP failed! "
				<< " (k = " << algo.state().k() << ")\n"
				<< "QPFailure description: " << qp_excpt.what() << "\n";
		}
		if( s.k() >= algo.algo_cntr().max_iter() ) {
			if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
				out	<< "\nThe maximum number of rSQP iterations\n"
					<< " have been exceeded so quit "
					<< " (k = " << algo.state().k() << ")\n";
			}
			algo.terminate(false);
			return false;
		}
		if( last_qp_failure_k_ == s.k() ) {
			if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
				out	<< "\nThe QP failed again even with a new reduced Hessian rHL_k!"
					<< " (k = " << algo.state().k() << ")\n"
					<< "We quit!\n";
			}
			throw;
		}
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
			out	<< "\nWiping out all memory for rHL and going back to reinitalize it ..."
				<< " (k = " << algo.state().k() << ")\n";
		}
		last_qp_failure_k_ = s.k(); // Remember this for later!
		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
			out	<< "Wipe out all update rHL_{k} for all k\n"
				<< "goto ReducedHessian\n";
		}
		s.rHL().set_all_not_updated();
		algo.do_step_next( ReducedHessian_name );
		return false;
	}
	return false;	// will never be executed.
}

void ReducedSpaceSQPPack::QPFailureReinitReducedHessian_Step::print_step(
	  const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "do null space step : " << typeid(null_space_step()).name() << std::endl;
	null_space_step().print_step(algo,step_poss,type,assoc_step_poss,out,L+"    ");
	out
		<< L << "end null space step\n"
		<< L << "if QPFailure was thrown then\n"
		<< L << "    if QP failed already then\n"
		<< L << "        rethrow QPFailure\n"
		<< L << "    end\n"
		<< L << "    if k > max_iter then\n"
		<< L << "        terminate the algorithm!\n"
		<< L << "    end\n"
		<< L << "    set all rHL_{k} to not updated\n"
		<< L << "    goto ReducedHessian\n"
		<< L << "end\n"
		;
}
