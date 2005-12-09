// ////////////////////////////////////////////////////////////////////////////
// NLPAlgoContainer.cpp
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

//#define RELEASE_TRACE

#include <assert.h>

#include <iostream>	// used for debugging the Release version.

#include "MoochoPack_NLPAlgoContainer.hpp"
#include "MoochoPack_NLPAlgoInterface.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "NLPInterfacePack_NLP.hpp"
#include "Teuchos_TestForException.hpp"

namespace {

void report_final_point( const MoochoPack::NLPAlgoState& s, const bool optimal, NLPInterfacePack::NLP* nlp )
{
	const AbstractLinAlgPack::size_type
		m  = nlp->m(),
		nb = nlp->num_bounded_x();
	const IterationPack::IterQuantityAccess<AbstractLinAlgPack::VectorMutable>
		&x_iq = s.x();
	if( x_iq.updated_k(0) ) {
		nlp->report_final_solution(
			x_iq.get_k(0)                                                      // x
			,( m  && s.lambda().updated_k(0) ) ? &s.lambda().get_k(0)  : NULL  // lambda
			,( nb && s.nu().updated_k(0)     ) ? &s.nu().get_k(0)      : NULL  // nu
			, optimal                                                          // optimal = false
			);
	}
}

} // end namespace

namespace MoochoPack {

// Overridden from rSQPAlgoClient interface

void NLPAlgoContainer::set_config(const config_ptr_t& config)
{
	algo_ = Teuchos::null; // Remove our reference to the current (configured?) algorithm.
	config_ = config;
}

NLPAlgoContainer::config_ptr_t&
NLPAlgoContainer::get_config()
{	
	return config_;
}

const NLPAlgoContainer::config_ptr_t&
NLPAlgoContainer::get_config() const
{	
	return config_;
}

NLPAlgoConfig&
NLPAlgoContainer::config()
{	
	return *config_;
}

const NLPAlgoConfig&
NLPAlgoContainer::config() const
{	
	return *config_;
}

NLPSolverClientInterface::EFindMinReturn
NLPAlgoContainer::find_min()
{
	config().init_algo(&algo());
	EFindMinReturn solve_return;
	try {
		solve_return = algo().dispatch();
	}
	catch(...) {
		report_final_point(algo().retrieve_state(),false,&nlp());
		throw;
	}
	report_final_point(algo().retrieve_state(),solve_return==NLPSolverClientInterface::SOLUTION_FOUND,&nlp());
	return solve_return;
}

void NLPAlgoContainer::configure_algorithm(std::ostream* trase_out)
{
	assert_valid_setup();
	if(!get_algo().get())
		config().config_algo_cntr(this,trase_out);
}

void NLPAlgoContainer::print_algorithm(std::ostream& out) const
{
	algo().interface_print_algorithm(out);
}

void NLPAlgoContainer::set_algo_timing( bool algo_timing )
{
	algo().interface_set_algo_timing(algo_timing);
}

bool NLPAlgoContainer::algo_timing() const
{
	return algo().interface_algo_timing();
}

void NLPAlgoContainer::print_algorithm_times(
	std::ostream& out ) const
{
	algo().interface_print_algorithm_times(out);
}

void NLPAlgoContainer::assert_valid_setup() const {
	TEST_FOR_EXCEPTION(
		!get_nlp().get(), NLPSolverClientInterface::InvalidSetup
		,"NLPAlgoContainer::assert_valid_setup() : The NLP object has not been set" );
	TEST_FOR_EXCEPTION(
		!get_track().get(), NLPSolverClientInterface::InvalidSetup
		,"NLPAlgoContainer::assert_valid_setup() : The AlgorithmTracker object has not been set" );
	TEST_FOR_EXCEPTION(
		!get_config().get(), NLPSolverClientInterface::InvalidSetup
		,"NLPAlgoContainer::assert_valid_setup() : The NLPAlgoConfig object has not been set" );
}

} // end namespace MoochoPack
