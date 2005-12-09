// ////////////////////////////////////////////////////////////////////////////
// MoochoPack_NLPAlgo.hpp
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

#ifndef RSQP_ALGO_H
#define RSQP_ALGO_H

#include "MoochoPack_NLPAlgoInterface.hpp"
#include "MoochoPack_NLPAlgoContainer.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "IterationPack_Algorithm.hpp"
#include "StandardAggregationMacros.hpp"

namespace MoochoPack {

///
/** rSQP Algorithm control class.
  */
class NLPAlgo
	: public NLPAlgoInterface
	, public IterationPack::Algorithm
{
public:

	/** @name Public Types */
	//@{

	//@}

	/// Constructs with no step, added_step, pre_step, post_step, state, or decomp_sys objects added.
	NLPAlgo();

	/// <<std aggr>> members for algo_cntr
	STANDARD_AGGREGATION_MEMBERS( NLPAlgoContainer, algo_cntr )

	/// <<std aggr>> members for nlp
	STANDARD_AGGREGATION_MEMBERS( NLP, nlp )

	///
	NLPAlgoState& rsqp_state()
	{	return dynamic_cast<NLPAlgoState&>(state()); }

	///
	const NLPAlgoState& rsqp_state() const
	{	return dynamic_cast<const NLPAlgoState&>(state()); }

	///
	void do_step_first(Algorithm::poss_type first_step_poss)
	{	first_step_poss_ = first_step_poss; }

	/** @name Overridden form rSQPAlgoInteface */
	//@{	
	
	///
	const NLPAlgoState& retrieve_state() const;

	///
	/** This is the main control function for the rSQP algorithm.
	  *
	  * This function basically just calls Algorithm::do_algorithm(...).
	  */
	NLPSolverClientInterface::EFindMinReturn dispatch();

	///
	void interface_print_algorithm(std::ostream& out) const;
	///
	void interface_set_algo_timing( bool algo_timing );
	///
	bool interface_algo_timing() const;
	///
	void interface_print_algorithm_times( std::ostream& out ) const;

	//@}

	/// overridden from Algorihth.

	///
	void print_algorithm(std::ostream& out) const;

protected:

	// First step to execute
	Algorithm::poss_type first_step_poss_;

};	// end class NLPAlgo

}	// end namespace MoochoPack

#endif	// RSQP_ALGO_H
