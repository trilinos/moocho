// ////////////////////////////////////////////////////////////////////////////
// MeritFunc_ModifiedL1LargerSteps_AddedStep.hpp
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

#ifndef MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_H
#define MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_H

#include "../rSQPAlgo_Step.h"
#include "ConstrainedOptPack/src/globalization/MeritFuncNLP.hpp"
#include "StandardCompositionMacros.hpp"
#include "StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

///
/** This function increases the penalty parameters of the modifed L1
  * merit function to allow for larger steps by taking advantage
  * of constraints that are reduced for a full step.
  */
class MeritFunc_ModifiedL1LargerSteps_AddedStep
	: public rSQPAlgo_Step
{
public:

	/// <<std comp>> members for merit_func
	STANDARD_COMPOSITION_MEMBERS(MeritFuncNLP,merit_func)

	///
	MeritFunc_ModifiedL1LargerSteps_AddedStep(
		  const merit_func_ptr_t& merit_func
		, value_type	eta
		, int			after_k_iter				= 3
		, value_type	obj_increase_threshold		= 1e-4
		, value_type	max_pos_penalty_increase	= 1.0
		, value_type	pos_to_neg_penalty_increase	= 1.0
		, value_type	incr_mult_factor			= 1e-4 );

	/// eta.
	STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,eta)

	/// after_k_iter.
	STANDARD_MEMBER_COMPOSITION_MEMBERS(int,after_k_iter)

	/// obj_increase_threshold.
	STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,obj_increase_threshold)

	/// max_pos_penalty_increase
	STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,max_pos_penalty_increase)

	/// pos_to_neg_penalty_increase
	STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,pos_to_neg_penalty_increase)

	/// incr_mult_factor
	STANDARD_MEMBER_COMPOSITION_MEMBERS(value_type,incr_mult_factor)

	// ///////////////////////////////
	// Overridden from AlgorithmStep

	///
	bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
		, poss_type assoc_step_poss);

	///
	void print_step( const Algorithm& algo, poss_type step_poss
		, IterationPack::EDoStepType type, poss_type assoc_step_poss
		, std::ostream& out, const std::string& leading_str ) const;

private:
	// not defined and not to be called.
	MeritFunc_ModifiedL1LargerSteps_AddedStep();
	MeritFunc_ModifiedL1LargerSteps_AddedStep(const MeritFunc_ModifiedL1LargerSteps_AddedStep&);
	MeritFunc_ModifiedL1LargerSteps_AddedStep& operator=(const MeritFunc_ModifiedL1LargerSteps_AddedStep&);
	
};	// end class MeritFunc_ModifiedL1LargerSteps_AddedStep

}	// end namespace MoochoPack 

#endif	// MERIT_FUNC_MODIFIED_L1_LARGER_STEPS_ADDED_STEP_H
