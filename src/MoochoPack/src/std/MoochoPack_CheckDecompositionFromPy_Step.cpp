// ////////////////////////////////////////////////////////////////////////////
// CheckDecompositionFromPy_Step.cpp
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

#include <typeinfo>
#include <limits>

#include "MoochoPack/src/std/CheckDecompositionFromPy_Step.hpp"
#include "MoochoPack/src/moocho_algo_conversion.hpp"
#include "IterationPack/src/print_algorithm_step.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/Vector.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/LinAlgOpPack.hpp"

namespace MoochoPack {

CheckDecompositionFromPy_Step::CheckDecompositionFromPy_Step(
	const new_decomp_strategy_ptr_t   &new_decomp_strategy
	,value_type                       max_decomposition_cond_change_frac
	)
	:new_decomp_strategy_(new_decomp_strategy)
	,max_cond_( 0.01 / std::numeric_limits<value_type>::epsilon() )
	,max_decomposition_cond_change_frac_( max_decomposition_cond_change_frac )
{
	reset();
}

void CheckDecompositionFromPy_Step::reset() {
	beta_min_ = std::numeric_limits<value_type>::max();
}

// Overridden

bool CheckDecompositionFromPy_Step::do_step(
	Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
	)
{
	NLPAlgo                &algo       = rsqp_algo(_algo);
	NLPAlgoState               &s          = algo.rsqp_state();
	const Range1D           equ_decomp  = s.equ_decomp();
	EJournalOutputLevel     olevel      = algo.algo_cntr().journal_output_level();
	std::ostream            &out        = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using IterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	bool select_new_decomposition = false;

	const value_type
		small_num = std::numeric_limits<value_type>::min(),
		beta = s.py().get_k(0).norm_inf() / (s.c().get_k(0).sub_view(equ_decomp)->norm_inf()+small_num);

	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		out	<< "\nbeta = ||py||/||c|| = " << beta << std::endl;
	}

	// Check to see if a new basis was selected or not
	IterQuantityAccess<index_type>
		&num_basis_iq = s.num_basis();
	if( num_basis_iq.updated_k(0) ) {
		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
			out	<< "\nnum_basis_k was updated so the basis changed so we will skip this check\n"
				<< "    reset min ||py||/||c|| to current value + 1\n";
		beta_min_ = beta + 1.0;
		return true;
	}

	if( beta + 1.0 < beta_min_ ) {
		beta_min_ = beta + 1.0;
	}
	else {
		if( (beta + 1.0)/ beta_min_ > max_decomposition_cond_change_frac() ) {
			if( (int)olevel >= (int)PRINT_BASIC_ALGORITHM_INFO ) {
				out	<< "Select a new decomposition"
						<< " (k = " << algo.state().k() << ")\n"
					<< "beta_change = ( ||py||/||c|| = " << beta
						<< " ) / ( min ||py||/||c|| = " << beta_min_ << " )\n"
					<< "beta_change = " << (beta/beta_min_) << " > max_decomposition_cond_change_frac = "
						<< max_decomposition_cond_change_frac() << std::endl;
			}
			select_new_decomposition = true;
		}
	}
	if( !select_new_decomposition && beta > max_cond() ) {
		if( (int)olevel >= (int)PRINT_BASIC_ALGORITHM_INFO ) {
			out	<< "\nConditioning of decomposition matrix is > " << beta
				<< " > max_cond = " << max_cond() << std::endl
				<< "Selecting a new decomposition ... "
				<< " (k = " << algo.state().k() << ")\n";
		}
		select_new_decomposition = true;
	}

	if(select_new_decomposition) {
		reset();
		return new_decomp_strategy().new_decomposition(algo,step_poss,type,assoc_step_poss);
	}

	return true;		
}

void CheckDecompositionFromPy_Step::print_step(
	const Algorithm& algo, poss_type step_poss
	,IterationPack::EDoStepType type, poss_type assoc_step_poss
	,std::ostream& out, const std::string& L ) const
{
	out
		<< L << "default: beta_min = inf\n"
		<< L << "         max_decomposition_cond_change_frac = " << max_decomposition_cond_change_frac() << std::endl
		<< L << "         max_cond = 0.01 * mach_eps\n"
		<< L << "beta = norm_inf(py_k) / (norm_inf(c_k(equ_decomp))+small_number)\n"
		<< L << "select_new_decomposition = false\n"
		<< L << "if num_basis_k is updated then\n"
		<< L << "  beta_min = beta + 1\n"
		<< L << "end\n"
		<< L << "if beta + 1 < beta_min then\n"
		<< L << "  beta_min = beta + 1\n"
		<< L << "else\n"
		<< L << "  if (beta + 1) / beta_min > max_decomposition_cond_change_frac then\n"
		<< L << "    select_new_decomposition = true\n"
		<< L << "  end\n"
		<< L << "end\n"
		<< L << "if beta > max_cond then\n"
		<< L << "  select_new_decomposition = true\n"
		<< L << "end\n"
		<< L << "if select_new_decomposition == true then\n"
		<< L << "  new decomposition selection : " << typeid(new_decomp_strategy()).name() << std::endl
		;
	new_decomp_strategy().print_new_decomposition(
		rsqp_algo(algo),step_poss,type,assoc_step_poss,out, L + "    " );
	out
		<< L << "  end new decomposition selection\n"
		<< L << "end\n"
		;
}

}	// end namespace MoochoPack 
