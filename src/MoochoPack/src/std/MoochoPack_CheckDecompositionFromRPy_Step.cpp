// ////////////////////////////////////////////////////////////////////////////
// CheckDecompositionFromRPy_Step.cpp
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

#include "ReducedSpaceSQPPack/include/std/CheckDecompositionFromRPy_Step.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/VectorWithOp.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"

namespace ReducedSpaceSQPPack {

CheckDecompositionFromRPy_Step::CheckDecompositionFromRPy_Step(
	const new_decomp_strategy_ptr_t   &new_decomp_strategy
	,value_type                       max_decomposition_cond_change_frac
	)
	:new_decomp_strategy_(new_decomp_strategy)
	,max_decomposition_cond_change_frac_(max_decomposition_cond_change_frac)
{
	reset();
}

void CheckDecompositionFromRPy_Step::reset() {
	beta_min_ = std::numeric_limits<value_type>::max();
}

// Overridden

bool CheckDecompositionFromRPy_Step::do_step( Algorithm& _algo, poss_type step_poss
	, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss )
{
	rSQPAlgo                &algo       = rsqp_algo(_algo);
	rSQPState               &s          = algo.rsqp_state();
	const Range1D           equ_decomp  = s.equ_decomp();
	EJournalOutputLevel     olevel      = algo.algo_cntr().journal_output_level();
	std::ostream            &out        = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	bool select_new_decomposition = false;

	// Compute: resid = (Gc(decomp)'*Y) * py + c(decomp)
	const VectorWithOp                  &py_k       = s.py().get_k(0);
	const VectorWithOp                  &c_k        = s.c().get_k(0);
	VectorWithOp::vec_ptr_t             c_decomp_k  = c_k.sub_view(equ_decomp);
	VectorWithOpMutable::vec_mut_ptr_t  resid       = c_decomp_k->space().create_member();

	// resid = R*py + c(equ_decomp)
	LinAlgOpPack::V_MtV( resid.get(), s.R().get_k(0), BLAS_Cpp::no_trans, py_k );
	LinAlgOpPack::Vp_V( resid.get(), *c_decomp_k );

	const value_type
		small_num    = std::numeric_limits<value_type>::min(),
		nrm_resid    = resid->norm_inf(),
		nrm_c_decomp = c_decomp_k->norm_inf(),
		beta         = nrm_resid / (nrm_c_decomp+small_num);

	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		out	<< "\nbeta = ||R*py_k + c_k(decomp)||inf / (||c_k(decomp)||inf + small_number)"
			<< "\n     = "<<nrm_resid<<" / ("<<nrm_c_decomp<<" + "<<small_num<<")"
			<< "\n     = " << beta << std::endl;
	}

	// Check to see if a new basis was selected or not
	IterQuantityAccess<index_type>
		&num_basis_iq = s.num_basis();
	if( num_basis_iq.updated_k(0) ) {
		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
			out	<< "\nnum_basis_k was updated so the basis changed so we will skip this check\n"
				<< "    reset min ||R*py+c||/||c|| to current value\n";
		beta_min_ = beta;
		return true;
	}
	
	if( beta != 0.0 ) {
		if( beta < beta_min_ ) {
			beta_min_ = beta;
		}
		else {
			if( beta / beta_min_ > max_decomposition_cond_change_frac() ) {
				if( (int)olevel >= (int)PRINT_BASIC_ALGORITHM_INFO ) {
					out
						<< "\nbeta_change = ( ||R*py+c||/||c|| = " << beta
						<< " ) / ( min ||R*py+c||/||c|| = " << beta_min_ << " )\n"
						<< "              = " << (beta/beta_min_) << " > max_decomposition_cond_change_frac = "
						<< max_decomposition_cond_change_frac()
						<< "\n\nSelecting a new decomposition"
						<< " (k = " << algo.state().k() << ") ...\n";
				}
				select_new_decomposition = true;
			}
		}
		if(select_new_decomposition) {
			reset();
			return new_decomp_strategy().new_decomposition(algo,step_poss,type,assoc_step_poss);
		}
	}
	return true;
}

void CheckDecompositionFromRPy_Step::print_step( const Algorithm& algo, poss_type step_poss
	, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Try to detect when the decomposition is becomming illconditioned\n"
		<< L << "default: beta_min = inf\n"
		<< L << "         max_decomposition_cond_change_frac = " << max_decomposition_cond_change_frac() << std::endl
		<< L << "beta = norm_inf(R*py_k + c_k(decomp)) / (norm_inf(c_k(decomp))+small_number)\n"
		<< L << "select_new_decomposition = false\n"
		<< L << "if num_basis_k is updated then\n"
		<< L << "  beta_min = beta\n"
		<< L << "end\n"
		<< L << "if beta < beta_min then\n"
		<< L << "  beta_min = beta\n"
		<< L << "else\n"
		<< L << "  if beta/ beta_min > max_decomposition_cond_change_frac then\n"
		<< L << "        select_new_decomposition = true\n"
		<< L << "    end\n"
		<< L << "end\n"
		<< L << "if select_new_decomposition == true then\n"
		<< L << "    new decomposition selection : " << typeid(new_decomp_strategy()).name() << std::endl
		;
	new_decomp_strategy().print_new_decomposition(
		rsqp_algo(algo),step_poss,type,assoc_step_poss,out, L + "    " );
	out
		<< L << "    end new decomposition selection\n"
		<< L << "end\n"
		;
}

}	// end namespace ReducedSpaceSQPPack 
