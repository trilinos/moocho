// ////////////////////////////////////////////////////////////////////////////
// DampenCrossTermStd_Step.cpp
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
#include <sstream>
#include <limits>

#include "../std/DampenCrossTermStd_Step.hpp"
#include "../rsqp_algo_conversion.hpp"
#include "GeneralIterationPack/src/print_algorithm_step.hpp"
#include "ConstrainedOptimizationPack/src/VectorWithNorms.h"
#include "SparseLinAlgPack/src/MatrixWithOpFactorized.hpp"
#include "DenseLinAlgPack/src/LinAlgOpPack.hpp"
#include "DenseLinAlgPack/src/DVectorClass.hpp"
#include "DenseLinAlgPack/src/DVectorOut.hpp"

ReducedSpaceSQPPack::DampenCrossTermStd_Step::DampenCrossTermStd_Step(const value_type& frac_descent)
	: frac_descent_(frac_descent)
{}

bool ReducedSpaceSQPPack::DampenCrossTermStd_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using SparseLinAlgPack::V_InvMtV;
	using DenseLinAlgPack::norm_inf;
	using DenseLinAlgPack::dot;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
	}

	if( s.w().updated_k(0) ) {

		// inv(rHL_k) * rGf_k
		const DVectorSlice rGf_k = s.rGf().get_k(0)();
		DVector Inv_rHL_rGf;
		V_InvMtV( &Inv_rHL_rGf, dynamic_cast<MatrixWithOpFactorized&>(s.rHL().get_k(0))
			, BLAS_Cpp::no_trans, rGf_k );
		
		const value_type
			small_num			= 1e-20,
			rGfT_Inv_rHL_rGf	= dot( Inv_rHL_rGf(), rGf_k ),					// rGf_k'*inv(rHL_k)*rGf_k
			rGfT_Inv_rHL_w		= dot( Inv_rHL_rGf(), s.w().get_k(0)() ),		// rGf_k'*inv(rHL_k)*w_k
			term				= -(1.0-frac_descent()) * (rGfT_Inv_rHL_rGf + 2*small_num)
									/ (rGfT_Inv_rHL_w + small_num);

		if( rGfT_Inv_rHL_w >= 0.0 ) {
			// We know that the descent property will be satisfied for all zeta_k > 0
			// so set zeta_k = 1
			s.zeta().set_k(0) = 1.0;
		}
		else {
			// For some zeta_k > 0 the descent property will be violated so we may have to
			// cut zeta_k back from 1.
			s.zeta().set_k(0) = std::_MIN( term, 1.0 );
		}

		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out	<< "\nterm1 = rGf_k'*inv(rHL_k)*rGf_k                = "	<<	rGfT_Inv_rHL_rGf;
			out	<< "\nterm2 = rGf_k'*inv(rHL_k)*w_k                  = "	<<	rGfT_Inv_rHL_w;
			out	<< "\n(1-frac_descent)*(term1+2*small)/(term2+small) = "	<<  term;
			out	<< "\nzeta_k                                         = "	<<  s.zeta().get_k(0)
				<< std::endl;
		}

		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
			out	<< "\ninv(rHL_k)*rGf_k = "	<<	Inv_rHL_rGf();
		}

		if( rGfT_Inv_rHL_rGf < 0.0 ) {
			std::ostringstream omsg;
			omsg
				<< "Error, rGf_k'*inv(rHL_k)*rGf_k = " << rGfT_Inv_rHL_rGf << " < 0.0 and therefore "
				<< "the reduced Hessian rHL_k can not be positive definite";
			if( (int)(olevel) >= (int)(PRINT_ALGORITHM_STEPS) ) {
				out << omsg.str();
			}
			throw std::runtime_error( std::string("DampenCrossTermStd_Step::do_step(...) : ")
										+ omsg.str() );
		}
	}

	return true;
}

void ReducedSpaceSQPPack::DampenCrossTermStd_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Compute the dampening parameter for the reduced QP cross term w_k\n"
		<< L << "default: frac_descent = " << frac_descent() << std::endl
		<< L << "if w_k is update then\n"
		<< L << "    find zeta_k s.t.\n"
		<< L << "        Gf_k'*Z_k*pz_k ~\n"
		<< L << "           - zeta_k * rGf_k'*inv(rHL_k)*w_k - rGf_k'*inv(rHL_k)*rGf_k\n"
		<< L << "             <= - frac_descent * rGf_k'*inv(rHL_k)*rGf_k\n"
		<< L << "end\n";
}
