// ////////////////////////////////////////////////////////////////////////////
// CrossTermExactStd_Step.cpp
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

#include "../std/CrossTermExactStd_Step.hpp"
#include "../moocho_algo_conversion.hpp"
#include "IterationPack/src/print_algorithm_step.hpp"
#include "ConstrainedOptPack/src/VectorWithNorms.h"
#include "AbstractLinAlgPack/src/MatrixOp.hpp"
#include "DenseLinAlgPack/src/LinAlgOpPack.hpp"
#include "DenseLinAlgPack/src/DVectorClass.hpp"
#include "DenseLinAlgPack/src/DVectorOut.hpp"

namespace LinAlgOpPack {
	using AbstractLinAlgPack::Vp_StMtV;
}

bool MoochoPack::CrossTermExactStd_Step::do_step(Algorithm& _algo
	, poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using LinAlgOpPack::V_MtV;
	using DenseLinAlgPack::norm_inf;

	NLPAlgo	&algo	= rsqp_algo(_algo);
	NLPAlgoState	&s		= algo.rsqp_state();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using IterationPack::print_algorithm_step;
		print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
	}

	// tmp = HL * Ypy
	DVector tmp;
	V_MtV( &tmp, s.HL().get_k(0), BLAS_Cpp::no_trans, s.Ypy().get_k(0)() );
	// w = Z' * tmp
	V_MtV( &s.w().set_k(0).v(), s.Z().get_k(0), BLAS_Cpp::trans, tmp() );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||w||inf = "	<< s.w().get_k(0).norm_inf() << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nw_k =\n" << s.w().get_k(0)();
	}

	return true;
}

void MoochoPack::CrossTermExactStd_Step::print_step( const Algorithm& algo
	, poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Evaluate the exact reduced QP cross term\n"
		<< L << "w_k = Z_k' * HL_k * Ypy_k\n";
}
