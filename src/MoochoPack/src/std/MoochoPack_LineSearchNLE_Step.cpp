// ////////////////////////////////////////////////////////////////////////////
// LineSearchNLE_Step.cpp
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

#include "MoochoPack/src/std/LineSearchNLE_Step.hpp"
#include "MoochoPack/src/MoochoPackExceptions.hpp"
#include "MoochoPack/src/moocho_algo_conversion.hpp"
#include "IterationPack/src/print_algorithm_step.hpp"
#include "ConstrainedOptPack/src/globalization/MeritFuncNLESqrResid.hpp"
#include "ConstrainedOptPack/src/globalization/MeritFuncCalcNLE.hpp"
#include "ConstrainedOptPack/src/globalization/MeritFuncCalc1DQuadratic.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorMutable.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorStdOps.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorOut.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/LinAlgOpPack.hpp"
#include "ThrowException.hpp"

namespace MoochoPack {

LineSearchNLE_Step::LineSearchNLE_Step(
	const direct_line_search_ptr_t& direct_line_search
	)
	:direct_line_search_(direct_line_search)
{}

bool LineSearchNLE_Step::do_step(
	Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
	,poss_type assoc_step_poss
	)
{
	using AbstractLinAlgPack::Vp_StV;
	using LinAlgOpPack::V_VpV;

	NLPAlgo	        &algo   = rsqp_algo(_algo);
	NLPAlgoState    &s      = algo.rsqp_state();
	NLP             &nlp    = algo.nlp();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using IterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	const size_type
		n  = nlp.n(),
		m  = nlp.m();

	THROW_EXCEPTION( m == 0, std::logic_error, "LineSearchNLE_Step::do_step(...) : Error!" );

	// /////////////////////////////////////////
	// Set references to iteration quantities
	//
	// Set k+1 first then go back to get k+0 to ensure
	// we have backward storage!

	IterQuantityAccess<value_type>
		&alpha_iq = s.alpha();
	IterQuantityAccess<VectorMutable>
		&x_iq    = s.x(),
		&d_iq    = s.d(),
		&c_iq    = s.c();
	
	VectorMutable        &x_kp1   = x_iq.get_k(+1);
	const Vector         &x_k     = x_iq.get_k(0);
	VectorMutable        &c_kp1   = c_iq.get_k(+1);
	const Vector         &c_k     = c_iq.get_k(0);
	const Vector         &d_k     = d_iq.get_k(0);
	value_type           &alpha_k = alpha_iq.get_k(0);

	// //////////////////////////////////////
	// Build the merit function
	
	ConstrainedOptPack::MeritFuncNLESqrResid  phi_c;

	// /////////////////////////////////////
	// Compute Dphi_k, phi_kp1 and phi_k

	// phi_k, phi_kp1
	const value_type  phi_k   = phi_c.value(c_k);
	value_type        phi_kp1 = phi_c.value(c_kp1);
	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		out << "\nBegin definition of NLP merit function phi_c.value(c(x)):\n";
		phi_c.print_merit_func( out, "    " );
		out	<< "end definition of the NLP merit funciton\n";
	}
	// Dphi_k
	phi_c.calc_deriv(c_k);
	const value_type
		Dphi_k  = phi_c.deriv();
	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		out	<< "\nDphi_k = "	<< Dphi_k << std::endl;
	}
	THROW_EXCEPTION(
		Dphi_k >= 0, LineSearchFailure
		,"LineSearchNLE_Step::do_step(...) : " 
		"Error, d_k is not a descent direction for the merit function "
		"since Dphi_k = " << Dphi_k << " >= 0" );

	// //////////////////////
	// Do the line search
	
	nlp.unset_quantities();
	nlp.set_c( &c_kp1 );
	ConstrainedOptPack::MeritFuncCalcNLE    phi_c_calc( &phi_c, &nlp );
	const Vector*                           xd[2] = { &x_k, &d_k };
	MeritFuncCalc1DQuadratic                phi_calc_1d( phi_c_calc, 1, xd, &x_kp1 );
	
	if( !direct_line_search().do_line_search(
			phi_calc_1d, phi_k, &alpha_k, &phi_kp1
			,( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS)
			   ? &out : static_cast<std::ostream*>(0)	)
			)
		)
	{
		// The line search failed!
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) )
			out
				<< "\nThe maximum number of linesearch iterations has been exceeded "
				<< "(k = " << algo.state().k() << ")\n"
				<< "(phi_k - phi_kp1)/phi_k = " << ((phi_k - phi_kp1)/phi_k)
				<< "\nso we will reject the step and declare a line search failure.\n";
		THROW_EXCEPTION(
			true, LineSearchFailure
			,"LineSearchNLE_Step::do_step(): Line search failure" );
	}

	nlp.unset_quantities();
	
	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		out	<< "\nalpha_k                = " << alpha_k;
		out << "\n||x_kp1||inf           = " << x_kp1.norm_inf();
		out << "\n||c_kp1||inf           = " << c_kp1.norm_inf();
		out << "\nphi_kp1                = " << phi_kp1;
		out << std::endl;
	}
	
	if( (int)olevel >= (int)PRINT_VECTORS ) {
		out << "\nx_kp1 =\n" << x_kp1;
		out << "\nc_kp1 =\n" << c_kp1;
	}

	return true;
}

void LineSearchNLE_Step::print_step(
	const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
	,std::ostream& out, const std::string& L
	) const
{
	out
		<< L << "*** Preform a line search for c(x_k + alpha_k*d_k) along the full space search direction d_k.\n"
		<< L << "ToDo: Fill this in!\n";
}

} // end namespace MoochoPack
