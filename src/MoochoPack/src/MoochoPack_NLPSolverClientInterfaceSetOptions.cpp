// ////////////////////////////////////////////////////////////////
// NLPSolverClientInterfaceSetOptions.cpp
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

#include <assert.h>
#include <math.h>

#include "MoochoPack/src/NLPSolverClientInterfaceSetOptions.hpp"
#include "StringToBool.hpp"

// Define the options
namespace {

	const int local_num_options = 10;

	enum local_EOptions {
        MAX_ITER,
        MAX_RUN_TIME,
        OPT_TOL,
        FEAS_TOL,
		COMP_TOL,
        STEP_TOL,
		JOURNAL_OUTPUT_LEVEL,
		JOURNAL_PRINT_DIGITS,
		CHECK_RESULTS,
		CALC_CONDITIONING
	};

	const char* local_SOptions[local_num_options]	= {
        ("max_iter"),
        ("max_run_time"),
        ("opt_tol"),
        ("feas_tol"),
		("comp_tol"),
        ("step_tol"),
		("journal_output_level"),
		("journal_print_digits"),
		("check_results"),
		("calc_conditioning")
	};

}

namespace MoochoPack {

NLPSolverClientInterfaceSetOptions::NLPSolverClientInterfaceSetOptions(
			  NLPSolverClientInterface* target
			, const char opt_grp_name[] )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  opt_grp_name, local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			NLPSolverClientInterface >( target )
{}

void NLPSolverClientInterfaceSetOptions::set_option(
	int option_num, const std::string& option_value )
{
	namespace ofsp = OptionsFromStreamPack;
	using ofsp::StringToBool;

	typedef NLPSolverClientInterface target_t;
	switch( (local_EOptions)option_num ) {
	    case MAX_ITER:
			target().max_iter(::abs(::atoi(option_value.c_str())));
			break;
		case MAX_RUN_TIME:
			target().max_run_time(::fabs(::atof(option_value.c_str())));
			break;
		case OPT_TOL:
			target().opt_tol(::fabs(::atof(option_value.c_str())));
			break;
		case FEAS_TOL:
			target().feas_tol(::fabs(::atof(option_value.c_str())));
			break;
	    case COMP_TOL:
		    target().comp_tol(::fabs(::atof(option_value.c_str())));
		    break;
		case STEP_TOL:
			target().step_tol(::fabs(::atof(option_value.c_str())));
			break;
		case JOURNAL_OUTPUT_LEVEL:
		{
			if( option_value == "PRINT_NOTHING" )
				target().journal_output_level(PRINT_NOTHING);
			else if( option_value == "PRINT_BASIC_ALGORITHM_INFO" )
				target().journal_output_level(PRINT_BASIC_ALGORITHM_INFO);
			else if( option_value == "PRINT_ALGORITHM_STEPS" )
				target().journal_output_level(PRINT_ALGORITHM_STEPS);
			else if( option_value == "PRINT_ACTIVE_SET" )
				target().journal_output_level(PRINT_ACTIVE_SET);
			else if( option_value == "PRINT_VECTORS" )
				target().journal_output_level(PRINT_VECTORS);
			else if( option_value == "PRINT_ITERATION_QUANTITIES" )
				target().journal_output_level(PRINT_ITERATION_QUANTITIES);
			else
				throw std::invalid_argument( "NLPSolverClientInterfaceSetOptions::set_option(...) : "
					"Error, incorrect value for \"journal_output_level\"." );
			break;
		}
	    case JOURNAL_PRINT_DIGITS:
			target().journal_print_digits(::abs(::atoi(option_value.c_str())));
			break;
		case CHECK_RESULTS:
			target().check_results(
				StringToBool( "check_results", option_value.c_str() )
				);
			break;
		case CALC_CONDITIONING:
			target().calc_conditioning(
				StringToBool( "calc_conditioning", option_value.c_str() )
				);
			break;
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace MoochoPack 
