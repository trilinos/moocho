// ////////////////////////////////////////////////////////////////
// DirectSparseSolverMA28SetOptions.cpp
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

#ifdef SPARSE_SOLVER_PACK_USE_MA28

#include <assert.h>

#include "DirectSparseSolverMA28SetOptions.hpp"
#include "StringToBool.hpp"

// Define the options
namespace {

	const int local_num_options = 8;

	enum local_EOptions {
		ESTIMATED_FILLIN_RATIO
	  	,U
		,GROW
		,TOL
		,NSRCH
		,LBIG
		,PRINT_MA28_OUTPUTS
		,OUTPUT_FILE_NAME
	};

	const char* local_SOptions[local_num_options]	= {
		"estimated_fillin_ratio"
	  	,"u"
		,"grow"
		,"tol"
		,"nsrch"
		,"lbig"
		,"print_ma28_outputs"
		,"output_file_name"
	};

}

namespace AbstractLinAlgPack {

DirectSparseSolverMA28SetOptions::DirectSparseSolverMA28SetOptions(
			DirectSparseSolverMA28* qp_solver )
	:	OptionsFromStreamPack::SetOptionsFromStreamNode(
			  "DirectSparseSolverMA28", local_num_options, local_SOptions )
		, OptionsFromStreamPack::SetOptionsToTargetBase<
			DirectSparseSolverMA28 >( qp_solver )
{}

void DirectSparseSolverMA28SetOptions::set_option(
	int option_num, const std::string& option_value )
{
	using OptionsFromStreamPack::StringToBool;

	switch( (local_EOptions)option_num ) {
		case ESTIMATED_FILLIN_RATIO: {
			target().estimated_fillin_ratio( ::atof( option_value.c_str() ) );
			break;
		}
		case U: {
			target().u( ::atof( option_value.c_str() ) );
			break;
		}
		case GROW: {
			target().grow( StringToBool( "grow", option_value.c_str() ) );
			break;
		}
		case TOL: {
			target().tol( ::atof( option_value.c_str() ) );
			break;
		}
		case NSRCH: {
			target().nsrch( ::atoi( option_value.c_str() ) );
			break;
		}
		case LBIG: {
			target().lbig( StringToBool( "lbig", option_value.c_str() ) );
			break;
		}
		case PRINT_MA28_OUTPUTS: {
			target().print_ma28_outputs( StringToBool( "grow", option_value.c_str() ) );
			break;
		}
		case OUTPUT_FILE_NAME: {
			if( option_value == "NONE" )
				target().output_file_name( "" );
			else
				target().output_file_name( option_value );
			break;
		}
		default:
			assert(0);	// Local error only?
	}
}

}	// end namespace AbstractLinAlgPack 

#endif // SPARSE_SOLVER_PACK_USE_MA28
