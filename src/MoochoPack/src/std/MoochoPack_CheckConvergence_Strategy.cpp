// ////////////////////////////////////////////////////////////////////////////
// CheckConvergence_Strategy.cpp
//
// Copyright (C) 2001
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

#include "MoochoPack/src/std/CheckConvergence_Strategy.hpp"
#include "MoochoMoreUtilities/src/StringToBool.hpp"

namespace MoochoPack {

///*******************************************
//  CheckConvergence_Strategy
///*******************************************


CheckConvergence_Strategy::CheckConvergence_Strategy(
  EOptErrorCheck opt_error_check,
  EScaleKKTErrorBy scale_opt_error_by,
  EScaleKKTErrorBy scale_feas_error_by,
  EScaleKKTErrorBy scale_comp_error_by,
  bool scale_opt_error_by_Gf
  )
	:
	opt_error_check_(opt_error_check),
	scale_opt_error_by_(scale_opt_error_by),
	scale_feas_error_by_(scale_feas_error_by),
	scale_comp_error_by_(scale_comp_error_by),
	scale_opt_error_by_Gf_(scale_opt_error_by_Gf) 
	{}


///*******************************************
//  CheckConvergence_StrategySetOptions
///*******************************************

// Define the options
namespace {

const int local_num_options = 4;

enum local_EOptions 
	{
	SCALE_OPT_ERROR_BY,
	SCALE_FEAS_ERROR_BY,
	SCALE_COMP_ERROR_BY,
	SCALE_OPT_ERROR_BY_GF
	};

const char* local_SOptions[local_num_options] = 
	{
	"scale_opt_error_by",
	"scale_feas_error_by",
	"scale_comp_error_by",
	"scale_opt_error_by_Gf",
	};

} // end namespace

CheckConvergence_StrategySetOptions::CheckConvergence_StrategySetOptions(
  CheckConvergence_Strategy* target,
  const char opt_grp_name[] )
	:  
	OptionsFromStreamPack::SetOptionsFromStreamNode(
	  opt_grp_name, local_num_options, local_SOptions ),
	OptionsFromStreamPack::SetOptionsToTargetBase<
	  CheckConvergence_Strategy >( target )
	{}


void CheckConvergence_StrategySetOptions::setOption(
  int option_num,
  const std::string& option_value )
	{
	using OptionsFromStreamPack::StringToBool;

	typedef CheckConvergence_Strategy target_t;
	switch( (local_EOptions)option_num ) 
		{
		case SCALE_OPT_ERROR_BY:
		case SCALE_FEAS_ERROR_BY:
		case SCALE_COMP_ERROR_BY:
			{
			const std::string &option = option_value.c_str();
			CheckConvergence_Strategy::EScaleKKTErrorBy scale_by = target_t::SCALE_BY_ONE;

			if( option == "SCALE_BY_ONE" )
				{ scale_by = target_t::SCALE_BY_ONE; }
			else if( option == "SCALE_BY_NORM_2_X" )
				{ scale_by = target_t::SCALE_BY_NORM_2_X; }
			else if( option == "SCALE_BY_NORM_INF_X" )
				{ scale_by = target_t::SCALE_BY_NORM_INF_X; }
			else
				{
				throw std::invalid_argument( "Error, incorrect value for "
											 "\"scale_kkt_error_by\".  Only the options "
											 "SCALE_BY_ONE, SCALE_BY_NORM_2_X, and SCALE_BY_NORM_INF_X "
											 "are available" );
				}


			if ((local_EOptions) option_num == SCALE_OPT_ERROR_BY)
				{
				target().scale_opt_error_by(scale_by);
				}
			else if ((local_EOptions) option_num == SCALE_FEAS_ERROR_BY)
				{
				target().scale_feas_error_by(scale_by);
				}
			else if ((local_EOptions) option_num == SCALE_COMP_ERROR_BY)
				{
				target().scale_comp_error_by(scale_by);
				}
			else
				{
				TEST_FOR_EXCEPTION( true,
								 std::logic_error,
								 "Unaccounted for option_num in CheckConvergence_Strategy.cpp"
				  );
				}

			break;
			}
		case SCALE_OPT_ERROR_BY_GF: 
			{
			target().scale_opt_error_by_Gf(
			  StringToBool( "scale_opt_error_by_Gf", option_value.c_str() ) );
			break;
			}
		default:
			assert(0);	// Local error only?
		}
	}

} // end namespace MoochoPack
