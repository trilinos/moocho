// ////////////////////////////////////////////////////////////////
// LineSearchFilter_StepSetOptions.cpp
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

#include "MoochoPack/src/std/LineSearchFilter_StepSetOptions.hpp"
#include "StringToBool.hpp"
#include "ThrowException.hpp"

// Define the options
namespace {

const int local_num_options = 9;

enum local_EOptions 
	{
		GAMMA_THETA
		,GAMMA_F
		,GAMMA_ALPHA
		,DELTA
		,S_THETA
		,S_F
		,THETA_SMALL_FACT
		,ETA_F
		,BACK_TRACK_FRAC
	};

const char* local_SOptions[local_num_options] = 
	{
		"gamma_theta"
		,"gamma_f"
		,"gamma_alpha"
		,"delta"
		,"s_theta"
		,"s_f"
		,"theta_small_fact"
		,"eta_f"
		,"back_track_frac"
	};

}

namespace MoochoPack {

LineSearchFilter_StepSetOptions::LineSearchFilter_StepSetOptions(
  LineSearchFilter_Step* target
  , const char opt_grp_name[] )
	:
	OptionsFromStreamPack::SetOptionsFromStreamNode(
	  opt_grp_name, local_num_options, local_SOptions ),
	OptionsFromStreamPack::SetOptionsToTargetBase< LineSearchFilter_Step >( target )
	{
	}

void LineSearchFilter_StepSetOptions::set_option( 
  int option_num, const std::string& option_value )
	{
	using OptionsFromStreamPack::StringToBool;
  
	typedef LineSearchFilter_Step target_t;
	switch( (local_EOptions)option_num ) 
		{
		case GAMMA_THETA:
			target().gamma_theta(::atof(option_value.c_str()));
			break;
		case GAMMA_F:
			target().gamma_f(::atof(option_value.c_str()));
			break;
		case GAMMA_ALPHA:
			target().gamma_alpha(::atof(option_value.c_str()));
			break;
		case DELTA:
			target().delta(::atof(option_value.c_str()));
			break;
		case S_THETA:
			target().s_theta(::atof(option_value.c_str()));
			break;
		case S_F:
			target().s_f(::atof(option_value.c_str()));
			break;
		case THETA_SMALL_FACT:
			target().theta_small_fact(::atof(option_value.c_str()));
			break;
		case ETA_F:
			target().eta_f(::atof(option_value.c_str()));
			break;
		case BACK_TRACK_FRAC:
			target().back_track_frac(::atof(option_value.c_str()));
			break;
		default:
			assert(0);	// Local error only?
		}
	}

}	// end namespace MoochoPack 
