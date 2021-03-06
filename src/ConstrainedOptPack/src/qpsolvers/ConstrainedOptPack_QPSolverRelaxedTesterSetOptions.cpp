// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <assert.h>
#include <math.h>

#include "ConstrainedOptPack_QPSolverRelaxedTesterSetOptions.hpp"

// Define the options
namespace {

  const int local_num_options = 6;

  enum local_EOptions {
    OPT_WARNING_TOL
    ,OPT_ERROR_TOL
    ,FEAS_WARNING_TOL
    ,FEAS_ERROR_TOL
    ,COMP_WARNING_TOL
    ,COMP_ERROR_TOL
  };

  const char* local_SOptions[local_num_options]	= {
      "opt_warning_tol"
      ,"opt_error_tol"
      ,"feas_warning_tol"
      ,"feas_error_tol"
      ,"comp_warning_tol"
      ,"comp_error_tol"
  };

}

namespace ConstrainedOptPack {

QPSolverRelaxedTesterSetOptions::QPSolverRelaxedTesterSetOptions(
        QPSolverRelaxedTester* target
      , const char opt_grp_name[] )
  :	OptionsFromStreamPack::SetOptionsFromStreamNode(
        opt_grp_name, local_num_options, local_SOptions )
    , OptionsFromStreamPack::SetOptionsToTargetBase<
      QPSolverRelaxedTester >( target )
{}

void QPSolverRelaxedTesterSetOptions::setOption(
  int option_num, const std::string& option_value )
{
  typedef QPSolverRelaxedTester target_t;
  switch( (local_EOptions)option_num ) {
      case OPT_WARNING_TOL:
      target().opt_warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case OPT_ERROR_TOL:
      target().opt_error_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case FEAS_WARNING_TOL:
      target().feas_warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case FEAS_ERROR_TOL:
      target().feas_error_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case COMP_WARNING_TOL:
      target().comp_warning_tol(std::fabs(std::atof(option_value.c_str())));
      break;
      case COMP_ERROR_TOL:
      target().comp_error_tol(std::fabs(std::atof(option_value.c_str())));
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);	// Local error only?
  }
}

}	// end namespace ConstrainedOptPack
