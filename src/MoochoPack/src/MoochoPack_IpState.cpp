// ////////////////////////////////////////////////////////////////////
// ipState.cpp
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
//
#include "ReducedSpaceSQPPack/src/ipState.hpp"
#include "AbstractLinAlgPack/src/MatrixSymDiagStd.hpp"

// Iteration Quantity Strings
extern const std::string ReducedSpaceSQPPack::barrier_parameter_name = "barrier_parameter";
extern const std::string ReducedSpaceSQPPack::barrier_obj_name = "barrier_obj";
extern const std::string ReducedSpaceSQPPack::grad_barrier_obj_name = "grad_barrier_obj";
extern const std::string ReducedSpaceSQPPack::e_tol_name = "e_tol";
extern const std::string ReducedSpaceSQPPack::Vu_name = "Vu";
extern const std::string ReducedSpaceSQPPack::Vl_name = "Vl";
extern const std::string ReducedSpaceSQPPack::invXu_name = "invXu";
extern const std::string ReducedSpaceSQPPack::invXl_name = "invXl";
extern const std::string ReducedSpaceSQPPack::rHB_name = "rHB";
extern const std::string ReducedSpaceSQPPack::B_name = "B";
extern const std::string ReducedSpaceSQPPack::Sigma_name = "Sigma";
extern const std::string ReducedSpaceSQPPack::w_sigma_name = "w_sigma";
extern const std::string ReducedSpaceSQPPack::dvl_name = "dvl";
extern const std::string ReducedSpaceSQPPack::dvu_name = "dvu";
extern const std::string ReducedSpaceSQPPack::alpha_vl_name = "alpha_vl";
extern const std::string ReducedSpaceSQPPack::alpha_vu_name = "alpha_vu";

namespace ReducedSpaceSQPPack {

ipState::ipState(
  const decomp_sys_ptr_t& decomp_sys
  ,const vec_space_ptr_t& space_x
  ,const vec_space_ptr_t& space_c
  ,const vec_space_ptr_t& space_h
  ,const vec_space_ptr_t& space_range
  ,const vec_space_ptr_t& space_null
  )
	:
	rSQPState(decomp_sys, space_x, space_c, space_h, space_range, space_null)
	{
	}

ipState::~ipState()
	{
	}

///********** Iteration Quantities **************

STATE_SCALAR_IQ_DEF(ipState, barrier_parameter, barrier_parameter_name)

STATE_SCALAR_IQ_DEF(ipState, barrier_obj, barrier_obj_name)

STATE_VECTOR_IQ_DEF(ipState, grad_barrier_obj, grad_barrier_obj_name, get_space_x(), VST_SPACE_X)

STATE_SCALAR_IQ_DEF(ipState, e_tol, e_tol_name)

STATE_IQ_DEF(ipState, MatrixSymDiagStd, Vu, Vu_name)

STATE_IQ_DEF(ipState, MatrixSymDiagStd, Vl, Vl_name)

STATE_IQ_DEF(ipState, MatrixSymDiagStd, invXu, invXu_name)

STATE_IQ_DEF(ipState, MatrixSymDiagStd, invXl, invXl_name)

STATE_IQ_DEF(ipState, MatrixSymOp, rHB, rHB_name)

STATE_IQ_DEF(ipState, MatrixSymOp, B, B_name)

STATE_IQ_DEF(ipState, MatrixSymDiagStd, Sigma, Sigma_name)

STATE_VECTOR_IQ_DEF(ipState, w_sigma, w_sigma_name, get_space_null(), VST_SPACE_NULL )  
STATE_VECTOR_IQ_DEF(ipState, dvl, dvl_name, get_space_x(), VST_SPACE_X)
STATE_VECTOR_IQ_DEF(ipState, dvu, dvu_name, get_space_x(), VST_SPACE_X)

STATE_SCALAR_IQ_DEF(ipState, alpha_vl, alpha_vl_name)
STATE_SCALAR_IQ_DEF(ipState, alpha_vu, alpha_vu_name)
} // end namespace ReducedSpaceSQPPack
