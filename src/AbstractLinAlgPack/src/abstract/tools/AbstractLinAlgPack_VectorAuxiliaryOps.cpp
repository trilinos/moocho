// ////////////////////////////////////////////////////////////////////
// VectorAuxiliaryOps.cpp
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

#include "AbstractLinAlgPack/src/abstract/tools/VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorMutable.hpp"
#include "RTOp_ROp_max.h"
#include "RTOp_ROp_max_near_feas_step.h"
#include "RTOp_ROp_max_rel_step.h"
#include "RTOp_ROp_num_bounded.h"
#include "RTOp_ROp_combined_nu_comp_err.h"
#include "RTOp_ROp_comp_err_with_mu.h"
#include "RTOp_ROp_fraction_to_boundary.h"
#include "RTOp_ROp_fraction_to_zero_boundary.h"
#include "RTOp_ROp_log_bound_barrier.h"
#include "RTOp_TOp_correct_multipliers.h"
#include "RTOp_ROp_max_inequ_viol.h"
#include "RTOp_TOp_multiplier_step.h"
#include "RTOp_TOp_force_in_bounds.h"
#include "RTOp_TOp_ele_wise_sqrt.h"
#include "RTOp_TOp_inv_of_difference.h"
#include "RTOp_TOp_max_vec_scalar.h"
#include "RTOp_TOp_max_abs_vec_scalar.h"
#include "RTOpCppC.hpp"
#include "ThrowException.hpp"

namespace {

// log_bound_barrier
static RTOpPack::RTOpC          log_bound_barrier_op;
static RTOpPack::ReductTarget   log_bound_barrier_targ;
// combined_nu_comp_err
static RTOpPack::RTOpC          combined_nu_comp_err_op;
static RTOpPack::ReductTarget   combined_nu_comp_err_targ;
// combined_nu_comp_err_lower
static RTOpPack::RTOpC          combined_nu_comp_err_lower_op;
static RTOpPack::ReductTarget   combined_nu_comp_err_lower_targ;
// combined_nu_comp_err_upper
static RTOpPack::RTOpC          combined_nu_comp_err_upper_op;
static RTOpPack::ReductTarget   combined_nu_comp_err_upper_targ;
// comp_err_with_mu
static RTOpPack::RTOpC          comp_err_with_mu_op;
static RTOpPack::ReductTarget   comp_err_with_mu_targ;
// maximum near feasible step
static RTOpPack::RTOpC          max_near_feas_step_op;
static RTOpPack::ReductTarget   max_near_feas_step_targ;
// fraction to boundary rule
static RTOpPack::RTOpC          fraction_to_boundary_op;
static RTOpPack::ReductTarget   fraction_to_boundary_targ;
// fraction to zero boundary rule
static RTOpPack::RTOpC          fraction_to_zero_boundary_op;
static RTOpPack::ReductTarget   fraction_to_zero_boundary_targ;
// maximum relative step
static RTOpPack::RTOpC          max_rel_step_op;
static RTOpPack::ReductTarget   max_rel_step_targ;
// number of bounded elements
static RTOpPack::RTOpC          num_bounded_op;
static RTOpPack::ReductTarget   num_bounded_targ;
// force in bounds
static RTOpPack::RTOpC          force_in_bounds_op;
// force in bounds with buffer
static RTOpPack::RTOpC          force_in_bounds_buffer_op;
// inv_of_difference
static RTOpPack::RTOpC          inv_of_difference_op;
// correct_multipliers
static RTOpPack::RTOpC          correct_lower_bound_multipliers_op;
static RTOpPack::RTOpC          correct_upper_bound_multipliers_op;
// multipliers step
static RTOpPack::RTOpC          lowerbound_multipliers_step_op;
static RTOpPack::RTOpC          upperbound_multipliers_step_op;
// element wise square root
static RTOpPack::RTOpC          ele_wise_sqrt_op;

// Simple class for an object that will initialize the RTOp_Server.
class init_rtop_server_t {
public:
	init_rtop_server_t() {
		// Operator and target obj for log_bound_barrier
		if(0>RTOp_ROp_log_bound_barrier_construct(&log_bound_barrier_op.op() ))
			assert(0);
		log_bound_barrier_op.reduct_obj_create(&log_bound_barrier_targ);
		// Operator and target obj for combined_nu_comp_err
		if(0>RTOp_ROp_combined_nu_comp_err_construct(&combined_nu_comp_err_op.op() ))
			assert(0);
		combined_nu_comp_err_op.reduct_obj_create(&combined_nu_comp_err_targ);
		// Operator and target obj for combined_nu_comp_err_lower
		if(0>RTOp_ROp_combined_nu_comp_err_one_only_construct(&combined_nu_comp_err_lower_op.op() ))
			assert(0);
		combined_nu_comp_err_lower_op.reduct_obj_create(&combined_nu_comp_err_lower_targ);
		// Operator and target obj for combined_nu_comp_err_upper
		if(0>RTOp_ROp_combined_nu_comp_err_one_only_construct(&combined_nu_comp_err_upper_op.op() ))
			assert(0);
		combined_nu_comp_err_upper_op.reduct_obj_create(&combined_nu_comp_err_upper_targ);
		// Operator and target obj for comp_err_with_mu
		if(0>RTOp_ROp_comp_err_with_mu_construct(0.0, 0.0, &comp_err_with_mu_op.op()))
			assert(0);
		comp_err_with_mu_op.reduct_obj_create(&comp_err_with_mu_targ);
		// Operator and target obj for max_near_feas_step
		if(0>RTOp_ROp_max_near_feas_step_construct(0.0,&max_near_feas_step_op.op() ))
			assert(0);
		max_near_feas_step_op.reduct_obj_create(&max_near_feas_step_targ);
		// Operator and target obj for max_rel_step
		if(0>RTOp_ROp_max_rel_step_construct(&max_rel_step_op.op() ))
			assert(0);
		max_rel_step_op.reduct_obj_create(&max_rel_step_targ);
		// Operator and target obj for fraction to boundary
		if(0>RTOp_ROp_fraction_to_boundary_construct(0.99,&fraction_to_boundary_op.op() ))
			assert(0);
		fraction_to_boundary_op.reduct_obj_create(&fraction_to_boundary_targ);
		// Operator and target obj for fraction to zero boundary
		if(0>RTOp_ROp_fraction_to_zero_boundary_construct(0.99, &fraction_to_zero_boundary_op.op()))
			assert(0);
		fraction_to_zero_boundary_op.reduct_obj_create(&fraction_to_zero_boundary_targ);
		// Operator and target obj for num_bounded
		if(0>RTOp_ROp_num_bounded_construct(0.0,&num_bounded_op.op() ))
			assert(0);
		num_bounded_op.reduct_obj_create(&num_bounded_targ);
		// Operator force_in_bounds
		if(0>RTOp_TOp_force_in_bounds_construct( &force_in_bounds_op.op() ))
			assert(0);
		// Operator force_in_bounds_buffer
		if(0>RTOp_TOp_force_in_bounds_buffer_construct( 0.01, 0.001, &force_in_bounds_buffer_op.op() ))
			assert(0);
		// Operator inv_of_difference
		if(0>RTOp_TOp_inv_of_difference_construct( 1.0,  &inv_of_difference_op.op()))
			assert(0);
		// correct_lower_bounds_multipliers
		if(0>RTOp_TOp_Correct_Multipliers_construct( -1e50, 0,  &correct_lower_bound_multipliers_op.op()))
			assert(0);
		// correct_upper_bounds_multipliers
		if(0>RTOp_TOp_Correct_Multipliers_construct( 1e50, 1,  &correct_upper_bound_multipliers_op.op()))
			assert(0);
		// lower_bounds_multipliers step
		if(0>RTOp_TOp_multiplier_step_construct( 1.0, -1.0,  &lowerbound_multipliers_step_op.op()))
			assert(0);
		// upper_bounds_multipliers step
		if(0>RTOp_TOp_multiplier_step_construct( 1.0, 1.0,  &upperbound_multipliers_step_op.op()))
			assert(0);
		// ele_wise_sqrt
		if(0>RTOp_TOp_ele_wise_sqrt_construct( &ele_wise_sqrt_op.op()))
			assert(0);
	}
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

AbstractLinAlgPack::value_type
AbstractLinAlgPack::max( const Vector& v )
{
	RTOpPack::RTOpC          op;
	RTOpPack::ReductTarget   reduct_obj;
	assert(0==RTOp_ROp_max_construct(&op.op()));
	op.reduct_obj_create(&reduct_obj);
	const Vector* vecs[1] = { &v };
	apply_op(op,1,vecs,0,NULL,reduct_obj.obj());
	return RTOp_ROp_max_val(reduct_obj.obj());
}

std::pair<AbstractLinAlgPack::value_type,AbstractLinAlgPack::value_type>
AbstractLinAlgPack::max_near_feas_step(
	const Vector& x, const Vector& d
	,const Vector& xl, const Vector& xu
	,value_type max_bnd_viol
	)
{
	const int num_vecs = 4;
	const Vector*
		vecs[num_vecs] = { &xl, &x, &d, &xu };
	if(0!=RTOp_ROp_max_near_feas_step_set_beta( max_bnd_viol, &max_near_feas_step_op.op() )) assert(0);
	max_near_feas_step_targ.reinit();
	apply_op(
		max_near_feas_step_op, num_vecs, vecs, 0, NULL
		,max_near_feas_step_targ.obj() );
	RTOp_ROp_max_near_feas_step_reduct_obj_t
		u = RTOp_ROp_max_near_feas_step_val(max_near_feas_step_targ.obj());
	return std::pair<value_type,value_type>(u.alpha_pos,u.alpha_neg);
}

AbstractLinAlgPack::value_type
AbstractLinAlgPack::max_rel_step(
	const Vector& x, const Vector& d
	)
{
	const int num_vecs = 2;
	const Vector*
		vecs[num_vecs] = { &x, &d };
	max_rel_step_targ.reinit();
	apply_op(
		max_rel_step_op, num_vecs, vecs, 0, NULL
		,max_rel_step_targ.obj() );
	return RTOp_ROp_max_rel_step_val(max_rel_step_targ.obj());
}


AbstractLinAlgPack::value_type
AbstractLinAlgPack::fraction_to_boundary(
	const value_type tau,
	const Vector& x,
	const Vector& d,
	const Vector& xl,
	const Vector& xu
	)
{
	if(0!=RTOp_ROp_fraction_to_boundary_init( tau, &fraction_to_boundary_op.op() )) assert(0);
	fraction_to_boundary_targ.reinit();

	const int num_vecs = 4;
	const Vector*
		vecs[num_vecs] = { &x, &d, &xl, &xu };

	apply_op(
		fraction_to_boundary_op, num_vecs, vecs, 0, NULL
		,fraction_to_boundary_targ.obj() );
	
	return RTOp_ROp_fraction_to_boundary_val(fraction_to_boundary_targ.obj());
}

AbstractLinAlgPack::value_type
AbstractLinAlgPack::fraction_to_zero_boundary(
  const value_type tau,
  const Vector& x,
  const Vector& d
  )
{
	if(0!=RTOp_ROp_fraction_to_zero_boundary_init( tau, &fraction_to_zero_boundary_op.op() )) assert(0);
	fraction_to_zero_boundary_targ.reinit();

	const int num_vecs = 2;
	const Vector*
		vecs[num_vecs] = { &x, &d };

	apply_op(
		fraction_to_zero_boundary_op, num_vecs, vecs, 0, NULL
		,fraction_to_zero_boundary_targ.obj() );

	return RTOp_ROp_fraction_to_zero_boundary_val(fraction_to_zero_boundary_targ.obj());
}

AbstractLinAlgPack::size_type
AbstractLinAlgPack:: num_bounded(
	const Vector& xl, const Vector& xu
	,value_type inf_bound
	)
{
	const int num_vecs = 2;
	const Vector*
		vecs[num_vecs] = { &xl, &xu };
	if(0!=RTOp_ROp_num_bounded_set_inf_bnd( inf_bound, &num_bounded_op.op() )) assert(0);
	num_bounded_targ.reinit();
	apply_op(
		num_bounded_op, num_vecs, vecs, 0, NULL
		,num_bounded_targ.obj() );
	return RTOp_ROp_num_bounded_val(num_bounded_targ.obj());
}

AbstractLinAlgPack::value_type
AbstractLinAlgPack::log_bound_barrier(
	const Vector    &x
	,const Vector   &xl
	,const Vector   &xu
	)
{
	log_bound_barrier_targ.reinit();
	const int num_vecs = 3;
	const Vector*
		vecs[num_vecs] = { &x, &xl, &xu };
	apply_op(
		log_bound_barrier_op, num_vecs, vecs, 0, NULL
		,log_bound_barrier_targ.obj()
		);

	return RTOp_ROp_log_bound_barrier_val(log_bound_barrier_targ.obj());
}

AbstractLinAlgPack::value_type
AbstractLinAlgPack::combined_nu_comp_err(
	const Vector     &v
	,const Vector    &x
	,const Vector   &xl
	,const Vector   &xu
	)
{
	combined_nu_comp_err_targ.reinit();;
	const int num_vecs = 4;
	const Vector*
		vecs[num_vecs] = {&v, &x, &xl, &xu };
	apply_op(
		combined_nu_comp_err_op, num_vecs, vecs, 0, NULL
		,combined_nu_comp_err_targ.obj()
		);
	return RTOp_ROp_combined_nu_comp_err_val(combined_nu_comp_err_targ.obj());
}

AbstractLinAlgPack::value_type
AbstractLinAlgPack::combined_nu_comp_err_lower(
	const Vector    &v
	,const Vector   &x
	,const Vector   &xl
	)
{
	combined_nu_comp_err_lower_targ.reinit();
	const int num_vecs = 3;
	const Vector*
		vecs[num_vecs] = {&v, &xl, &x};
	apply_op(
		combined_nu_comp_err_lower_op, num_vecs, vecs, 0, NULL
		,combined_nu_comp_err_lower_targ.obj()
		);
	return RTOp_ROp_combined_nu_comp_err_one_only_val(combined_nu_comp_err_lower_targ.obj());
}


AbstractLinAlgPack::value_type
AbstractLinAlgPack::combined_nu_comp_err_upper(
	const Vector    &v
	,const Vector   &x
	,const Vector   &xu
	)
{
	combined_nu_comp_err_upper_targ.reinit();
	const int num_vecs = 3;
	const Vector*
		vecs[num_vecs] = {&v, &xu, &x};
	apply_op(
		combined_nu_comp_err_upper_op, num_vecs, vecs, 0, NULL
		,combined_nu_comp_err_upper_targ.obj()
		);
	return RTOp_ROp_combined_nu_comp_err_one_only_val(combined_nu_comp_err_upper_targ.obj());
}

AbstractLinAlgPack::value_type
AbstractLinAlgPack::IP_comp_err_with_mu(
  const value_type mu
  ,const value_type inf_bound
  ,const Vector &x
  ,const Vector &xl
  ,const Vector &xu
  ,const Vector &vl
  ,const Vector &vu
  )
{
	if(0!=RTOp_ROp_comp_err_with_mu_init(mu, inf_bound, &comp_err_with_mu_op.op())) assert(0);
	comp_err_with_mu_targ.reinit();
	const int num_vecs = 5;
	const Vector*
		vecs[num_vecs] = {&x, &xl, &xu, &vl, &vu};
	apply_op(
		comp_err_with_mu_op, num_vecs, vecs, 0, NULL
		,comp_err_with_mu_targ.obj()
		);
	return RTOp_ROp_comp_err_with_mu_val(comp_err_with_mu_targ.obj());
}

bool AbstractLinAlgPack::max_inequ_viol(
	const AbstractLinAlgPack::Vector   &v
	,const AbstractLinAlgPack::Vector  &vL
	,const AbstractLinAlgPack::Vector  &vU
	,AbstractLinAlgPack::size_type     *max_viol_i
	,AbstractLinAlgPack::value_type    *max_viol
	,AbstractLinAlgPack::value_type    *v_i
	,int                               *bnd_type
	,AbstractLinAlgPack::value_type    *vLU_i
	)
{
	RTOpPack::RTOpC          op;
	RTOpPack::ReductTarget   reduct_obj;
	RTOp_ROp_max_inequ_viol_construct(&op.op());
	op.reduct_obj_create(&reduct_obj);
	const int num_vecs = 3;
	const Vector*
		vecs[num_vecs] = { &v, &vL, &vU };
	apply_op(
		op, num_vecs, vecs, 0, NULL
		,reduct_obj.obj()
		);
	const RTOp_ROp_max_inequ_viol_reduct_obj_t
		ro = RTOp_ROp_max_inequ_viol_val(reduct_obj.obj());
	*max_viol_i = ro.max_viol_i;
	*max_viol   = ro.max_viol;
	*v_i        = ro.v_i;
	*bnd_type   = ro.bnd_type;
	*vLU_i      = ro.vLU_i;
	return *max_viol_i > 0.0;
}

void AbstractLinAlgPack::force_in_bounds(
	const Vector& xl, const Vector& xu
	,VectorMutable* x
	)
{
#ifdef _DEBUG
	THROW_EXCEPTION(x==NULL,std::logic_error,"force_in_bounds(...), Error");
#endif
	const Vector*  vecs[2]      = { &xl, &xu };
	VectorMutable* targ_vecs[1] = { x };
	apply_op(force_in_bounds_op,2,vecs,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}


void AbstractLinAlgPack::force_in_bounds_buffer(
  const value_type rel_push,
  const value_type abs_push,
  const Vector& xl, 
  const Vector& xu,
  VectorMutable* x
  )
{
	if(0!=RTOp_TOp_force_in_bounds_buffer_init( rel_push, abs_push, &force_in_bounds_buffer_op.op())) assert(0);
	const Vector*  vecs[2]      = { &xl, &xu };
	VectorMutable* targ_vecs[1] = { x };
	apply_op(force_in_bounds_buffer_op,2,vecs,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}


void AbstractLinAlgPack::inv_of_difference(
  const value_type      alpha
  ,const Vector   &v0
  ,const Vector   &v1
  ,VectorMutable  *z
  )
{
	if(0!=RTOp_TOp_inv_of_difference_init( alpha, &inv_of_difference_op.op())) assert(0);
	const Vector*  vecs[2]      = { &v0, &v1 };
	VectorMutable* targ_vecs[1] = { z };
	apply_op(inv_of_difference_op,2,vecs,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

void AbstractLinAlgPack::correct_lower_bound_multipliers(
  const Vector       &xl
  ,const value_type  inf_bound_limit
  ,VectorMutable     *vl
  )
{
	if(0!=RTOp_TOp_Correct_Multipliers_init( inf_bound_limit, 0, &correct_lower_bound_multipliers_op.op())) assert(0);
	const Vector*  vecs[1]      = { &xl };
	VectorMutable* targ_vecs[1] = { vl };
	apply_op(correct_lower_bound_multipliers_op,1,vecs,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

void AbstractLinAlgPack::correct_upper_bound_multipliers(
  const Vector       &xu
  ,const value_type  inf_bound_limit
  ,VectorMutable     *vu
  )
{
	if(0!=RTOp_TOp_Correct_Multipliers_init( inf_bound_limit, 1, &correct_upper_bound_multipliers_op.op())) assert(0);
	const Vector*  vecs[1]      = { &xu };
	VectorMutable* targ_vecs[1] = { vu };
	apply_op(correct_upper_bound_multipliers_op,1,vecs,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

void AbstractLinAlgPack::lowerbound_multipliers_step(
  const value_type mu,
  const Vector& invXl,
  const Vector& vl,
  const Vector& d_k,
  VectorMutable* dvl_k
  )
{
	if(0!=RTOp_TOp_multiplier_step_init(mu, -1.0, &lowerbound_multipliers_step_op.op())) assert(0);
	const Vector*  vecs[]      = { &invXl, &vl, &d_k };
	VectorMutable* targ_vecs[] = { dvl_k };
	apply_op(lowerbound_multipliers_step_op,3,vecs,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

void AbstractLinAlgPack::upperbound_multipliers_step(
  const value_type mu,
  const Vector& invXu,
  const Vector& vu,
  const Vector& d_k,
  VectorMutable* dvu_k
  )
{
 	if(0!=RTOp_TOp_multiplier_step_init(mu, 1.0, &upperbound_multipliers_step_op.op())) assert(0);
	const Vector*  vecs[] = { &invXu, &vu, &d_k };
	VectorMutable* targ_vecs[] = { dvu_k };
	apply_op(upperbound_multipliers_step_op,3,vecs,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

void AbstractLinAlgPack::ele_wise_sqrt(
  VectorMutable* z
  )
{
	VectorMutable* targ_vecs[] = { z };
	apply_op(ele_wise_sqrt_op,0,NULL,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);  
}

void AbstractLinAlgPack::max_vec_scalar(
	value_type        min_ele
	,VectorMutable    *y
	)
{
	RTOpPack::RTOpC op;
	if(0!=RTOp_TOp_max_vec_scalar_construct(min_ele,&op.op())) assert(0);
	VectorMutable* targ_vecs[] = { y };
	apply_op(op,0,NULL,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}

void AbstractLinAlgPack::max_abs_vec_scalar(
	value_type        min_ele
	,VectorMutable    *y
	)
{
	RTOpPack::RTOpC op;
	if(0!=RTOp_TOp_max_abs_vec_scalar_construct(min_ele,&op.op())) assert(0);
	VectorMutable* targ_vecs[] = { y };
	apply_op(op,0,NULL,1,targ_vecs,RTOp_REDUCT_OBJ_NULL);
}
