// /////////////////////////////////////////////
// ExampleNLPFirstOrderDirectRTOps.c
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
#include <malloc.h>

#include "ExampleNLPFirstOrderDirectRTOps.h"
#include "RTOpPack/src/RTOp_obj_null_vtbl.h"
#include "RTOpPack/src/RTOp_obj_index_vtbl.h"

// Implementation for RTOp_TOp_explnlp2_c_eval

static int explnlp2_c_eval_apply_op(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
	, const int num_vecs, const struct RTOp_SubVector vecs[]
	, const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
	, RTOp_ReductTarget targ_obj )
{
	// c
	size_t                 sub_dim;
	RTOp_value_type        *c_val;
	ptrdiff_t              c_val_s;
	// xD
	const RTOp_value_type  *xD_val;
	ptrdiff_t              xD_val_s;
	// xI
	const RTOp_value_type  *xI_val;
	ptrdiff_t              xI_val_s;

	register RTOp_index_type  k;

	//
	// Validate the input
	//
	if( num_vecs != 2 || vecs == NULL )
		return RTOp_ERR_INVALID_NUM_VECS;
	if( num_targ_vecs != 1 || targ_vecs == NULL )
		return RTOp_ERR_INVALID_NUM_TARG_VECS;
	if( targ_vecs[0].sub_dim != vecs[0].sub_dim
		|| targ_vecs[0].sub_dim != vecs[1].sub_dim
		|| targ_vecs[0].global_offset != vecs[0].global_offset
		|| targ_vecs[0].global_offset != vecs[1].global_offset )
		return RTOp_ERR_INCOMPATIBLE_VECS;
	if( vecs[0].indices || vecs[0].sub_nz == 0 || vecs[1].indices || vecs[1].sub_nz == 0 )
		return RTOp_ERR_UNSUPPORTED_VEC_TYPE;

	//
	// Get pointers to data
	//

	// c
	sub_dim       = targ_vecs[0].sub_dim;
	c_val         = targ_vecs[0].values;
	c_val_s       = targ_vecs[0].values_stride;
	// xD
	xD_val         = vecs[0].values;
	xD_val_s       = vecs[0].values_stride;
	// xI
	xI_val         = vecs[1].values;
	xI_val_s       = vecs[1].values_stride;

	//
	// Compute c(j) = xI(i) * ( xD(i) - 1 ) - 10 * xD(i)
	//

	if( c_val_s == 1 && xD_val_s == 1 && xI_val_s == 1 ) {
		// Slightly faster loop for unit stride vectors
		for( k = 0; k < sub_dim; ++k, ++xI_val )
			*c_val++ = (*xD_val++) * (*xI_val - 1.0) - 10.0 * (*xI_val);
	}
	else {
		// More general implementation for non-unit strides
		for( k = 0; k < sub_dim; ++k, c_val+=c_val_s, xD_val+=xD_val_s, xI_val+=xI_val_s )
			*c_val = (*xD_val) * (*xI_val - 1.0) - 10.0 * (*xI_val);
	}

	return 0; // success?
}

const char RTOp_TOp_explnlp2_c_eval_name[] = "TOp_explnlp2_c_eval";

const struct RTOp_RTOp_vtbl_t RTOp_TOp_explnlp2_c_eval_vtbl =
{
	&RTOp_obj_null_vtbl  // Use a null object for instance data
	,&RTOp_obj_null_vtbl // use null type for target object
	,NULL
	,explnlp2_c_eval_apply_op
	,NULL
	,NULL
};

int RTOp_TOp_explnlp2_c_eval_construct( struct RTOp_RTOp* op )
{
	op->obj_data  = NULL;
	op->vtbl      = &RTOp_TOp_explnlp2_c_eval_vtbl;
	return 0; // success?
}

int RTOp_TOp_explnlp2_c_eval_destroy( struct RTOp_RTOp* op )
{
	op->obj_data  = NULL;
	op->vtbl      = NULL;
	return 0; // success?
}

// Implementation for RTOp_TOp_explnlp2_calc_py_D

static int explnlp2_calc_py_D_apply_op(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
	, const int num_vecs, const struct RTOp_SubVector vecs[]
	, const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
	, RTOp_ReductTarget targ_obj )
{
	size_t                 sub_dim;
	// xD
	const RTOp_value_type  *xD_val;
	ptrdiff_t              xD_val_s;
	// xI
	const RTOp_value_type  *xI_val;
	ptrdiff_t              xI_val_s;
	// c
	const RTOp_value_type  *c_val;
	ptrdiff_t              c_val_s;
	// d
	RTOp_value_type        *d_val;
	ptrdiff_t              d_val_s;
	// py
	RTOp_value_type        *py_val;
	ptrdiff_t              py_val_s;

	register RTOp_index_type  k;
	int                       all_unit_stride = 0;
	RTOp_value_type           denom;

	// task
	int task;
	assert(obj_data);
	task = *(int*)obj_data;
	assert(0 <= task && task <= 2);

	//
	// Validate the input
	//
	if( ( (task == 0 || task == 1) && num_vecs != 2 )
		|| ( (task == 2) && num_vecs != 3 )
		|| vecs == NULL )
		return RTOp_ERR_INVALID_NUM_VECS;
	if( ( (task == 0 || task == 1) && num_targ_vecs != 1 )
		|| ( (task == 2) && num_targ_vecs != 2 )
		|| targ_vecs == NULL )
		return RTOp_ERR_INVALID_NUM_TARG_VECS;
	if( targ_vecs[0].sub_dim != vecs[0].sub_dim
		|| targ_vecs[0].sub_dim != vecs[1].sub_dim
		|| ( task == 2 && ( targ_vecs[0].sub_dim != vecs[2].sub_dim ) )
		|| ( task == 2 && ( targ_vecs[0].sub_dim != targ_vecs[1].sub_dim ) )
		|| targ_vecs[0].global_offset != vecs[0].global_offset
		|| targ_vecs[0].global_offset != vecs[1].global_offset
		|| ( task == 2 && (targ_vecs[0].global_offset != vecs[2].global_offset ) )
		|| ( task == 2 && ( targ_vecs[0].global_offset != targ_vecs[1].global_offset ) ) )
		return RTOp_ERR_INCOMPATIBLE_VECS;
	if( vecs[0].indices || vecs[0].sub_nz == 0 || vecs[1].indices || vecs[1].sub_nz == 0
		|| ( task == 2 && (vecs[2].indices || vecs[2].sub_nz == 0) ) )
		return RTOp_ERR_UNSUPPORTED_VEC_TYPE;

	//
	// Get pointers to data
	//

	sub_dim = vecs[0].sub_dim;

	k = 0;
	// xD
	xD_val         = vecs[k].values;
	xD_val_s       = vecs[k].values_stride;
	++k;
	if( task == 1 || task == 2 ) {
		// xI
		xI_val         = vecs[k].values;
		xI_val_s       = vecs[k].values_stride;
		++k;
	}
	if( task == 0 || task == 2 ) {
		// c
		c_val         = vecs[k].values;
		c_val_s       = vecs[k].values_stride;
		++k;
	}
	k = 0;
	if( task == 1 || task == 2 ) {
		// d
		d_val         = targ_vecs[k].values;
		d_val_s       = targ_vecs[k].values_stride;
		++k;
	}
	if( task == 0 || task == 2 ) {
		// py
		py_val         = targ_vecs[k].values;
		py_val_s       = targ_vecs[k].values_stride;
		++k;
	}

	// Determine if all the vectors have unit stride!
	all_unit_stride = 1;
	for( k = 0; k < num_vecs && !all_unit_stride; ++k )
		if( vecs[k].values_stride != 1 )
			all_unit_stride = 0;
	for( k = 0; k < num_targ_vecs && !all_unit_stride; ++k )
		if( targ_vecs[k].values_stride != 1 )
			all_unit_stride = 0;

	//
	// Compute py and/or D
	//

	if( all_unit_stride) {
		if(task == 0) {
			// Compute py only
			for( k = 0; k < sub_dim; ++k )
				*py_val++ = *c_val++ / ( 1.0 - *xI_val++ );
		}
		if(task == 1) {
			// Compute D only
			for( k = 0; k < sub_dim; ++k )
				*d_val++ = ( *xD_val++ - 10.0 ) / ( 1.0 - *xI_val++ );
		}
		if(task == 2) {
			// Compute py and D
			for( k = 0; k < sub_dim; ++k ) {
				denom = ( 1.0 - *xI_val++ );
				*d_val++ = ( *xD_val++ - 10.0 ) / denom;
				*py_val++ = *c_val++ / denom;
			}
		}
	}
	else {
		assert(0); // ToDo: Implement if needed!
	}

	return 0; // success?
}

const char RTOp_TOp_explnlp2_calc_py_D_name[] = "TOp_explnlp2_calc_py_D";

const struct RTOp_RTOp_vtbl_t RTOp_TOp_explnlp2_calc_py_D_vtbl =
{
	&RTOp_obj_index_vtbl  // Use an index object (task) for instance data
	,&RTOp_obj_null_vtbl  // use null type for target object
	,NULL
	,explnlp2_calc_py_D_apply_op
	,NULL
	,NULL
};

int RTOp_TOp_explnlp2_calc_py_D_construct( int task, struct RTOp_RTOp* op )
{
	int result;
#ifdef RTOp_DEBUG
	assert( 0 <= task && task <= 2 );
#endif	
	op->obj_data  = NULL;
	op->vtbl      = &RTOp_TOp_explnlp2_calc_py_D_vtbl;
	result = op->vtbl->obj_data_vtbl->obj_create( NULL, NULL, &op->obj_data );
	if(result != 0) return result;
#ifdef RTOp_DEBUG
	assert(op->obj_data);
#endif
	*((int*)op->obj_data) = task;
	return 0; // success?
}

int RTOp_TOp_explnlp2_calc_py_D_set_task( int task, struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
	assert( 0 <= task && task <= 2 );
	assert(op->obj_data);
#endif
	*((int*)op->obj_data) = task;
	return 0; // success?
}

int RTOp_TOp_explnlp2_calc_py_D_destroy( struct RTOp_RTOp* op )
{
	int result;
	result = op->vtbl->reduct_vtbl->obj_free( NULL, NULL, &op->obj_data );
	if(result != 0) return result;
	op->obj_data  = NULL;
	op->vtbl      = NULL;
	return 0; // success?
}
