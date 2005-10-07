/*
// @HEADER
// ***********************************************************************
// 
//      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 6/27/2002 at 20:41 */
/* */

#define max(a,b) ( (a) > (b) ? (a) : (b) )
#define min(a,b) ( (a) < (b) ? (a) : (b) )

#include "RTOp_TOp_inv_of_difference.h"
#include "RTOp_obj_value_vtbl.h"  /* vtbl for operator object instance data */



/* Implementation functions for RTOp_RTOp */

static int RTOp_TOp_inv_of_difference_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget reduct_obj )
{
  /* */
  /* Declare local variables */
  /* */

    /* Access to the operator object instance data */
    RTOp_value_type *alpha = (RTOp_value_type*)obj_data;
    /* Vector data */
    RTOp_index_type           sub_dim;
    /* z0 */
    RTOp_value_type           *z0_val;
    ptrdiff_t                 z0_val_s;
    /* v0 */
    const RTOp_value_type     *v0_val;
    ptrdiff_t                 v0_val_s;
    /* v1 */
    const RTOp_value_type     *v1_val;
    ptrdiff_t                 v1_val_s;

    register RTOp_index_type  k;

  /* */
  /* Validate the input */
  /* */
    if( num_vecs != 2 || ( num_vecs && vecs == NULL ) )
        return RTOp_ERR_INVALID_NUM_VECS;
    if( num_targ_vecs != 1 || ( num_targ_vecs && targ_vecs == NULL ) )
        return RTOp_ERR_INVALID_NUM_TARG_VECS;
    if( /* Validate sub_dim */
        vecs[1].sub_dim != vecs[0].sub_dim
        || targ_vecs[0].sub_dim != vecs[0].sub_dim
        )
        return RTOp_ERR_INCOMPATIBLE_VECS;
    assert(obj_data);


  /* */
  /* Get pointers to data */
  /* */
    sub_dim       = vecs[0].sub_dim;
    /* z0 */
    z0_val        = targ_vecs[0].values;
    z0_val_s      = targ_vecs[0].values_stride;
    /* v0 */
    v0_val        = vecs[0].values;
    v0_val_s      = vecs[0].values_stride;
    /* v1 */
    v1_val        = vecs[1].values;
    v1_val_s      = vecs[1].values_stride;


  /* */
  /* Apply the operator: */
  /* */
    /*    element-wise transformation : z0 = alpha/(v0-v1); */
    /* */
    for( k = 0; k < sub_dim; ++k, v0_val += v0_val_s, v1_val += v1_val_s, z0_val += z0_val_s )
    {
        /* Element-wise transformation */
/*        if (fabs((*v0_val)-(*v1_val)) > 5.562685e-309) */
    { (*z0_val) = (*alpha)/((*v0_val)-(*v1_val)); }
/*  else */
/*    { (*z0_val) = 1.797693e+308; } */
    }

  return 0; /* success? */
}

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_TOp_inv_of_difference_vtbl =
{
  &RTOp_obj_value_vtbl
  ,&RTOp_obj_null_vtbl
  ,"TOp_inv_of_difference"
  ,NULL
  ,RTOp_TOp_inv_of_difference_apply_op
  ,NULL
  ,NULL
};

/* Class specific functions */

int RTOp_TOp_inv_of_difference_construct( RTOp_value_type alpha,  struct RTOp_RTOp* op )
{
#ifdef RTOp_DEBUG
  assert(op);
#endif
  op->obj_data  = NULL;
  op->vtbl      = &RTOp_TOp_inv_of_difference_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  return RTOp_TOp_inv_of_difference_init(alpha,op);
}

int RTOp_TOp_inv_of_difference_destroy( struct RTOp_RTOp* op )
{
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->obj_data  = NULL;
  op->vtbl      = NULL;
  return 0;
}

int RTOp_TOp_inv_of_difference_init( RTOp_value_type alpha, struct RTOp_RTOp* op )
{
    RTOp_value_type *ptr_alpha = (RTOp_value_type*)op->obj_data;
    *ptr_alpha = alpha;
    return 0;
}
