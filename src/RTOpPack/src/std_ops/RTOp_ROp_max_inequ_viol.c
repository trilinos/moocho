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

#define max(a,b) ( (a) > (b) ? (a) : (b) )
#define min(a,b) ( (a) < (b) ? (a) : (b) )

#include "RTOp_ROp_max_inequ_viol.h"
#include "RTOp_obj_null_vtbl.h"
#include "RTOp_obj_value_vtbl.h"
#include "RTOp_obj_free_free.h"
#include "RTOp_get_reduct_op.hpp"

/* */
/* Implementation functions */
/* */

/* Functions for the reduction target object */

static int get_targ_type_num_entries(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void* obj_data
  ,int* num_values
  ,int* num_indexes
  ,int* num_chars
  )
{
  *num_values  = 3;
  *num_indexes = 2;
  *num_chars   = 0;
  return 0;
}

static int targ_obj_reinit(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* obj_data
  , RTOp_ReductTarget targ_obj )
{
  struct RTOp_ROp_max_inequ_viol_reduct_obj_t
    *targ = (struct RTOp_ROp_max_inequ_viol_reduct_obj_t*)targ_obj;
  targ->max_viol   = 0.0; /* None violated yet! */
  targ->v_i        = 0.0; /* arbitrary */
  targ->vLU_i      = 0.0; /* arbitrary */
  targ->max_viol_i = 0;   /* None violated yet! */
  targ->bnd_type   = -2;  /* Invalid type (just in case used by accident) */
  return 0;
}

static int targ_obj_create(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* obj_data
  ,RTOp_ReductTarget* targ_obj
  )
{
  const int mem_size = sizeof(struct RTOp_ROp_max_inequ_viol_reduct_obj_t);
  *targ_obj = malloc( mem_size );
  return targ_obj_reinit(vtbl,obj_data,*targ_obj);
}

static int targ_extract_state(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void *       obj_data
  ,void *             reduct_obj
  ,int                num_values
  ,RTOp_value_type    value_data[]
  ,int                num_indexes
  ,RTOp_index_type    index_data[]
  ,int                num_chars
  ,RTOp_char_type     char_data[]
  )
{
  struct RTOp_ROp_max_inequ_viol_reduct_obj_t *targ = NULL;
  assert( reduct_obj );
  assert( num_values  == 3 );
  assert( num_indexes == 3 );
  assert( num_chars   == 0 );
  targ = (struct RTOp_ROp_max_inequ_viol_reduct_obj_t*)reduct_obj;
  value_data[0] = targ->max_viol;
  value_data[1] = targ->v_i;
  value_data[2] = targ->vLU_i;
  index_data[0] = targ->max_viol_i;
  index_data[1] = targ->bnd_type;
  return 0;
}

static int targ_load_state(
  const struct RTOp_obj_type_vtbl_t* vtbl
  ,const void*             obj_data
  ,int                     num_values
  ,const RTOp_value_type   value_data[]
  ,int                     num_indexes
  ,const RTOp_index_type   index_data[]
  ,int                     num_chars
  ,const RTOp_char_type    char_data[]
  ,void **                 reduct_obj
  )
{
  struct RTOp_ROp_max_inequ_viol_reduct_obj_t *targ = NULL;
  assert( *reduct_obj );
  assert( num_values  == 3 );
  assert( num_indexes == 2 );
  assert( num_chars   == 0 );
  targ = (struct RTOp_ROp_max_inequ_viol_reduct_obj_t*)*reduct_obj;
  targ->max_viol   = value_data[0];
  targ->v_i        = value_data[1];
  targ->vLU_i      = value_data[2];
  targ->max_viol_i = index_data[0];
  targ->bnd_type   = index_data[1];
  return 0;
}

static const struct RTOp_obj_type_vtbl_t targ_obj_vtbl =
{
  get_targ_type_num_entries
  ,targ_obj_create
  ,targ_obj_reinit
  ,RTOp_obj_free_free
  ,targ_extract_state
  ,targ_load_state
};

/* Other functions */

static int RTOp_ROp_max_inequ_viol_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  ,const int num_vecs, const struct RTOp_SubVector vecs[]
  ,const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  ,RTOp_ReductTarget targ_obj
  )
{
  /* */
  /* Declare local variables */
  /* */

  /* targ */
  struct RTOp_ROp_max_inequ_viol_reduct_obj_t
    *targ = NULL;
  /* global_off */
  size_t                 global_offset;
  /* sub_dim */
  size_t                 sub_dim;
  /* v */
  const RTOp_value_type  *v_val = NULL;
  ptrdiff_t              v_val_s;
  /* vL */
  const RTOp_value_type  *vL_val = NULL;
  ptrdiff_t              vL_val_s;
  /* vU */
  const RTOp_value_type  *vU_val = NULL;
  ptrdiff_t              vU_val_s;

  register size_t  k;
  RTOp_index_type  i;
  RTOp_value_type  v_scale;
  RTOp_value_type  violL;
  RTOp_value_type  violU;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 3 )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 0 )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  if(    vecs[0].global_offset != vecs[1].global_offset
    || vecs[0].sub_dim       != vecs[1].sub_dim
    || vecs[0].global_offset != vecs[2].global_offset
    || vecs[0].sub_dim       != vecs[2].sub_dim )
    return RTOp_ERR_INCOMPATIBLE_VECS;

  /* */
  /* Get pointers to the data */
  /* */

  /* targ */
  targ            = (struct RTOp_ROp_max_inequ_viol_reduct_obj_t*)targ_obj;
  /* global_off */
  global_offset   = vecs[0].global_offset;
  /* sub_dim */
  sub_dim         = vecs[0].sub_dim;
  /* v */
  v_val           = vecs[0].values;
  v_val_s         = vecs[0].values_stride;
  /* vL */
  vL_val          = vecs[1].values;
  vL_val_s        = vecs[1].values_stride;
  /* vU */
  vU_val          = vecs[2].values;
  vU_val_s        = vecs[2].values_stride;

  /* */
  /* Perform the reduction operation. */
  /* */

  i = global_offset + 1;
  for( k = 0; k < sub_dim; ++k, ++i, v_val += v_val_s, vL_val += vL_val_s, vU_val += vU_val_s ) {
    v_scale = 1.0 / (1.0 + fabs(*v_val));
    /* (vL - v)*v_scale */
    violL = (*vL_val - *v_val) * v_scale;
    /* (v - vU)*v_scale */
    violU = (*v_val - *vU_val) * v_scale;
    /* Perform the reduction */
    if(
      ( max(violL,violU) > targ->max_viol )
      ||
      ( max(violL,violU) == targ->max_viol && i < targ->max_viol_i )
      )
     {
      targ->bnd_type = ( violL > 0.0
                ? ( *vL_val == *vU_val
                  ? 0  /* EQUALITY */
                  : -1 /* LOWER */
                  )
                : +1 /* UPPER */
                );
      targ->max_viol   = ( targ->bnd_type <= 0 ? violL : violU );
      targ->v_i        = *v_val;
      targ->vLU_i      = ( targ->bnd_type <= 0 ? *vL_val : *vU_val );
      targ->max_viol_i = i;
    }
  }

  return 0; /* success? */
}

static int reduce_reduct_objs(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data /* Can be NULL! */
  , RTOp_ReductTarget in_reduct_obj, RTOp_ReductTarget inout_reduct_obj )
{
  const struct RTOp_ROp_max_inequ_viol_reduct_obj_t
    *i_targ = (const struct RTOp_ROp_max_inequ_viol_reduct_obj_t*)in_reduct_obj;
  struct RTOp_ROp_max_inequ_viol_reduct_obj_t
    *io_targ = (struct RTOp_ROp_max_inequ_viol_reduct_obj_t*)inout_reduct_obj;

  if(
    ( i_targ->max_viol > io_targ->max_viol )
    ||
    ( i_targ->max_viol == io_targ->max_viol && i_targ->max_viol_i < io_targ->max_viol_i )
    )
  {
    io_targ->max_viol     = i_targ->max_viol;
    io_targ->v_i          = i_targ->v_i;
    io_targ->vLU_i        = i_targ->vLU_i;
    io_targ->max_viol_i   = i_targ->max_viol_i;
    io_targ->bnd_type     = i_targ->bnd_type;
  }
  return 0;
}

INSERT_GET_REDUCT_OP_FUNCS(
  3,2,0,RTOp_ROp_max_inequ_viol_reduct_obj_t,reduce_reduct_objs
  ,targ_load_state,targ_extract_state
  ,external_reduct_op,get_reduct_op)

const struct RTOp_RTOp_vtbl_t RTOp_ROp_max_inequ_viol_vtbl =
{
  &RTOp_obj_null_vtbl
  ,&targ_obj_vtbl
  ,"ROp_max_inequ_viol"
  ,NULL
  ,RTOp_ROp_max_inequ_viol_apply_op
  ,reduce_reduct_objs
  ,get_reduct_op
};

/* Class specific functions */

int RTOp_ROp_max_inequ_viol_construct( struct RTOp_RTOp* op )
{
  op->vtbl     = &RTOp_ROp_max_inequ_viol_vtbl;
  op->obj_data = NULL;
  return 0; /* success? */
}

int RTOp_ROp_max_inequ_viol_destroy( struct RTOp_RTOp* op )
{
  op->vtbl     = NULL;
  op->obj_data = NULL;
  return 0; /* success? */
}

struct RTOp_ROp_max_inequ_viol_reduct_obj_t
RTOp_ROp_max_inequ_viol_val(RTOp_ReductTarget targ_obj)
{
  return *(struct RTOp_ROp_max_inequ_viol_reduct_obj_t*)targ_obj;
}
