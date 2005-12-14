/*
// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

/* */
/* RAB: 1/22/01: Warning, this implementation based */
/* on the struct may not be portable in a heterogeneous */
/* environment!  The user has been warned! */
/* */

#include <stdlib.h>

#include "RTOp_TOp_random_vector.h"
#include "RTOp_obj_null_vtbl.h"
#include "RTOp_obj_free_free.h"

/* Functions for the reduction/transformation object */

struct RTOp_TOp_random_vector_bnd_t {
  RTOp_value_type  l;
  RTOp_value_type  u;
};

static int get_op_type_num_entries(
  const struct RTOp_obj_type_vtbl_t*   vtbl
  ,const void* obj_data
  ,int* num_values
  ,int* num_indexes
  ,int* num_chars
  )
{
  *num_values  = 2;
  *num_indexes = 0;
  *num_chars   = 0;
  return 0;
}

static int op_create(
  const struct RTOp_obj_type_vtbl_t* vtbl, const void* instance_data
  , RTOp_ReductTarget* obj )
{
	*obj = malloc(sizeof(struct RTOp_TOp_random_vector_bnd_t));
	return 0;
}

static int extract_op_state(
  const struct RTOp_obj_type_vtbl_t*   vtbl
  ,const void *       dummy
  ,void *             obj_data
  ,int                num_values
  ,RTOp_value_type    value_data[]
  ,int                num_indexes
  ,RTOp_index_type    index_data[]
  ,int                num_chars
  ,RTOp_char_type     char_data[]
  )
{
  const struct RTOp_TOp_random_vector_bnd_t *bnd = NULL;
  assert(obj_data);
  assert( num_values  == 2 );
  assert( num_indexes == 0 );
  assert( num_chars   == 0 );
  bnd = (const struct RTOp_TOp_random_vector_bnd_t*)obj_data;
  value_data[0] = bnd->l;
  value_data[1] = bnd->u;
  return 0;
}

static int load_op_state(
  const struct RTOp_obj_type_vtbl_t*   vtbl
  ,const void *            dummy
  ,int                     num_values
  ,const RTOp_value_type   value_data[]
  ,int                     num_indexes
  ,const RTOp_index_type   index_data[]
  ,int                     num_chars
  ,const RTOp_char_type    char_data[]
  ,void **                 obj_data
  )
{
  struct RTOp_TOp_random_vector_bnd_t *bnd = NULL;
  assert( obj_data );
  assert( num_values  == 2 );
  assert( num_indexes == 0 );
  assert( num_chars   == 0 );
  if(*obj_data == NULL)
    *obj_data = malloc(sizeof(struct RTOp_TOp_random_vector_bnd_t));
  bnd = (struct RTOp_TOp_random_vector_bnd_t*)*obj_data;
  bnd->l = value_data[0];
  bnd->u = value_data[1];
  return 0;
}

static struct RTOp_obj_type_vtbl_t  instance_obj_vtbl =
{
  get_op_type_num_entries
  ,op_create
  ,NULL
  ,RTOp_obj_free_free
  ,extract_op_state
  ,load_op_state
};

/* Implementation functions for RTOp_RTOp */

static int RTOp_TOp_random_vector_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  /* */
  /* Get pointers to data */
  /* */

  /* lbnd, ubnd */
  const struct RTOp_TOp_random_vector_bnd_t *bnd = (const struct RTOp_TOp_random_vector_bnd_t*)obj_data;

  /* z */
  RTOp_index_type  z_sub_dim;
  RTOp_value_type  *z_val      = NULL;
  ptrdiff_t        z_val_s;

  register RTOp_index_type k;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 0 || vecs != NULL )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 1 || targ_vecs == NULL )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;

  /* z */
  z_sub_dim     = targ_vecs[0].sub_dim;
  z_val         = targ_vecs[0].values;
  z_val_s       = targ_vecs[0].values_stride;

  /* */
  /* Assign the elements to random values */
  /* */

  for( k = 0; k < z_sub_dim; ++k, z_val += z_val_s )
    *z_val = bnd->l + ((RTOp_value_type)rand())/RAND_MAX * (bnd->u - bnd->l);

  return 0; /* success? */
}

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_TOp_random_vector_vtbl =
{
  &instance_obj_vtbl
  ,&RTOp_obj_null_vtbl /* use null type for target object */
  ,"TOp_random_vector"
  ,NULL /* use default from reduct_vtbl */
  ,RTOp_TOp_random_vector_apply_op
  ,NULL
  ,NULL
};

/* Class specific functions */

int RTOp_TOp_random_vector_construct( RTOp_value_type lbnd, RTOp_value_type ubnd
  , struct RTOp_RTOp* op )
{
  struct RTOp_TOp_random_vector_bnd_t *bnd = NULL;
  op->vtbl      = &RTOp_TOp_random_vector_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  bnd = (struct RTOp_TOp_random_vector_bnd_t*)op->obj_data;
  bnd->l = lbnd;
  bnd->u = ubnd;
  return 0; /* success? */
}

int RTOp_TOp_random_vector_destroy( struct RTOp_RTOp* op )
{
  free( op->obj_data );
  op->obj_data  = NULL;
  op->vtbl      = NULL;
  return 0; /* success? */
}

int RTOp_TOp_random_vector_set_bounds( RTOp_value_type lbnd, RTOp_value_type ubnd
  , struct RTOp_RTOp* op )
{
  struct RTOp_TOp_random_vector_bnd_t *bnd = (struct RTOp_TOp_random_vector_bnd_t*)op->obj_data;
  bnd->l = lbnd;
  bnd->u = ubnd;
  return 0; /* success? */
}
