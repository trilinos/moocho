/* @HEADER
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

/* ///////////////////////////////////////////// */
/* RTOp_TOp_max_abs_vec_scalar.h */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 7/13/2002 at 13:48 */
/* */

#ifndef RTOp_TOp_max_abs_vec_scalar_H
#define RTOp_TOp_max_abs_vec_scalar_H

#include "RTOp.h"
#include "RTOp_obj_null_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_TOp_max_abs_vec_scalar.h
 *
 \verbatim


element-wise transformation:
    z0 = max(fabs(z0),min_ele);

 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * ToDo: Write the documentation for this class!
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_TOp_max_abs_vec_scalar_vtbl;

/* Constructor */
int RTOp_TOp_max_abs_vec_scalar_construct( RTOp_value_type min_ele,  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_TOp_max_abs_vec_scalar_destroy( struct RTOp_RTOp* op );

/* Initialize the state of the operator object */
int RTOp_TOp_max_abs_vec_scalar_init( RTOp_value_type min_ele, struct RTOp_RTOp* op );



/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_TOp_max_abs_vec_scalar_H */