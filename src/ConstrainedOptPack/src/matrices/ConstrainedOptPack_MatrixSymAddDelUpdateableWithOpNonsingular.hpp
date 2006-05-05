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

#ifndef MATRIX_SYM_ADD_DEL_UPDATEABLE_WITH_OP_NONSINGULAR_H
#define MATRIX_SYM_ADD_DEL_UPDATEABLE_WITH_OP_NONSINGULAR_H

#include "ConstrainedOptPack_Types.hpp"

namespace ConstrainedOptPack {

///
/** Interface for updating a symmetric matrix and its factorization 
 * by adding and deleting rows and columns and preforming operations with it.
 *
 * ToDo: Finish documentation.
 */
class MatrixSymAddDelUpdateableWithOpNonsingular {
public:

  ///
  virtual ~MatrixSymAddDelUpdateableWithOpNonsingular() {}
  ///
  virtual const MatrixSymOpNonsing& op_interface() const = 0;
  ///
  virtual MatrixSymAddDelUpdateable& update_interface() = 0;
  ///
  virtual const MatrixSymAddDelUpdateable& update_interface() const = 0;

};	// end class MatrixSymAddDelUpdateableWithOpNonsingular

}	// namespace ConstrainedOptPack 

#endif	// MATRIX_SYM_ADD_DEL_UPDATEABLE_WITH_OP_NONSINGULAR_H
