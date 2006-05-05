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

#ifndef CHECK_DECOMPOSITION_FROM_R_PY_STEP_H
#define CHECK_DECOMPOSITION_FROM_R_PY_STEP_H

#include "MoochoPack_NewDecompositionSelection_Strategy.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

///
/** Check if the decomposition is going singular and if it is select a new decomposition.
 *
 * This steps checks if the decomposition is going singular if the computation
 * for the range space step looks like it is becomming more inaccurate.
 *
 * In particular we want to know how cond(C) is changing.  To do this we
 * will monitor increases in the error in solving the equation:
 \verbatim

   R*py + c(equ_decomp) = 0
 \endverbatim
 * Therefore we will check for increases in the ratio:
 \verbatim

   ratio = ||R*py + c(equ_decomp)||inf / ||c(equ_decomp)||inf
 \endverbatim
 *
 * If this ratio goes up dramatically, then this is a tail tell
 * sign that <tt>R</tt> is becomming illconditioned.  See the algorithm
 * printout for a more detailed description what what is going on.
 */
class CheckDecompositionFromRPy_Step
  : public IterationPack::AlgorithmStep // doxygen needs full name
{
public:

  /// <<std comp>> members for Decomposition Select Strategy object.
  STANDARD_COMPOSITION_MEMBERS( NewDecompositionSelection_Strategy, new_decomp_strategy )

  ///
  /** Set the maximum change in the relative error in the range space
     * step before the selection of a new decomposition is triggered.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_decomposition_cond_change_frac )

  ///
  CheckDecompositionFromRPy_Step(
     const new_decomp_strategy_ptr_t    &new_decomp_strategy
    ,value_type                         max_decomposition_cond_change_frac = 1e+4
    );

  /// Call the reset initialization of all defaults.
  void reset();

  /** @name Overridden from AlgorithmStep */
  //@{
  ///
  bool do_step(Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss);
  ///
  void print_step( const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
    , poss_type assoc_step_poss, std::ostream& out, const std::string& leading_str ) const;
  //@}

private:

  value_type	beta_min_;

  // Not defined and not to be called
  CheckDecompositionFromRPy_Step();

};	// end class CheckDecompositionFromRPy_Step

}	// end namespace MoochoPack 

#endif	// CHECK_DECOMPOSITION_FROM_R_PY_STEP_H
