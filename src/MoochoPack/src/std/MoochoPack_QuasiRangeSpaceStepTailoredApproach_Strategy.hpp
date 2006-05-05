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

#ifndef QUASI_RANGE_SPACE_STEP_TAILORED_APPROACH_STRATEGY_H
#define QUASI_RANGE_SPACE_STEP_TAILORED_APPROACH_STRATEGY_H

#include "MoochoPack_QuasiRangeSpaceStep_Strategy.hpp"

namespace MoochoPack {

///
/** Strategy class for computing a quasi range space step for the
 * tailored approach NLP interface.
 */
class QuasiRangeSpaceStepTailoredApproach_Strategy : public QuasiRangeSpaceStep_Strategy {
public:

  // ////////////////////////////////////////////////////////////
  // Overridden from QuasiRangeSpaceStep_Strategy

  ///
  /** Calls the NLPrSQPTailoredApproach iterface to compute the step.
   *
   * ToDo: Finish documentation!
   */
    bool solve_quasi_range_space_step(
    std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
    ,const DVectorSlice& xo, const DVectorSlice& c_xo, DVectorSlice* v
      );

  ///
  void print_step( std::ostream& out, const std::string& leading_str ) const;

}; // end class QuasiRangeSpaceStepTailoredApproach_Strategy

} // end namespace MoochoPack

#endif // QUASI_RANGE_SPACE_STEP_TAILORED_APPROACH_STRATEGY_H
