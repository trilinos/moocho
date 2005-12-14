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

#ifndef ALGORITHM_TRACK_TESTING_H
#define ALGORITHM_TRACK_TESTING_H

#include "IterationPack_AlgorithmTracker.hpp"

namespace IterationPack {

///
/** Testing class.
  */
class AlgorithmTrackTesting : public AlgorithmTracker {
public:

	AlgorithmTrackTesting(const ostream_ptr_t& journal_out) : AlgorithmTracker(journal_out)
	{}

	// Overriden
	
	///
	void output_iteration(const Algorithm& algo) const;

	///
	void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;

};	// end class AlgorithmTrackTesting

}	// end namespace IterationPack 

#endif // ALGORITHM_TRACK_TESTING_H
