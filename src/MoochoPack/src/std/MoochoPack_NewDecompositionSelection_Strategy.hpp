// ////////////////////////////////////////////////////////////////////////////
// NewDecompositionSelection_Strategy.h
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

#ifndef NEW_DECOMPOSITION_SELECTION_STRATEGY_H
#define NEW_DECOMPOSITION_SELECTION_STRATEGY_H

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackTypes.h"
#include "GeneralIterationPack/include/Algorithm.h"

namespace ReducedSpaceSQPPack {

///
/** Abstract interface for an object that directs the selection of a new
 * decomposition.
 */
class NewDecompositionSelection_Strategy {
public:

	///
	virtual ~NewDecompositionSelection_Strategy() {}

	///
	virtual bool new_decomposition(
		rSQPAlgo& algo, Algorithm::poss_type step_poss
		,GeneralIterationPack::EDoStepType type, Algorithm::poss_type assoc_step_poss
		) = 0;

	///
	virtual void print_new_decomposition(
		const rSQPAlgo& algo, Algorithm::poss_type step_poss
		,GeneralIterationPack::EDoStepType type, Algorithm::poss_type assoc_step_poss
		,std::ostream& out, const std::string& leading_str
		) const = 0;

};	// end class NewDecompositionSelection_Strategy

}	// end namespace ReducedSpaceSQPPack 

#endif	// NEW_DECOMPOSITION_SELECTION_STRATEGY_H