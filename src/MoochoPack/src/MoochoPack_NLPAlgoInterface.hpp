// ////////////////////////////////////////////////////////////////////////////
// rSQPAlgoInterface.hpp
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

#ifndef RSQP_ALGO_INTERFACE_H
#define RSQP_ALGO_INTERFACE_H

#include "ReducedSpaceSQPPackTypes.hpp"
#include "rSQPSolverClientInterface.hpp"

namespace ReducedSpaceSQPPack {

///
/** Interface \c rSQPAlgoContainer uses to access \c rSQPAlgo.
 *
 * This interface helps avoid dangerous usage stategies for an \c rSQPAlgo object.
 */
class rSQPAlgoInterface {
public:

	///
	virtual ~rSQPAlgoInterface() {}

	///
	/** Print the algorithm description.
	 */
	virtual void interface_print_algorithm(std::ostream& out) const = 0;

	///
	/** Start the iterations.
	  *
	  * This function returns true if the solution was found and false
	  * if the maximum number of iterations was reached before the
	  * solution was found.
	  */
	virtual rSQPSolverClientInterface::EFindMinReturn dispatch() = 0;

	///
	/** Return the state object.
	 */
	virtual const rSQPState& retrieve_state() const = 0;

	/** @name Algorithm timing */
	//@{

	///
	virtual void interface_set_algo_timing( bool algo_timing ) = 0;
	///
	virtual bool interface_algo_timing() const = 0;
	///
	virtual void interface_print_algorithm_times( std::ostream& out ) const = 0;

	//@}

};	// end class rSQPAlgoInterface

}	// end namespace ReducedSpaceSQPPack

#endif	// RSQP_ALGO_INTERFACE_H
