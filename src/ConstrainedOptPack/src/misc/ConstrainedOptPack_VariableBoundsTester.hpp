// ///////////////////////////////////////////////////////////////////////////
// VariableBoundsTester.h
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

#ifndef VARIABLE_BOUNDS_TESTER_H
#define VARIABLE_BOUNDS_TESTER_H

#include "ConstrainedOptimizationPackTypes.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ConstrainedOptimizationPack {

///
/** Tests that a set of variables are within their bounds.
  *
  \verbatim
    xL <= x <= xU
  \endverbatim
  *
  * The relative error for each comparison is
  * rel_err(i) = (xL(i)-x(i))/(1+||x||inf)
  * or rel_err(i) = (x(i)-xU(i))/(1+||x||inf).
  * If rel_err(i) >= error_tol, then the tests will be terminated immediately.
  * All of the rel_err(i) >= warning_tol will be printed.
  */
class VariableBoundsTester {
public:

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol )

	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol )

	///
	VariableBoundsTester(
		  value_type	warning_tol		= 1e-10
		, value_type	error_tol		= 1e-5
		);

	///
	virtual ~VariableBoundsTester() {}

	///
	/** Check that the variables are within bounds.
	 *
	 * @param  print_all_warnings
	 *              [in] If true, then all errors greater than warning_tol will
	 *              be printed.
	 * @param  xL   [in] Sparse lower bound vector (xL.size()==x.size())
	 * @param  xL_name
	 *              [in] The name of the vector xL (null terminated string).
	 * @param  xU   [in] Sparse upper bound vector (xU.size()==x.size())
	 * @param  xU_name
	 *              [in] The name of the vector xU (null terminated string).
	 * @param  x    [in] Variable to test that it is in bounds.
	 * @param  x_name
	 *              [in] The name of the vector x (null terminated string).
	 *
	 * @return \c true if all of the errors are greater than the error tolerances,
	 * otherwise it returns \c false
	 */
	virtual bool check_in_bounds(
		std::ostream* out, bool print_all_warnings, bool print_vectors
		,const VectorWithOp& xL, const char xL_name[]
		,const VectorWithOp& xU, const char xU_name[]
		,const VectorWithOp& x,  const char x_name[]
		);

};	// end class VariableBoundsTester

}	// end namespace ConstrainedOptimizationPackTypes

#endif	// VARIABLE_BOUNDS_TESTER_H
