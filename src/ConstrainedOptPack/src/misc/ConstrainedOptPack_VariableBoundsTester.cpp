// //////////////////////////////////////////////////////////////////////////////
// VariableBoundsTester.cpp
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

#include "ConstrainedOptimizationPack/src/VariableBoundsTester.hpp"
#include "AbstractLinAlgPack/src/VectorSpace.hpp"
#include "AbstractLinAlgPack/src/VectorMutable.hpp"
#include "AbstractLinAlgPack/src/VectorAuxiliaryOps.hpp"

namespace ConstrainedOptimizationPack {

// public

VariableBoundsTester::VariableBoundsTester(
	value_type   warning_tol
	,value_type  error_tol
	)
	:warning_tol_(warning_tol)
	,error_tol_(error_tol)
{}

bool VariableBoundsTester::check_in_bounds(
	std::ostream* out, bool print_all_warnings, bool print_vectors
	,const Vector& xL, const char xL_name[]
	,const Vector& xU, const char xU_name[]
	,const Vector& x,  const char x_name[]
	)
{
	using AbstractLinAlgPack::max_near_feas_step;

	if(out)
		*out
			<< "\n*** Checking that variables are in bounds\n";

	VectorSpace::vec_mut_ptr_t zero = x.space().create_member(0.0);
	std::pair<value_type,value_type>
		u = max_near_feas_step( x, *zero, xL, xU, warning_tol() );
	if(u.first < 0.0) {
		if(out)
			*out << "\nWarning! the variables " << xL_name << " <= " << x_name << " <= " << xU_name
				<< " are out of bounds by more than warning_tol = "	<< warning_tol() << "\n";
		u = max_near_feas_step( x, *zero, xL, xU, error_tol() );
		if(u.first < 0.0) {
			if(out)
				*out << "\nError! the variables " << xL_name << " <= " << x_name << " <= " << xU_name
					<< " are out of bounds by more than error_tol = " << error_tol() << "\n";
			return false;
		}
	}
	return true;
}

}	// end namespace ConstrainedOptimizationPack
