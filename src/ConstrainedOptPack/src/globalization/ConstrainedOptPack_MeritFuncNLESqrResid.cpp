// //////////////////////////////////////////////////////////////////////////////////
// MeritFuncNLESqrResid.cpp
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

#include "MeritFuncNLESqrResid.hpp"
#include "DenseLinAlgPack/src/DVectorClass.hpp"
#include "DenseLinAlgPack/src/DVectorOp.hpp"

namespace ConstrainedOptPack {

MeritFuncNLESqrResid::MeritFuncNLESqrResid()
	: deriv_(0.0)
{}

value_type MeritFuncNLESqrResid::calc_deriv( const DVectorSlice& c_k )
{
	using DenseLinAlgPack::dot;
	return deriv_ = - dot(c_k,c_k);
}

// Overridden from MeritFuncNLP

value_type MeritFuncNLESqrResid::value(const DVectorSlice& c) const
{
	using DenseLinAlgPack::dot;
	return 0.5 * dot(c,c);
}

value_type MeritFuncNLESqrResid::deriv() const
{
	return deriv_;
}

void MeritFuncNLESqrResid::print_merit_func(std::ostream& out
	, const std::string& L ) const
{
	out
		<< L << "*** Define a square of constraint residuals merit funciton\n"
		<< L << "*** (assumes Gc_k'*d_k + c_k = 0):\n"
		<< L << "phi(c) = 1/2 * dot(c,c)\n"
		<< L << "Dphi(x_k,d_k) = - dot(c_k,c_k)\n";
}

}	// end namespace ConstrainedOptPack 
