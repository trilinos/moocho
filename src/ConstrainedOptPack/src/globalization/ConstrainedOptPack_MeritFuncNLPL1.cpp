// /////////////////////////////////////////////////////////////////////////////
// MeritFuncNLPL1.cpp
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

#include "ConstrainedOptPack/src/globalization/MeritFuncNLPL1.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/Vector.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorStdOps.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/VectorAuxiliaryOps.hpp"
#include "ThrowException.hpp"
#include "dynamic_cast_verbose.hpp"

namespace ConstrainedOptPack {

MeritFuncNLPL1::MeritFuncNLPL1()
	: deriv_(0.0), mu_(0.0)
{}

// Overridden from MeritFuncNLP

MeritFuncNLP& MeritFuncNLPL1::operator=(const MeritFuncNLP& merit_func)
{
	using DynamicCastHelperPack::dyn_cast;
	const MeritFuncNLPL1 &merit_func_l1 = dyn_cast<const MeritFuncNLPL1>(merit_func);
	if(this == &merit_func_l1)
		return *this; // assignment to self
	this->deriv_   = merit_func_l1.deriv_;
	this->mu_      = merit_func_l1.mu_;
	return *this;
}

value_type MeritFuncNLPL1::value(
	value_type       f
	,const Vector    *c
	,const Vector    *h
	,const Vector    *hl
	,const Vector    *hu
	) const
{
	value_type phi_val_h = 0.0;
	if(h) {
		size_type   max_viol_i;
		value_type  max_viol, h_i, hlu_i;
		int         bnd_type;
		AbstractLinAlgPack::max_inequ_viol(*h,*hl,*hu,&max_viol_i,&max_viol,&h_i,&bnd_type,&hlu_i);
		if(max_viol_i) phi_val_h += mu_ * fabs(h_i - hlu_i);
	}
	return f + ( c ? mu_ * c->norm_1() : 0.0) + phi_val_h;
}

value_type MeritFuncNLPL1::deriv() const
{
	return deriv_;
}

void MeritFuncNLPL1::print_merit_func(std::ostream& out
	, const std::string& L ) const
{
	out
		<< L << "*** Define L1 merit funciton (assumes Gc_k'*d_k + c_k = 0):\n"
		<< L << "phi(f,c) = f + mu_k*( norm(c,1) + max_viol( hl <= h <= hu ) )\n"
		<< L << "Dphi(x_k,d_k) = Gf_k' * d_k - mu*( norm(c_k,1)  + max_viol( hl <= h_k <= hu ) )\n";
}

// Overridden from MeritFuncNLPDirecDeriv

value_type MeritFuncNLPL1::calc_deriv(
	const Vector    &Gf_k
	,const Vector   *c_k
	,const Vector   *h_k
	,const Vector   *hl
	,const Vector   *hu
	,const Vector   &d_k
	)
{
	using AbstractLinAlgPack::dot;
	THROW_EXCEPTION(
		h_k || hl || hu, std::logic_error
		,"MeritFuncNLPL1::value(...) : Error! general inequalities are not supported yet" );
	return deriv_ = dot( Gf_k, d_k ) - ( c_k ? mu_ * c_k->norm_1() : 0.0 );
}

// Overridden from MeritFuncPenaltyParam

void MeritFuncNLPL1::mu(value_type mu)
{
	mu_ = mu;
}

value_type MeritFuncNLPL1::mu() const
{
	return mu_;
}

}	// end namespace ConstrainedOptPack 
