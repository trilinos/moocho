// //////////////////////////////////////////
// NLPWBCounterExample.cpp
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

#include <math.h>

#include "NLPWBCounterExample.hpp"
#include "DenseLinAlgPack/src/PermVecMat.hpp"
#include "DenseLinAlgPack/src/LinAlgOpPack.hpp"
#include "ThrowException.hpp"

namespace NLPInterfacePack {

// Constructors / initializers

NLPWBCounterExample::NLPWBCounterExample(
	value_type    a
	,value_type   b
	,value_type   x1_init
	,value_type   x2_init
	,value_type   x3_init
	,bool         nlp_selects_basis
	,bool         linear_obj
	)
	:is_initialized_(false)
	,nlp_selects_basis_(nlp_selects_basis)
	,basis_selection_was_given_(false)
	,linear_obj_(linear_obj)
	,n_orig_(3)
	,m_orig_(2)
	,a_(a)
	,b_(b)
{
#ifdef _DEBUG	
	const char msg_err_head[] = "NLPWBCounterExample::NLPWBCounterExample(...) : Error";
	THROW_EXCEPTION( b <= 0,       std::invalid_argument, msg_err_head<<"!" );
	THROW_EXCEPTION( a + b*b == 0, std::invalid_argument, msg_err_head<<"!" );
	THROW_EXCEPTION( x2_init < 0,  std::invalid_argument, msg_err_head<<"!" );
	THROW_EXCEPTION( x3_init < 0,  std::invalid_argument, msg_err_head<<"!" );
#endif
	// Set the number of nonzeros in the Jacobian of the equality constraints
	Gc_orig_nz_ = 4;
	// Resize the vectors for the initial guess and variable bounds
	xinit_orig_.resize(n_orig_);
	xl_orig_.resize(n_orig_);
	xu_orig_.resize(n_orig_);
	// Set the inital guess and the variable bounds
	const value_type inf = NLP::infinite_bound();
	xinit_orig_(1) = x1_init; xl_orig_(1) = -inf; xu_orig_(1) = +inf;
	xinit_orig_(2) = x2_init; xl_orig_(2) =  0.0; xu_orig_(2) = +inf;
	xinit_orig_(3) = x3_init; xl_orig_(3) =  0.0; xu_orig_(3) = +inf;
}

// Overridden public members from NLP

void NLPWBCounterExample::initialize(bool test_setup)
{
	if(is_initialized_) {
		NLPSerialPreprocessExplJac::initialize(test_setup);
		return;
	}
	// Nothing to initialize?
	NLPSerialPreprocessExplJac::initialize(test_setup);
	is_initialized_ = true;
}

bool NLPWBCounterExample::is_initialized() const
{
	return is_initialized_;
}

value_type NLPWBCounterExample::max_var_bounds_viol() const
{
	return +1e+20; // Functions defined everywhere!
}

// Overridden protected methods from NLPSerialPreprocess

bool NLPWBCounterExample::imp_nlp_has_changed() const
{
	return !is_initialized_;
}

size_type NLPWBCounterExample::imp_n_orig() const
{
	return n_orig_;
}

size_type NLPWBCounterExample::imp_m_orig() const
{
	return m_orig_;
}

size_type NLPWBCounterExample::imp_mI_orig() const
{
	return 0;
}

const DVectorSlice NLPWBCounterExample::imp_xinit_orig() const
{
	return xinit_orig_();
}

bool NLPWBCounterExample::imp_has_var_bounds() const
{
	return true;
}

const DVectorSlice NLPWBCounterExample::imp_xl_orig() const
{
	return xl_orig_();
}

const DVectorSlice NLPWBCounterExample::imp_xu_orig() const
{
	return xu_orig_();
}

const DVectorSlice NLPWBCounterExample::imp_hl_orig() const
{
	assert(0);  // Should never be called
 	return xinit_orig_();
}

const DVectorSlice NLPWBCounterExample::imp_hu_orig() const
{
	assert(0);  // Should never be called
 	return xinit_orig_();
}

void NLPWBCounterExample::imp_calc_f_orig(
	const DVectorSlice           &x_full
	,bool                        newx
	,const ZeroOrderInfoSerial   &zero_order_info
	) const
{
	DVectorSlice x = x_full(1,n_orig_);
	*zero_order_info.f = ( linear_obj_ ? x(1) : 0.5*x(1)*x(1) );
}

void NLPWBCounterExample::imp_calc_c_orig(
	const DVectorSlice           &x_full
	,bool                        newx
	,const ZeroOrderInfoSerial   &zero_order_info
	) const
{
	DVectorSlice x = x_full(1,n_orig_); 
	assert(zero_order_info.c);
	DVector      &c = *zero_order_info.c;
	//
	c(1) = x(1)*x(1) - x(2) + a_;
	c(2) = x(1)      - x(3) - b_;
}

void NLPWBCounterExample::imp_calc_h_orig(
	const DVectorSlice           &x_full
	,bool                        newx
	,const ZeroOrderInfoSerial   &zero_order_info
	) const
{
	assert(0); // Should never be called
}

void NLPWBCounterExample::imp_calc_Gf_orig(
	const DVectorSlice           &x_full
	,bool                        newx
	,const ObjGradInfoSerial     &obj_grad_info
	) const
{
	DVectorSlice  x     = x_full(1,n_orig_); 
	DVector       &Gf   = *obj_grad_info.Gf;
	Gf(1) = (linear_obj_ ? 1.0 : x(1) );
	Gf(2) = 0.0;
	Gf(3) = 0.0;
}

bool NLPWBCounterExample::imp_get_next_basis(
	IVector      *var_perm_full
	,IVector     *equ_perm_full
	,size_type   *rank_full
	,size_type   *rank
	)
{
#ifdef _DEBUG
	assert(var_perm_full);
	assert(equ_perm_full);
	assert(rank_full);
	assert(rank);
#endif
	if(basis_selection_was_given_)
		return false; // Already gave this basis selection.
	//
	// Select x(2) ans x(3) as the basic variables (sorted!)
	//
	var_perm_full->resize(n_orig_);
	equ_perm_full->resize(m_orig_);
	(*var_perm_full)(1) = 2;  // The basis variables
	(*var_perm_full)(2) = 3;  // ""
	(*var_perm_full)(3) = 1;  // The nonbasis variable
	DenseLinAlgPack::identity_perm( equ_perm_full ); // Gc_orig is full rank
	// Set the rank (full rank)
	*rank_full                 = m_orig_;
	*rank                      = m_orig_;
	basis_selection_was_given_ = true;
	return true;
}

void NLPWBCounterExample::imp_report_orig_final_solution(
	const DVectorSlice      &x_orig
	,const DVectorSlice     *lambda_orig
	,const DVectorSlice     *lambdaI_orig
	,const DVectorSlice     *nu_orig
	,bool                   is_optimal
	) const
{
	// ToDo: Do something with the final soltuion?
}

bool NLPWBCounterExample::nlp_selects_basis() const
{
	return nlp_selects_basis_;
}

// Overridden protected methods from NLPSerialPreprocessExplJac

size_type NLPWBCounterExample::imp_Gc_nz_orig() const
{
	return Gc_orig_nz_;
}

size_type NLPWBCounterExample::imp_Gh_nz_orig() const
{
	return 0;
}

void NLPWBCounterExample::imp_calc_Gc_orig(
	const DVectorSlice& x_full, bool newx
	, const FirstOrderExplInfo& first_order_expl_info
	) const
{
	DVectorSlice x = x_full(1,n_orig_); 
	// Get references/pointers to data for Gc to be computed/updated.
	index_type
		&Gc_nz = *first_order_expl_info.Gc_nz;
	value_type
		*Gc_v = &(*first_order_expl_info.Gc_val)[0];
	index_type
		*Gc_i = ( first_order_expl_info.Gc_ivect
					  ? &(*first_order_expl_info.Gc_ivect)[0] : NULL ),
		*Gc_j = ( first_order_expl_info.Gc_jvect
					  ? &(*first_order_expl_info.Gc_jvect)[0] : NULL );
	assert( (Gc_i != NULL) == (Gc_j != NULL) );
	// Set up the nonzero structure of Gc_orig (sorted by variable and then by constraint)
	size_type nz = 0;
	if( Gc_i ) {
		// c(1)                             // c(2)
		Gc_i[nz] = 1; Gc_j[nz] = 1; ++nz;   Gc_i[nz] = 1; Gc_j[nz] = 2; ++nz;  // x(1)
		Gc_i[nz] = 2; Gc_j[nz] = 1; ++nz;                                      // x(2)
		                                    Gc_i[nz] = 3; Gc_j[nz] = 2; ++nz;  // x(3)
	
	    assert(nz == Gc_orig_nz_);
	}
	// Fill in the nonzero values of Gc_orig (must have the same order as structure!)
	nz = 0;
	// c(1)                    // c(2)
	Gc_v[nz] = 2*x(1); ++nz;   Gc_v[nz] = +1.0; ++nz;  // x(1)
	Gc_v[nz] =   -1.0; ++nz;                           // x(2)
	                           Gc_v[nz] = -1.0; ++nz;  // x(3)
    assert(nz == Gc_orig_nz_);
}

void NLPWBCounterExample::imp_calc_Gh_orig(
	const DVectorSlice& x_full, bool newx
	, const FirstOrderExplInfo& first_order_expl_info
	) const
{
	assert(0); // Should never be called
}

}	// end namespace NLPInterfacePack
