// //////////////////////////////////////////
// ExampleNLPBanded.cpp
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

#include "ExampleNLPBanded.h"
#include "LinAlgPack/include/PermVecMat.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "ThrowException.h"

namespace NLPInterfacePack {

// Constructors / initializers

ExampleNLPBanded::ExampleNLPBanded(
	size_type     nD
	,size_type    nI
	,size_type    bw
	,size_type    mU
	,size_type    mI
	,value_type   xo
	,value_type   xl
	,value_type   xu
	,value_type   hl
	,value_type   hu
	,bool         nlp_selects_basis
	,value_type   diag_scal
	,value_type   diag_vary
	,bool         sym_basis
	)
	:is_initialized_(false)
	,nlp_selects_basis_(nlp_selects_basis)
	,basis_selection_was_given_(false)
	,multi_calc_(false)
	,nD_(nD)
	,nI_(nI)
	,bw_(bw)
	,mU_(mU)
	,mI_(mI)
	,diag_scal_(diag_scal)
	,diag_vary_(diag_vary)
	,fu_( sym_basis ? 3 : 6 )
{
#ifdef _DEBUG	
	const char msg_err_head[] = "ExampleNLPBanded::ExampleNLPBanded(...) : Error";
	THROW_EXCEPTION(
		nI > nD, std::invalid_argument
		,msg_err_head<<"!" );
	THROW_EXCEPTION(
		bw < 1 || nD < bw, std::invalid_argument
		,msg_err_head<<"!" );
	THROW_EXCEPTION(
		mU != 0, std::invalid_argument
		,msg_err_head<<", can't handle undecomposed equalities yet!" );
#endif
	Gc_orig_nz_ = nD_ * ( 1 + 2*(bw_-1) + nI_ ); // Overestimate, ToDo: compute the exact value!
	Gh_orig_nz_ = 2*mI_;
	//
	xinit_orig_.resize(nD_ + nI_);
	xl_orig_.resize(xinit_orig_.dim());
	xu_orig_.resize(xinit_orig_.dim());
	hl_orig_.resize(mI_);
	hu_orig_.resize(mI_);
	co_orig_.resize(nD_ + mU_);
	//
	xinit_orig_ = xo;
	xl_orig_    = xl;
	xu_orig_    = xu;
	if( mI_ ) {
		hl_orig_    = hl;
		hu_orig_    = hu;
	}
	co_orig_    = 0.0;
}

// Overridden public members from NLP

void ExampleNLPBanded::initialize(bool test_setup)
{
	if(is_initialized_) {
		NLPSerialPreprocessExplJac::initialize(test_setup);
		return;
	}
	// Nothing to initialize?
	NLPSerialPreprocessExplJac::initialize(test_setup);
	is_initialized_ = true;
}

bool ExampleNLPBanded::is_initialized() const
{
	return is_initialized_;
}

value_type ExampleNLPBanded::max_var_bounds_viol() const
{
	return +1e+20; // Functions defined everywhere!
}

void ExampleNLPBanded::set_multi_calc(bool multi_calc) const
{
	multi_calc_ = multi_calc;
}

bool ExampleNLPBanded::multi_calc() const
{
	return multi_calc_;
}

// Overridden from NLPVarReductPerm

bool ExampleNLPBanded::nlp_selects_basis() const
{
	return nlp_selects_basis_;
}

// Overridden protected methods from NLPSerialPreprocess */

bool ExampleNLPBanded::imp_nlp_has_changed() const
{
	return !is_initialized_;
}

size_type ExampleNLPBanded::imp_n_orig() const
{
	return nD_ + nI_;
}

size_type ExampleNLPBanded::imp_m_orig() const
{
	return nD_ + mU_;
}

size_type ExampleNLPBanded::imp_mI_orig() const
{
	return mI_;
}

const VectorSlice ExampleNLPBanded::imp_xinit_orig() const
{
	return xinit_orig_();
}

bool ExampleNLPBanded::imp_has_var_bounds() const
{
	return false; // ToDo: Add bounds!
}

const VectorSlice ExampleNLPBanded::imp_xl_orig() const
{
	return xl_orig_();
}

const VectorSlice ExampleNLPBanded::imp_xu_orig() const
{
	return xu_orig_();
}

const VectorSlice ExampleNLPBanded::imp_hl_orig() const
{
	return hl_orig_();
}

const VectorSlice ExampleNLPBanded::imp_hu_orig() const
{
	return hu_orig_();
}

void ExampleNLPBanded::imp_calc_f_orig(
	const VectorSlice            &x_full
	,bool                        newx
	,const ZeroOrderInfoSerial   &zero_order_info
	) const
{
	inform_new_point(newx);
	const VectorSlice x_orig = x_full(1,imp_n_orig());
	*zero_order_info.f = ( 1.0 / 2.0 ) * LinAlgPack::dot( x_orig, x_orig );
}

void ExampleNLPBanded::imp_calc_c_orig(
	const VectorSlice            &x_full
	,bool                        newx
	,const ZeroOrderInfoSerial   &zero_order_info
	) const
{
	inform_new_point(newx);
	if(c_orig_updated_)
		return; // c(x) is already computed in *zero_order_info.c
	Vector
		&c  = *zero_order_info.c;
	const size_type
		num_I_per_D = nD_ / nI_,  // Integer division (rounds down)
		I_remainder = nD_ % nI_;
	size_type j = 0;
	const value_type
		ds_alpha = nD_ > 1 ? diag_scal_ * (diag_vary_ - 1.0)/(nD_ - 1.0) : 0.0,
		ds_beta = diag_scal_;
	for( size_type q_i = 1; q_i <= nI_; ++q_i ) {
		const size_type num_I_per_D_local = num_I_per_D + ( q_i <= I_remainder ? 1 : 0 );
		for( size_type q_k = 0; q_k < num_I_per_D_local; ++q_k ) {
			++j;
			const size_type
				klu = ( j - bw_     >= 0   ? bw_-1 : j-1   ),
				kuu = ( j + bw_ - 1 <= nD_ ? bw_-1 : nD_-j );
			const value_type
				ds_j = ds_alpha * (j-1) + ds_beta;
			value_type
				&c_j = (c(j) = ds_j * x_full(j));
			{for( size_type k = 1; k <= klu; ++k ) {
				c_j -= (3.0 / k) * x_full(j-k);
			}}
			{for( size_type k = 1; k <= kuu; ++k ) {
				c_j -= (fu_ / k) * x_full(j+k);
			}}
			const value_type
				term = x_full(nD_ + q_i) + 1;
			c_j *= (term * term);
			c_j += co_orig_(j);
		}
	}
	c_orig_updated_ = true;
}

void ExampleNLPBanded::imp_calc_h_orig(
	const VectorSlice            &x_full
	,bool                        newx
	,const ZeroOrderInfoSerial   &zero_order_info
	) const
{
	inform_new_point(newx);
	Vector
		&h  = *zero_order_info.h;
	const size_type
		num_I_per_D = nD_ / nI_,  // Integer division (rounds down)
		I_remainder = nD_ % nI_;
	size_type jI = 0;
	for( size_type q_i = 1; q_i <= nI_; ++q_i ) {
		const value_type  x_q = x_full(nD_ + q_i);
		const size_type   num_I_per_D_local = num_I_per_D + ( q_i <= I_remainder ? 1 : 0 );
		for( size_type q_k = 0; q_k < num_I_per_D_local; ++q_k ) {
			++jI;
			if( jI > mI_ ) goto EXIT_LOOP;
			h(jI) = x_full(jI) - x_q;
		}
	}
EXIT_LOOP:
	jI; // Must have a statement here
}

void ExampleNLPBanded::imp_calc_Gf_orig(
	const VectorSlice            &x_full
	,bool                        newx
	,const ObjGradInfoSerial     &obj_grad_info
	) const
{
	inform_new_point(newx);
	const Range1D var_orig(1,imp_n_orig());
	(*obj_grad_info.Gf)(var_orig) = x_full(var_orig);
}

bool ExampleNLPBanded::imp_get_next_basis(
	IVector      *var_perm
	,IVector     *equ_perm
	,size_type   *rank
	)
{
	if(basis_selection_was_given_)
		return false; // Already gave this basis selection.
	// Select the first nD variables as basis variables which gives
	// a nice banded matrix for the basis matrix C.
	// Also, if the general inequality constraints are begin
	// converted to equalities with slacks, make the slack variables
	// basic variables also (put them first).
#ifdef _DEBUG
	assert(var_perm);
	assert(equ_perm);
	assert(rank);
#endif
	const bool         cinequtoequ = this->convert_inequ_to_equ();
	const size_type    n_orig = nD_ + nI_;
	const size_type    n_full = n_orig + ( cinequtoequ ? mI_ : 0 );
	const size_type    m_full = nD_    + ( cinequtoequ ? mI_ : 0 );
	var_perm->resize(n_full);
	equ_perm->resize(m_full);
	if( cinequtoequ && mI_ ) {
		index_type k;
		// basic variables
		for( k = 1; k <= mI_; ++k )         // Put slacks first
			(*var_perm)(k) = n_orig + k;
		for( k = 1; k <= nD_; ++k )         // Followed by xD
			(*var_perm)(mI_ + k) = k;
		// non-basic variables
		for( k = 1; k <= nI_; ++k )         // Followed by nI
			(*var_perm)(mI_+nD_ + k) = nD_ + k;
	}
	else {
		LinAlgPack::identity_perm( var_perm );
	}
	LinAlgPack::identity_perm( equ_perm );
	*rank = m_full;
	basis_selection_was_given_ = true;
	return true;
}

void ExampleNLPBanded::imp_report_orig_final_solution(
	const VectorSlice      &x_orig
	,const VectorSlice     *lambda_orig
	,const VectorSlice     *lambdaI_orig
	,const VectorSlice     *nu_orig
	,bool                  is_optimal
	) const
{
	// ToDo: Do something with the final soltuion?
}

// Overridden protected methods from NLPSerialPreprocessExplJac */

size_type ExampleNLPBanded::imp_Gc_nz_orig() const
{
	return Gc_orig_nz_;
}

size_type ExampleNLPBanded::imp_Gh_nz_orig() const
{
	return Gh_orig_nz_;
}

void ExampleNLPBanded::imp_calc_Gc_orig(
	const VectorSlice& x_full, bool newx
	, const FirstOrderExplInfo& first_order_expl_info
	) const
{
	inform_new_point(newx);
	// Compute c(x) if not already
	this->imp_calc_c_orig( x_full, newx, zero_order_orig_info() );
	Vector
		&c = *first_order_expl_info.c; // This must not be NULL!
	// Get references/pointers to data for Gc to be computed/updated.
	index_type
		&Gc_nz = *first_order_expl_info.Gc_nz;
	value_type
		*Gc_val = &(*first_order_expl_info.Gc_val)[0];
	index_type
		*Gc_ivect = ( first_order_expl_info.Gc_ivect
					  ? &(*first_order_expl_info.Gc_ivect)[0] : NULL ),
		*Gc_jvect = ( first_order_expl_info.Gc_jvect
					  ? &(*first_order_expl_info.Gc_jvect)[0] : NULL );
	assert( (Gc_ivect != NULL) == (Gc_jvect != NULL) );
	// Set nonzeros for Gc (in sorted compressed column format w.r.t., i.e. grouped by constraints)
	const size_type
		num_I_per_D = nD_ / nI_,  // Integer division (rounds down)
		I_remainder = nD_ % nI_;
	Gc_nz = 0;
	size_type j = 0;
	const value_type
		ds_alpha = nD_ > 1 ? diag_scal_ * (diag_vary_ - 1.0)/(nD_ - 1.0) : 0.0,
		ds_beta = diag_scal_;
	for( size_type q_i = 1; q_i <= nI_; ++q_i ) {
		const value_type
			x_q = x_full(nD_ + q_i),
			x_q_term = (x_q + 1) * (x_q + 1);
		const size_type num_I_per_D_local = num_I_per_D + ( q_i <= I_remainder ? 1 : 0 );
		for( size_type q_k = 0; q_k < num_I_per_D_local; ++q_k ) {
			++j;
			const size_type
				klu = ( j - bw_     >= 0   ? bw_-1 : j-1   ),
				kuu = ( j + bw_ - 1 <= nD_ ? bw_-1 : nD_-j );
			const value_type
				ds_j = ds_alpha * (j-1) + ds_beta;
			//
			{for( index_type k = klu; k >= 1; --k ) {
				++Gc_nz;
				*Gc_val++ = -3.0 / k * x_q_term;
				if(Gc_ivect) {
					*Gc_ivect++ = j - k;
					*Gc_jvect++ = j;
				}
			}}
			//
			++Gc_nz;
			*Gc_val++ = ds_j * x_q_term;
			if(Gc_ivect) {
				*Gc_ivect++ = j;
				*Gc_jvect++ = j;
			}
			//
			{for( index_type k = 1; k <= kuu; ++k ) {
				++Gc_nz;
				*Gc_val++ = -fu_ / k * x_q_term;
				if(Gc_ivect) {
					*Gc_ivect++ = j + k;
					*Gc_jvect++ = j;
				}
			}}
			//
			++Gc_nz;
			*Gc_val++ = 2.0 * (c(j) - co_orig_(j)) / (x_q + 1);
			if(Gc_ivect) {
				*Gc_ivect++ = nD_ + q_i;
				*Gc_jvect++ = j;
			}
		}
	}
}

void ExampleNLPBanded::imp_calc_Gh_orig(
	const VectorSlice& x_full, bool newx
	, const FirstOrderExplInfo& first_order_expl_info
	) const
{
	inform_new_point(newx);
	// Get references/pointers to data for Gh to be computed/updated.
	index_type
		&Gh_nz = *first_order_expl_info.Gh_nz;
	value_type
		*Gh_val = &(*first_order_expl_info.Gh_val)[0];
	index_type
		*Gh_ivect = ( first_order_expl_info.Gh_ivect
					  ? &(*first_order_expl_info.Gh_ivect)[0] : NULL ),
		*Gh_jvect = ( first_order_expl_info.Gh_jvect
					  ? &(*first_order_expl_info.Gh_jvect)[0] : NULL );
	assert( (Gh_ivect != NULL) == (Gh_jvect != NULL) );
	// Set nonzeros for Gh (in sorted compressed column format w.r.t., i.e. grouped by constraints)
	const size_type
		num_I_per_D = nD_ / nI_,  // Integer division (rounds down)
		I_remainder = nD_ % nI_;
	Gh_nz = 0;
	size_type jI = 0;
	for( size_type q_i = 1; q_i <= nI_; ++q_i ) {
		const size_type   nD_q_i = nD_ + q_i;
		const value_type  x_q = x_full(nD_q_i);
		const size_type   num_I_per_D_local = num_I_per_D + ( q_i <= I_remainder ? 1 : 0 );
		for( size_type q_k = 0; q_k < num_I_per_D_local; ++q_k ) {
			++jI;
			if( jI > mI_ ) goto EXIT_LOOP;
			// w.r.t. x(jI)
			++Gh_nz;
			*Gh_val++ = 1.0;
			if(Gh_ivect) {
				*Gh_ivect++ = jI;
				*Gh_jvect++ = jI;
			}
			// w.r.t. x(nD+q(jI))
			++Gh_nz;
			*Gh_val++ = -1.0;
			if(Gh_ivect) {
				*Gh_ivect++ = nD_q_i;
				*Gh_jvect++ = jI;
			}
		}
	}
EXIT_LOOP:
	jI; // Must have a statement here
}

// private

void ExampleNLPBanded::assert_is_initialized() const
{
	assert(0); //  ToDo: Implemenet!
}

void ExampleNLPBanded::inform_new_point(bool newx) const
{
	if(newx) {
		c_orig_updated_ = false;
	}
}

}	// end namespace NLPInterfacePack
