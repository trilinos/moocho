// ///////////////////////////////////////////////////////////////////////
// NLP.cpp
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

#include <limits>

#include "NLPInterfacePack/src/abstract/interfaces/NLP.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorMutable.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorSpace.hpp"
#include "Teuchos_TestForException.hpp"

namespace {
const char name_f[] = "f";
const char name_c[] = "c";
const char name_c_breve[] = "c_breve";
const char name_h_breve[] = "h_breve";
NLPInterfacePack::NLP::options_ptr_t  null_options = Teuchos::null;
} // end namespace

namespace NLPInterfacePack {

// static

value_type NLP::infinite_bound()
{
//	return std::numeric_limits<value_type>::max();
	return 1e+50;
}

// constructors

NLP::NLP()
{}

// destructor

NLP::~NLP()
{}

void NLP::set_options( const options_ptr_t& options )
{}

const NLP::options_ptr_t&
NLP::get_options() const
{
	return null_options;
}

void NLP::initialize(bool test_setup)
{
	num_f_evals_ = num_c_evals_ = 0;
}

// dimensionality

size_type
NLP::n() const
{
	return this->space_x()->dim();
}

size_type 
NLP::m() const
{
	VectorSpace::space_ptr_t spc = this->space_c();
	return spc.get() ? spc->dim() : 0;
}

// initial guess

void NLP::get_init_lagrange_mult(
	VectorMutable*   lambda
	,VectorMutable*  nu
	) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( lambda  && this->m()  == 0,            std::logic_error, "" );
	TEST_FOR_EXCEPTION( nu      && this->num_bounded_x() == 0, std::logic_error, "" );
#endif
	if(lambda) {
#ifdef _DEBUG
		TEST_FOR_EXCEPTION( !this->space_c()->is_compatible(lambda->space()), VectorSpace::IncompatibleVectorSpaces, "" );
#endif
		*lambda = 0.0;
	}
	if(nu) {
#ifdef _DEBUG
		TEST_FOR_EXCEPTION( !this->space_x()->is_compatible(nu->space()), VectorSpace::IncompatibleVectorSpaces, "" );
#endif
		*nu = 0.0;
	}
}

// <<std comp>> members for f

void NLP::set_f(value_type* f)
{
	first_order_info_.f = f;
}

value_type* NLP::get_f()
{
	return StandardCompositionRelationshipsPack::get_role_name(first_order_info_.f, false, name_f);
}

value_type& NLP::f()
{
	return StandardCompositionRelationshipsPack::role_name(first_order_info_.f, false, name_f);
}

const value_type& NLP::f() const
{
	return StandardCompositionRelationshipsPack::role_name(first_order_info_.f, false, name_f);
}

// <<std comp>> members for c

void NLP::set_c(VectorMutable* c)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
	TEST_FOR_EXCEPTION( c && !this->space_c()->is_compatible(c->space()), VectorSpace::IncompatibleVectorSpaces, "" );
#endif
	first_order_info_.c = c;
}

VectorMutable* NLP::get_c()
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::get_role_name(first_order_info_.c, false, name_c);
}

VectorMutable& NLP::c()
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(first_order_info_.c, false, name_c);
}

const Vector& NLP::c() const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(first_order_info_.c, false, name_c);
}

void NLP::unset_quantities()
{
	first_order_info_.f = NULL;
	first_order_info_.c = NULL;
	first_order_info_.h = NULL;
	first_order_info_breve_.f = NULL;
	first_order_info_breve_.c = NULL;
	first_order_info_breve_.h = NULL;
}

// calculations

void NLP::calc_f(const Vector& x, bool newx) const
{
	StandardCompositionRelationshipsPack::assert_role_name_set(first_order_info_.f, "NLP::calc_f()", name_f);
	imp_calc_f(x,newx,zero_order_info());
	num_f_evals_++;
}

void NLP::calc_c(const Vector& x, bool newx) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	StandardCompositionRelationshipsPack::assert_role_name_set(first_order_info_.c, "NLP::calc_c()", name_c);
	imp_calc_c(x,newx,zero_order_info());
	num_c_evals_++;
}

void NLP::report_final_solution(
	const Vector&    x
	,const Vector*   lambda
	,const Vector*   nu
	,bool            optimal
	) const
{
	// The default behavior is just to ignore the solution!
}

size_type NLP::num_f_evals() const
{
	return num_f_evals_;
}

size_type NLP::num_c_evals() const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0, std::logic_error, "" );
#endif
	return num_c_evals_;
}

// General inequalities and slack variables

size_type NLP::ns() const
{
	vec_space_ptr_t space_h_breve = this->space_h_breve();
	return space_h_breve.get() ? space_h_breve->dim() : 0;
}

NLP::vec_space_ptr_t NLP::space_c_breve() const
{
	return this->space_c();
}

NLP::vec_space_ptr_t NLP::space_h_breve() const
{
	return Teuchos::null;
}

const Vector& NLP::hl_breve() const
{
	TEST_FOR_EXCEPTION(
		true, std::logic_error
		,"NLP::hl_breve(): Error, this method must be overridden if space_h_breve is defined" );
}

const Vector& NLP::hu_breve() const
{
	TEST_FOR_EXCEPTION(
		true, std::logic_error
		,"NLP::hl_breve(): Error, this method must be overridden if space_h_breve is defined" );
}

void NLP::set_c_breve(VectorMutable* c_breve)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() - this->ns() == 0, std::logic_error, "" );
	TEST_FOR_EXCEPTION( c_breve && !this->space_c_breve()->is_compatible(c_breve->space()), VectorSpace::IncompatibleVectorSpaces, "" );
#endif
	first_order_info_breve_.c = c_breve;
}

VectorMutable* NLP::get_c_breve()
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() - this->ns() == 0, std::logic_error, "" );
#endif
	return first_order_info_breve_.c;
}

VectorMutable& NLP::c_breve()
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() - this->ns() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(first_order_info_breve_.c, false, name_c_breve);
}

const Vector& NLP::c_breve() const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() - this->ns() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(first_order_info_breve_.c, false, name_c_breve);
}

void NLP::set_h_breve(VectorMutable* h_breve)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() - this->ns() == 0, std::logic_error, "" );
	TEST_FOR_EXCEPTION( h_breve && !this->space_h_breve()->is_compatible(h_breve->space()), VectorSpace::IncompatibleVectorSpaces, "" );
#endif
	first_order_info_breve_.c = h_breve;
}

VectorMutable* NLP::get_h_breve()
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() - this->ns() == 0, std::logic_error, "" );
#endif
	return first_order_info_breve_.h;
}

VectorMutable& NLP::h_breve()
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() - this->ns() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(first_order_info_breve_.c, false, name_h_breve);
}

const Vector& NLP::h_breve() const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() - this->ns() == 0, std::logic_error, "" );
#endif
	return StandardCompositionRelationshipsPack::role_name(first_order_info_breve_.c, false, name_h_breve);
}

const Permutation& NLP::P_var() const
{
	assert(0);
//	if(!P_var_.get()) P_var_ = Teuchos::rcp(new PermutationSerial(this->space_x());
	return *P_var_;
}

const Permutation& NLP::P_equ() const
{
	assert(0);
//	if(!P_equ_.get()) P_equ = Teuchos::rcp(new PermutationSerial(this->space_c());
	return *P_equ_;
}

void NLP::calc_c_breve(const Vector& x, bool newx) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->m() == 0 || this->ns() > 0, std::logic_error, "" );
#endif
	StandardCompositionRelationshipsPack::assert_role_name_set(first_order_info_breve_.c, "NLP::calc_c_breve()", name_c_breve);
	imp_calc_c_breve(x,newx,zero_order_info_breve());
	num_c_evals_++;
}

void NLP::calc_h_breve(const Vector& x, bool newx) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( this->ns() == 0, std::logic_error, "" );
#endif
	StandardCompositionRelationshipsPack::assert_role_name_set(first_order_info_breve_.h, "NLP::calc_h_breve()", name_h_breve);
	imp_calc_c_breve(x,newx,zero_order_info_breve());
	num_c_evals_++;
}

// protected

void NLP::imp_calc_c_breve(
	const Vector           &x
	,bool                  newx
	,const ZeroOrderInfo   &zero_order_info_breve
	) const
{
	imp_calc_c(x,newx,zero_order_info_breve);
}

void NLP::imp_calc_h_breve(
	const Vector           &x
	,bool                  newx
	,const ZeroOrderInfo   &zero_order_info_breve
	) const
{
	TEST_FOR_EXCEPTION(
		true, std::logic_error
		,"NLP::hl_breve(): Error, this method must be overridden if space_h_breve is defined" );
}

} // namespace NLPInterfacePack
