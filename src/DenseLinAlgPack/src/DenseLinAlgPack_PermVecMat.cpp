// ///////////////////////////////////////////////////////////////////////////////////////
// PermVecMat.cpp
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

#include "LinAlgPack/include/PermVecMat.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/IVector.h"

#ifdef _DEBUG   // Debug only!
bool LinAlgPack::PermVecMat_print = false;
#include <iostream>
#include "LinAlgPack/include/PermOut.h"
#include "LinAlgPack/include/VectorOut.h"
#endif

// Local assert function
namespace {
inline void i_assert_perm_size(size_t size1, size_t size2)
{
#ifdef LINALGPACK_CHECK_RHS_SIZES
	if(size1 != size2)
		throw std::length_error("The size of the permutation vector is not correct");
#endif
}

} // end namespace

void LinAlgPack::identity_perm(IVector* perm) {
	if(!perm->size())
		throw std::length_error("LinAlgPack::identity_perm(): perm must be sized");
	IVector::iterator	itr_perm		= perm->begin();
	for(size_type i = 1; i <= perm->size(); ++i)
		*itr_perm++ = i;
}

void LinAlgPack::inv_perm(const IVector& perm, IVector* inv_perm) {
	inv_perm->resize(perm.size());
	for(size_type i = 1; i <= perm.size(); ++i)
		(*inv_perm)(perm(i)) = i;
}

void LinAlgPack::perm_ele(const IVector& perm, VectorSlice* vs)
{
	i_assert_perm_size(vs->dim(),perm.size());
	Vector tmp_v(vs->dim());
	Vector::iterator		v_itr		= tmp_v.begin(),
							v_itr_end	= tmp_v.end();
	IVector::const_iterator perm_itr = perm.begin();
	// Copy the elements in the permed order into the temp vector
	for(; v_itr != v_itr_end; ++v_itr, ++perm_itr)
		*v_itr = (*vs)(*perm_itr);
		
	// Copy the permed vector back
	(*vs) = tmp_v;
}

void LinAlgPack::perm_ele(const VectorSlice& x, const IVector& perm, VectorSlice* y)
{
	i_assert_perm_size(x.dim(),perm.size());
	i_assert_perm_size(y->dim(),perm.size());

	IVector::const_iterator
		perm_itr	= perm.begin();
	VectorSlice::iterator
		y_itr		= y->begin(),
		y_end		= y->end();
	while(y_itr != y_end)
		*y_itr++ = x(*perm_itr++);
}

void LinAlgPack::inv_perm_ele(const VectorSlice& y, const IVector& perm, VectorSlice* x)
{
	i_assert_perm_size(y.dim(),perm.size());
	i_assert_perm_size(x->dim(),perm.size());
#ifdef _DEBUG
	if( PermVecMat_print ) {
		std::cerr
			<< "enter inv_perm_ele(y,perm,x):\n"
			<< "before:\n"
			<< "y =\n" << y
			<< "perm =\n" << perm
			<< "x =\n" << *x;
	}
#endif	
	VectorSlice::const_iterator
		y_itr		= y.begin(),
		y_end		= y.end();
	IVector::const_iterator
		perm_itr	= perm.begin();
	while(y_itr != y_end)
		(*x)(*perm_itr++) = *y_itr++;
#ifdef _DEBUG
	if( PermVecMat_print ) {
		std::cerr
			<< "inv_perm_ele(y,perm,x):\n"
			<< "after:\n"
			<< "x =\n" << *x
			<< "exit inv_perm_ele(...) ...\n";
	}
#endif	
}

void LinAlgPack::perm_rows(const IVector& row_perm, GenMatrixSlice* gms)
{
	i_assert_perm_size(gms->rows(),row_perm.size());
	GenMatrix tmp_gm(gms->rows(),gms->cols());
	GenMatrixSlice::size_type rows = gms->rows(), i;
	// Copy the rows in the correct order into the temp matrix.
	for(i = 1; i <= rows; ++i)
		tmp_gm.row(i) = gms->row(row_perm(i));
	// Copy the permed matrix back
	(*gms) = tmp_gm;
}

void LinAlgPack::perm_cols(const IVector& col_perm, GenMatrixSlice* gms)
{
	i_assert_perm_size(gms->cols(),col_perm.size());
	GenMatrix tmp_gm(gms->rows(),gms->cols());
	GenMatrixSlice::size_type cols = gms->cols(), i;
	// Copy the columns in the correct order into the temp matrix.
	for(i = 1; i <= cols; ++i)
		tmp_gm.col(i) = gms->col(col_perm(i));
	// Copy the permed matrix back
	(*gms) = tmp_gm;
}

void LinAlgPack::perm_rows_cols(const IVector& row_perm, const IVector& col_perm
	, GenMatrixSlice* gms)
{
	i_assert_perm_size(gms->rows(),row_perm.size());
	i_assert_perm_size(gms->cols(),col_perm.size());
	GenMatrix tmp_gm(gms->rows(),gms->cols());
	GenMatrixSlice::size_type rows = gms->rows(), cols = gms->cols(), i;
	// Copy the rows in the correct order into the temp matrix.
	for(i = 1; i <= rows; ++i)
		tmp_gm.row(i) = gms->row(row_perm(i));
	// Copy the columns in the correct order back into matrix.
	for(i = 1; i <= cols; ++i)
		gms->col(i) = tmp_gm.col(col_perm(i));
}
