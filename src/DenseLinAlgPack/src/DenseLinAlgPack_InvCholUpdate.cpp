// /////////////////////////////////////////////////////////////////
// InvCholUpdate.cpp
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

#include "InvCholUpdate.h"
#include "GenMatrixClass.h"
#include "VectorOp.h"

namespace {
// sign function
inline int sign(double v) {
	if(v > 0.0)
		return 1;
	else if(v == 0.0)
		return 0;
	else
		return -1;
}

}

void LinAlgPack::update_chol_factor(GenMatrixSlice* pM, VectorSlice* pu
	, const VectorSlice& v)
{
	GenMatrixSlice	&M = *pM;
	VectorSlice		&u = *pu;

	assert_gms_square(M);
	if(M.rows() != u.dim() || u.dim() != v.dim())
		throw std::length_error("update_chol_factor(): Sizes of matrix and vectors do not match");

	// Algorithm A3.4.1a in Dennis and Schabel
	
	// 0. Set lower diagonal to zero
	M.diag(-1) = 0.0;

	// 1,2 Find the largest k such that u(k) != 0.0
	size_type n = M.rows(), k = n;
	{
		VectorSlice::reverse_iterator r_itr_u = u.rbegin();
		while(*r_itr_u == 0.0 && k > 1) --k;
	}
	
	// 3.
	{
		VectorSlice::reverse_iterator
			r_itr_u_i		= u(1,k-1).rbegin(),	// iterator for u(i), i = k-1,...,1
			r_itr_u_ip1		= u(2,k).rbegin();		// iterator for u(i), i = k,...,2
		for(size_type i = k-1; i > 0 ; ++r_itr_u_i, ++r_itr_u_ip1, --i) {
			value_type u_i = *r_itr_u_i, u_ip1 = *r_itr_u_ip1;
			jacobi_rotate(&M, i, u_i, - u_ip1);	// 3.1
			*r_itr_u_i = (u_i == 0.0) ? ::fabs(u_ip1) : ::sqrt(u_i * u_i + u_ip1 * u_ip1); // 3.2
		}
	}

	// 4. M.row(1) += u(1) * v 
	LinAlgPack::Vp_StV(&M.row(1), u(1), v);

	// 5.
	{
		VectorSlice::const_iterator
			itr_M_i_i	= M.diag().begin(),		// iterator for M(i,i), for i = 1,...,k-1
			itr_M_ip1_i	= M.diag(-1).begin();	// iterator for M(i+1,i), for i = 1,...,k-1
		for(size_type i = 1; i < k; ++i)
			jacobi_rotate(&M, i, *itr_M_i_i++, - *itr_M_ip1_i++);
	}

/*  This code does not work

	size_type n = M.rows();

	// 0.
	{
		for(size_type i = 2; i <= n; ++i)
			M(i,i-1) = 0;
	}

	// 1.
	size_type k = n;
	
	// 2.
	while(u(k) == 0.0 && k > 1) --k;

	// 3.
	{
		for(size_type i = k-1; i >= 1; --i) {
			jacobi_rotate(M, i, u(i), - u(i+1));
			if(u(i) == 0)
				u(i) = ::fabs(u(i+1));
			else
				u(i) = ::sqrt(u(i) * u(i) + u(i+1) * u(i+1));
		}
	}

	// 4.
	{
		for(size_type j = 1; j <= n; ++j)
			M(1,j) = M(1,j) + dot(u(1),v(j));
	}

	// 5.
	{
		for(size_type i = 1; i <= k-1; ++i)
			jacobi_rotate(M, i, M(i,i), - M(i+1,i));
	}



*/


}

void LinAlgPack::jacobi_rotate(GenMatrixSlice* pM, size_type row_i, value_type alpha
	, value_type beta)
{
	GenMatrixSlice	&M = *pM;

	assert_gms_square(M);
		
	// Algorithm A3.4.1a in Dennis and Schabel
	
	// 1.
	value_type c, s;
	if(alpha == 0.0) {
		c = 0.0;
		s = sign(beta);
	}
	else {
		value_type den = ::sqrt(alpha * alpha + beta * beta);
		c = alpha / den;
		s = beta / den;
	}

	// 2.
	size_type i = row_i, n = M.rows();

	// Use iterators instead of element access
	VectorSlice::iterator
		itr_M_i		= M.row(i)(i,n).begin(),	// iterator for M(i,j), for j = i,...,n
		itr_M_i_end	= M.row(i)(i,n).end(),
		itr_M_ip1	= M.row(i+1)(i,n).begin();	// iterator for M(i+1,j), for j = i,...,n

	while(itr_M_i != itr_M_i_end) {
		value_type y = *itr_M_i, w = *itr_M_ip1;
		*itr_M_i++		= c * y - s * w;
		*itr_M_ip1++	= s * y + c * w;
	}

/*  This code does not work

	size_type n = M.rows(), i = row_i;

	// 1.
	value_type c, s;
	if(alpha == 0) {
		c = 0.0;
		s = ::fabs(beta);
	}
	else {
		value_type den = ::sqrt(alpha*alpha + beta*beta);
		c = alpha / den;
		s = beta / den;
	}

	// 2.
	{
		for(size_type j = i; j <= n; ++j) {
			size_type y = M(i,j), w = M(i+1,j);
			M(i,j) = c*y - s*w;
			M(i+1,j) = s*y + c*w;
		}
	}



*/

}
