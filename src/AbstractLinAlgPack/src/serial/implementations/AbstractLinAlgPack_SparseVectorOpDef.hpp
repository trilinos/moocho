// //////////////////////////////////////////////////////////////////////
// SparseVectorOpDef.h
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
//
// The declarations for these functions is in the file SparseVectorOpDecl.h
// but because of a bug with the MS VC++ 5.0 compiler you can not use
// namespace qualification with definitions of previously declared
// nonmember template funcitons.  By not including the declarations
// and by including this file for automatic instantiation, then
// if the function prototypes are not the same then a compile
// time error is more likely to occur.  Otherwise you could have
// to settle for a compile-time warning that the funciton has
// not been defined or a link-time error that the definition
// could not be found which will be the case when explicit
// instantiation is used.

// ToDo: 6/9/98 Finish upgrade

#ifndef SPARSE_VECTOR_OP_DEF_H
#define SPARSE_VECTOR_OP_DEF_H

#include "SparseLinAlgPackTypes.h"
#include "AbstractLinAlgPack/src/SparseVectorClass.h"
#include "LinAlgPack/src/VectorOp.h"
#include "LinAlgPack/src/GenMatrixAsTriSym.h"	// also included in SparseVectorOpDef.h
#include "LinAlgPack/src/GenMatrixClass.h"
#include "LinAlgPack/src/LinAlgPackAssertOp.h"

namespace {
template< class T >
inline
T my_my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
template< class T >
inline
T my_my_min( const T& v1, const T& v2 ) { return v1 < v2 ? v1 : v2; }
} // end namespace

namespace SparseLinAlgPack {

using LinAlgPack::VopV_assert_sizes;
using LinAlgPack::Vp_V_assert_sizes;
using LinAlgPack::Vp_MtV_assert_sizes;

using LinAlgPack::row;
using LinAlgPack::col;

namespace SparseVectorUtilityPack {
template<class T_SpVec>
value_type imp_dot2_V_V_SV(const VectorSlice& vs1, const VectorSlice& vs2, const T_SpVec& sv);
}

// result = dot(vs_rhs1,sv_rhs2)
template<class T_SpVec>
value_type dot_V_SV(const VectorSlice& vs_rhs1, const T_SpVec& sv_rhs2) {
	VopV_assert_sizes(vs_rhs1.dim(),sv_rhs2.dim());
	value_type result = 0.0;
	typename T_SpVec::difference_type offset = sv_rhs2.offset();
	for(typename T_SpVec::const_iterator iter = sv_rhs2.begin(); iter != sv_rhs2.end(); ++iter)
		result += vs_rhs1(iter->index()+offset) * iter->value();
	return result;
}

// result = dot(sv_rhs1,vs_rhs2).  Just call the above in reverse order
template<class T_SpVec>
value_type dot_SV_V(const T_SpVec& sv_rhs1, const VectorSlice& vs_rhs2) {
	return dot_V_SV(vs_rhs2,sv_rhs1);
}

// result = ||sv_rhs||1
template<class T_SpVec>
value_type norm_1_SV(const T_SpVec& sv_rhs) {
	typename T_SpVec::element_type::value_type result = 0.0;
	for(typename T_SpVec::const_iterator iter = sv_rhs.begin(); iter != sv_rhs.end(); ++iter)
		result += ::fabs(iter->value());
	return result;
}

// result = ||sv_rhs||2
template<class T_SpVec>
value_type norm_2_SV(const T_SpVec& sv_rhs) {
	typename T_SpVec::element_type::value_type result = 0.0;
	for(typename T_SpVec::const_iterator iter = sv_rhs.begin(); iter != sv_rhs.end(); ++iter)
		result += (iter->value()) * (iter->value());
	return result;
}

// result = ||sv_rhs||inf
template<class T_SpVec>
value_type norm_inf_SV(const T_SpVec& sv_rhs) {
	typename T_SpVec::element_type::value_type result = 0.0;
	for(typename T_SpVec::const_iterator iter = sv_rhs.begin(); iter != sv_rhs.end(); ++iter)
		result = my_my_max(result,::fabs(iter->value()));
	return result;
}

// result = max(sv_rhs)
template<class T_SpVec>
value_type max_SV(const T_SpVec& sv_rhs) {
	typename T_SpVec::element_type::value_type result = 0.0;
	for(typename T_SpVec::const_iterator iter = sv_rhs.begin(); iter != sv_rhs.end(); ++iter)
		result = my_my_max(iter->value(),result);
	return result;
}

// result = min(sv_rhs)
template<class T_SpVec>
value_type min_SV(const T_SpVec& sv_rhs) {
	typename T_SpVec::element_type::value_type result = 0.0;
	for(typename T_SpVec::const_iterator iter = sv_rhs.begin(); iter != sv_rhs.end(); ++iter)
		result = my_my_min(result,iter->value());
	return result;
}

// vs_lhs += alpha * sv_rhs (BLAS xAXPY)
template<class T_SpVec>
void Vt_S( T_SpVec* sv_lhs, value_type alpha )
{
	if( alpha == 1.0 ) return;
	for(typename T_SpVec::iterator iter = sv_lhs->begin(); iter != sv_lhs->end(); ++iter)
		iter->value() *= alpha;
}

// vs_lhs += alpha * sv_rhs (BLAS xAXPY)
template<class T_SpVec>
void Vp_StSV(VectorSlice* vs_lhs, value_type alpha, const T_SpVec& sv_rhs)
{
	Vp_V_assert_sizes(vs_lhs->dim(),sv_rhs.dim());
	typename T_SpVec::difference_type offset = sv_rhs.offset();
	for(typename T_SpVec::const_iterator iter = sv_rhs.begin(); iter != sv_rhs.end(); ++iter)
		(*vs_lhs)(iter->index() + offset) += alpha * iter->value();
}

// vs_lhs += alpha * op(gms_rhs1) * sv_rhs2 (BLAS xGEMV) (time = O(sv_rhs2.nz() * vs_lhs.dim())
template<class T_SpVec>
void Vp_StMtSV(VectorSlice* pvs_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const T_SpVec& sv_rhs2)
{
#ifdef _WINDOWS
	using LinAlgPack::Vp_StV;	// MS VC++ 6.0 needs help with the name lookups
#endif
	VectorSlice& vs_lhs = *pvs_lhs;

	Vp_MtV_assert_sizes(vs_lhs.dim(),gms_rhs1.rows(),gms_rhs1.cols(),trans_rhs1
						, sv_rhs2.dim());
	
	// Perform the operation by iterating through the sparse vector and performing
	// all of the operations on it.
	//
	// For sparse element e we do the following:
	//
	// vs_lhs += alpha * e.value() * gms_rhs1.col(e.index());

	typename T_SpVec::difference_type offset = sv_rhs2.offset();

	for(typename T_SpVec::const_iterator sv_rhs2_itr = sv_rhs2.begin(); sv_rhs2_itr != sv_rhs2.end(); ++sv_rhs2_itr)
		LinAlgPack::Vp_StV( &vs_lhs, alpha * sv_rhs2_itr->value()
            , col( gms_rhs1, trans_rhs1, sv_rhs2_itr->index() + offset ) );
}

// vs_lhs += alpha * op(tri_rhs1) * sv_rhs2 (BLAS xTRMV)
template<class T_SpVec>
void Vp_StMtSV(VectorSlice* pvs_lhs, value_type alpha, const tri_gms& tri_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const T_SpVec& sv_rhs2)
{
	VectorSlice &vs_lhs = *pvs_lhs;

	Vp_MtV_assert_sizes(vs_lhs.dim(),tri_rhs1.rows(),tri_rhs1.cols(),trans_rhs1
						, sv_rhs2.dim());

	// Get the effective matrix
	BLAS_Cpp::Uplo effective_uplo;
	if( (tri_rhs1.uplo() == BLAS_Cpp::lower && trans_rhs1 == BLAS_Cpp::no_trans) ||
		(tri_rhs1.uplo() == BLAS_Cpp::upper && trans_rhs1 == BLAS_Cpp::trans)		   )
	{
		effective_uplo = BLAS_Cpp::lower;
	}
	else {	// must be effective upper
		effective_uplo = BLAS_Cpp::upper;
	}

	size_type n = tri_rhs1.gms().rows();	// should be same as cols()

	// Implement the operation by looping through the sparse vector only once
	// and performing the row operations.  This gives a time = O(n * sv_rhs2.nz())
	typename T_SpVec::difference_type offset = sv_rhs2.offset();
	for(typename T_SpVec::const_iterator sv_itr = sv_rhs2.begin(); sv_itr != sv_rhs2.end(); ++sv_itr)
	{
		size_type j = sv_itr->index() + offset;

		// For the nonzero element j = sv_itr->index() we perfom the following
		// operations.
		//
		// Lower:
		//		[\]		[\ 0 0 0]	[\]
		//		[#]	 +=	[\ # 0 0] * [#] jth element
		//		[#]		[\ # \ 0]	[\]
		//		[#]		[\ # \ \]	[\]
		//				  jth
		//				  col
		//
		// Upper:
		//		[#]		[\ # \ \]	[\]
		//		[#]	 +=	[0 # \ \] * [#] jth element
		//		[\]		[0 0 \ \]	[\]
		//		[\]		[0 0 0 \]	[\]
		//				  jth
		//				  col
		//
		// If we were told that is it is unit diagonal then we will adjust
		// accordingly.

		size_type j_adjusted = j;	// will be adjusted for unit diagonal
		
		switch(effective_uplo) {
			case BLAS_Cpp::lower: {
				if(tri_rhs1.diag() == BLAS_Cpp::unit)
				{
					// Make the adjustment for unit diaganal
					++j_adjusted;
					vs_lhs(j) += alpha * sv_itr->value();	// diagonal element is one
				}
				// vs_lhs(j,n) = vs_lhs(j,n) + alpha * sv_itr->value() * tri_rhs1.col(j)(j,n)
				if(j_adjusted <= n)
				{
					LinAlgPack::Vp_StV( &vs_lhs(j_adjusted,n), alpha * sv_itr->value()
						,col(tri_rhs1.gms(),trans_rhs1,j)(j_adjusted,n) );
				}
				break;
			}
			case BLAS_Cpp::upper: {
				if(tri_rhs1.diag() == BLAS_Cpp::unit)
				{
					// Make the adjustment for unit diaganal
					--j_adjusted;
					vs_lhs(j) += alpha * sv_itr->value();	// diagonal element is one
				}
				// vs_lhs(1,j) = vs_lhs(1,j) + alpha * sv_itr->value() * tri_rhs1.col(j)(1,j)
				if(j_adjusted > 0)
				{
					LinAlgPack::Vp_StV( &vs_lhs(1,j_adjusted), alpha * sv_itr->value()
						,col(tri_rhs1.gms(),trans_rhs1,j)(1,j_adjusted) );
				}
				break;
			}
		}
	}
}

// vs_lhs += alpha * op(sym_rhs1) * sv_rhs2 (BLAS xSYMV)
template<class T_SpVec>
void Vp_StMtSV(VectorSlice* pvs_lhs, value_type alpha, const sym_gms& sym_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const T_SpVec& sv_rhs2)
{
	VectorSlice& vs_lhs = *pvs_lhs;

	Vp_MtV_assert_sizes(vs_lhs.dim(),sym_rhs1.rows(),sym_rhs1.cols(),trans_rhs1
						, sv_rhs2.dim());

	size_type size = sv_rhs2.dim();
	switch(sym_rhs1.uplo()) {
		case BLAS_Cpp::lower: {
			VectorSlice::iterator vs_lhs_itr; size_type i;
			for(vs_lhs_itr = vs_lhs.begin(), i = 1; i <= size; ++i)
			{
				if(i < size) {
					*vs_lhs_itr++ +=
						alpha *
						SparseVectorUtilityPack::imp_dot2_V_V_SV(
										 sym_rhs1.gms().row(i)(1,i)
										,sym_rhs1.gms().col(i)(i+1,size)
										,sv_rhs2);
				}
				else
					*vs_lhs_itr++ += alpha *
						dot_V_SV(sym_rhs1.gms().row(i),sv_rhs2);
			}
			break;			
		}
		case BLAS_Cpp::upper: {
			VectorSlice::iterator vs_lhs_itr; size_type i;
			for(vs_lhs_itr = vs_lhs.begin(), i = 1; i <= size; ++i)
			{
				if(i > 1) {
					*vs_lhs_itr++ +=
						alpha *
						SparseVectorUtilityPack::imp_dot2_V_V_SV(
										 sym_rhs1.gms().col(i)(1,i-1)
										,sym_rhs1.gms().row(i)(i,size)
										,sv_rhs2);
				}
				else
				*vs_lhs_itr++ += alpha * dot_V_SV(sym_rhs1.gms().row(i),sv_rhs2);
			}
			break;			
		}
	}	
}

namespace SparseVectorUtilityPack {

// Implementation for the product of a concatonated dense vector with a
// sparse vector.  Used for symetric matrix mulitplication.
// In Matlab notation: result = [vs1' , vs2' ] * sv
// where split = vs1.dim(), vs2.dim() == sv.dim() - split
//
// time = O(sv.nz()), space = O(1)
//
template<class T_SpVec>
value_type imp_dot2_V_V_SV(const VectorSlice& vs1, const VectorSlice& vs2, const T_SpVec& sv)
{
	size_type split = vs1.dim();
	value_type result = 0;
	typename T_SpVec::difference_type offset = sv.offset();
	for(typename T_SpVec::const_iterator sv_itr = sv.begin(); sv_itr != sv.end(); ++sv_itr) {
		typename T_SpVec::element_type::indice_type curr_indice = sv_itr->index()+offset;
		if(curr_indice <= split)
			result += vs1(curr_indice) * sv_itr->value();
		else
			result += vs2(curr_indice - split) * sv_itr->value();
	}
	return result;
}

}	// end namespace SparseVectorUtilityPack

}	// end namespace SparseLinAlgPack

#endif // SPARSE_VECTOR_OP_DEF_H
