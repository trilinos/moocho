// ///////////////////////////////////////////////////////////////////////////
// MatrixVectorTemplateOpDef.hpp
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
// Definitions of template functions declared in MatrixVectorTemplateOp.hpp.

#ifndef MATRIX_VECTOR_TEMPLATE_OP_DEF_H
#define MATRIX_VECTOR_TEMPLATE_OP_DEF_H

#include "MatrixVectorTemplateOp.hpp"
#include "DenseLinAlgPack/src/DMatrixClass.hpp"

// ///////////////////////////////////
// Matrix assignment

namespace {
// Typedef for vector returning functions (row or col but not diagonal)
typedef AbstractLinAlgPack::DVectorSlice (AbstractLinAlgPack::DMatrixSlice::*Pvec_func)
	(AbstractLinAlgPack::DMatrixSlice::size_type);
// Implement sparse matrix, dense matrix assignment.  Sizes are not checked.
template<class T_Matrix>
void imp_assign(AbstractLinAlgPack::DMatrixSlice& gms_lhs, const T_Matrix& gm_rhs
	, BLAS_Cpp::Transp trans_rhs)
{
	// If trans, then copy col into row, otherwise copy col into col.
	Pvec_func vec_func;
	if(trans_rhs == BLAS_Cpp::no_trans)		vec_func = &AbstractLinAlgPack::DMatrixSlice::col;
	else									vec_func = &AbstractLinAlgPack::DMatrixSlice::row;
	for(int k = 1; k <= gm_rhs.cols(); ++k)
		AbstractLinAlgPack::assign((gms_lhs.*vec_func)(k), gm_rhs.col(k));
}
} // end namespace

// Definitions of template functions for matrix-matrix assignment

/// gm_lhs = T_M (templated matrix type T_M)
template<class T_Matrix>
void AbstractLinAlgPack::assign(DMatrix& gm_lhs, const T_Matrix& gm_rhs
	, BLAS_Cpp::Transp trans_rhs)
{
	DenseLinAlgPack::resize_gm_lhs(gm_lhs,gm_rhs.rows(),gm_rhs.cols(),trans_rhs);
	DMatrixSlice gms_lhs = gm_lhs;
	imp_assign(gms_lhs,gm_rhs,trans_rhs);
}

/// gms_lhs = T_M (templated matrix type T_M)
template<class T_Matrix>
void AbstractLinAlgPack::assign(DMatrixSlice& gms_lhs, const T_Matrix& gm_rhs
	, BLAS_Cpp::Transp trans_rhs)
{
	DenseLinAlgPack::assert_gms_lhs(gms_lhs,gm_rhs.rows(),gm_rhs.cols(),trans_rhs);
	imp_assign(gms_lhs,gm_rhs,trans_rhs);
}

// //////////////////////////////////
// Matrix-DVector multiplication

namespace {
// Throw an exeption of the rhs arguments do not match
template<class T_Matrix>
void imp_assert_V_MtV_rhs_sizes(const T_Matrix& gm_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const AbstractLinAlgPack::DVectorSlice& vs_rhs2)
{
	typename T_Matrix::size_type
		cols = (trans_rhs1 == BLAS_Cpp::no_trans) ? gm_rhs1.cols() : gm_rhs1.rows();

	if(cols != vs_rhs2.size())
		throw std::length_error("V_MtV: The sizes of the rhs expression do not match");
}

// Implementation of matrix-vector multiply (no transpose).  Sizes are not checked.
template<class T_Matrix>
void imp_V_MtV_no_trans(AbstractLinAlgPack::DVectorSlice& vs_lhs, const T_Matrix& gm_rhs1
	, const AbstractLinAlgPack::DVectorSlice& vs_rhs2)
{
	typedef typename T_Matrix::size_type size_type;
	size_type rows = gm_rhs1.rows();
	AbstractLinAlgPack::DVectorSlice::iterator itr_v_lhs = vs_lhs.begin();
	for(size_type i = 1; i <= rows; ++i)
		*itr_v_lhs++ = AbstractLinAlgPack::dot(vs_rhs2,gm_rhs1.row(i));
}
// Implementation of matrix-vector multiply (transpose).  Sizes are not checked.
template<class T_Matrix>
void imp_V_MtV_trans(AbstractLinAlgPack::DVectorSlice& vs_lhs, const T_Matrix& gm_rhs1
	, const AbstractLinAlgPack::DVectorSlice& vs_rhs2)
{
	typedef typename T_Matrix::size_type size_type;
	size_type cols = gm_rhs1.cols();
	AbstractLinAlgPack::DVectorSlice::iterator itr_v_lhs = vs_lhs.begin();
	for(size_type j = 1; j <= cols; ++j)
		*itr_v_lhs++ = AbstractLinAlgPack::dot(vs_rhs2,gm_rhs1.col(j));
}

} // end namespace

// Definitions of template functions for matrix-vector multiplication

template<class T_Matrix>
void AbstractLinAlgPack::V_MtV(DVector& v_lhs, const T_Matrix& gm_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const DVectorSlice& vs_rhs2)
{
	imp_assert_V_MtV_rhs_sizes(gm_rhs1,trans_rhs1,vs_rhs2);
	v_lhs.resize( (trans_rhs1==BLAS_Cpp::no_trans) ? gm_rhs1.rows() : gm_rhs1.cols() );
	DVectorSlice vs_lhs = v_lhs;
	if(trans_rhs1 == BLAS_Cpp::no_trans)
		imp_V_MtV_no_trans(vs_lhs,gm_rhs1,vs_rhs2);
	else
		imp_V_MtV_trans(vs_lhs,gm_rhs1,vs_rhs2);
}

template<class T_Matrix>
void AbstractLinAlgPack::V_MtV(DVectorSlice& v_lhs, const T_Matrix& gm_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const DVectorSlice& vs_rhs2)
{
	imp_assert_V_MtV_rhs_sizes(gm_rhs1,trans_rhs1,vs_rhs2);
	DenseLinAlgPack::assert_resize_vs_lhs(v_lhs, (trans_rhs1==BLAS_Cpp::no_trans) ? gm_rhs1.rows() : gm_rhs1.cols());
	if(trans_rhs1 == BLAS_Cpp::no_trans)
		imp_V_MtV_no_trans(v_lhs,gm_rhs1,vs_rhs2);
	else
		imp_V_MtV_trans(v_lhs,gm_rhs1,vs_rhs2);
}

#endif // MATRIX_VECTOR_TEMPLATE_OP_DEF_H
