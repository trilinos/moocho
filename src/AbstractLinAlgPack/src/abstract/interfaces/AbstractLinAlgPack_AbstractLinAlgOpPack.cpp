// ////////////////////////////////////////////////////////////////
// AbstractLinAlgOpPack.cpp
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

#include "AbstractLinAlgPack/src/LinAlgOpPackDef.h"

// Level 1 BLAS for Matrices

// M_lhs = op(M_rhs).

void LinAlgOpPack::assign(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs)
{
	Mp_M_assert_compatibility( M_lhs, BLAS_Cpp::no_trans, M_rhs, trans_rhs );
	M_lhs->zero_out();
	Mp_StM(M_lhs,1.0,M_rhs,trans_rhs);
}

// M_lhs = alpha * op(M_rhs).

void LinAlgOpPack::M_StM(MatrixWithOp* M_lhs, value_type alpha, const MatrixWithOp& M_rhs, BLAS_Cpp::Transp trans_rhs)
{
	Mp_M_assert_compatibility( M_lhs, BLAS_Cpp::no_trans, M_rhs, trans_rhs );
	M_lhs->zero_out();
	Mp_StM(M_lhs,alpha,M_rhs,trans_rhs);
}

// M_lhs = op(M_rhs1) + op(M_rhs2).

void LinAlgOpPack::M_MpM(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	Mp_M_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1);
	MopM_assert_compatibility(M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
	M_lhs->zero_out();
	Mp_M(M_lhs,M_rhs1,trans_rhs1);
	Mp_M(M_lhs,M_rhs2,trans_rhs2);
}

// M_lhs = op(M_rhs1) - op(M_rhs2).

void LinAlgOpPack::M_MmM(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	Mp_M_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1);
	MopM_assert_compatibility(M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
	M_lhs->zero_out();
	Mp_M(M_lhs,M_rhs1,trans_rhs1);
	Mp_StM(M_lhs,-1.0,M_rhs2,trans_rhs2);
}

// M_lhs = alpha * op(M_rhs1) + op(m_rhs2).

void LinAlgOpPack::M_StMpM(MatrixWithOp* M_lhs, value_type alpha, const MatrixWithOp& M_rhs1, BLAS_Cpp::Transp trans_rhs1
	, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	Mp_M_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1);
	MopM_assert_compatibility(M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
	assign(M_lhs,M_rhs2,trans_rhs2);
	Mp_StM(M_lhs,alpha,M_rhs1,trans_rhs1);
}

// Level 3 BLAS

// M_lhs = alpha * op(M_rhs1) * op(M_rhs2).

void LinAlgOpPack::M_StMtM(MatrixWithOp* M_lhs, value_type alpha, const MatrixWithOp& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	Mp_MtM_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
	Mp_StMtM(M_lhs,alpha,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2,0.0);
}

// M_lhs = op(M_rhs1) * op(M_rhs2).

void LinAlgOpPack::M_MtM(MatrixWithOp* M_lhs, const MatrixWithOp& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const MatrixWithOp& M_rhs2, BLAS_Cpp::Transp trans_rhs2)
{
	Mp_MtM_assert_compatibility(M_lhs,BLAS_Cpp::no_trans,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2);
	Mp_StMtM(M_lhs,1.0,M_rhs1,trans_rhs1,M_rhs2,trans_rhs2,0.0);
}
