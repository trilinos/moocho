// /////////////////////////////////////////////////////////////////////
// MatrixSymNonsingularSerial.cpp
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

#include <assert.h>

#include "AbstractLinAlgPack/src/serial/interfaces/MatrixSymNonsingularSerial.hpp"
#include "AbstractLinAlgPack/src/serial/interfaces/MatrixSymWithOpGetGMSSymMutable.hpp"
#include "AbstractLinAlgPack/src/serial/interfaces/MatrixWithOpSerial.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/EtaVector.hpp"
#include "DenseLinAlgPack/src/DMatrixClass.hpp"
#include "DenseLinAlgPack/src/DMatrixOp.hpp"
#include "DenseLinAlgPack/src/DMatrixAsTriSym.hpp"
#include "DenseLinAlgPack/src/LinAlgOpPack.hpp"
#include "DenseLinAlgPack/src/DenseLinAlgPackAssertOp.hpp"
#include "dynamic_cast_verbose.hpp"

namespace LinAlgOpPack {
	using AbstractLinAlgPack::Vp_StMtV;
	using AbstractLinAlgPack::Mp_StMtM;
}

namespace AbstractLinAlgPack {

void MatrixSymNonsingularSerial::M_StMtInvMtM(
	  DMatrixSliceSym* S, value_type a, const MatrixWithOpSerial& B
	, BLAS_Cpp::Transp B_trans, EMatrixDummyArg ) const
{
	using BLAS_Cpp::trans;
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans_not;
	using AbstractLinAlgPack::M_StInvMtM;
	using DenseLinAlgPack::nonconst_tri_ele;
	using DenseLinAlgPack::tri_ele;
	using DenseLinAlgPack::assign;
	using LinAlgOpPack::M_StMtM;
	//
	// S = a * op(B) * inv(M) * op(B')
	//
	// We will form S won column at a time:
	//
	// S(:,j) = a * op(B) * inv(M) * op(B') * e(j)
	//
	// for j = 1 ... op(B').cols()
	//   t1 = op(B')*e(j)
	//   t2 = inv(M)*t1
	//   t3 = a*op(B)*t2
	//   S(:,j) = t3
	//
	// Above we only need to set the lower (lower triangle stored)
	// or upper (upper triangle stored) part of S(:,k)
	//
	DenseLinAlgPack::MtM_assert_sizes( rows(), cols(), no_trans
		, B.rows(), B.cols(), trans_not(B_trans) );
	DenseLinAlgPack::Mp_MtM_assert_sizes( S->rows(), S->cols(), no_trans
		, B.rows(), B.cols(), B_trans
		, B.rows(), B.cols(), trans_not(B_trans) );

	DVector t1, t2, t3; // ToDo: Use temp workspace!
	const size_type
		opBT_cols = BLAS_Cpp::cols( B.cols(), B.rows(), B_trans ),
		m         = S->rows();
	for( size_type j = 1; j <= m; ++j ) {
		EtaVector e_j(j,opBT_cols);                               // e(j)
		LinAlgOpPack::V_MtV( &t1, B, trans_not(B_trans), e_j() ); // t1 = op(B')*e(j)
		AbstractLinAlgPack::V_InvMtV( &t2, *this, no_trans, t1() ); // t2 = inv(M)*t1
		LinAlgOpPack::V_StMtV( &t3, a, B, B_trans, t2() );        // t3 = a*op(B)*t2
		Range1D
			rng = ( S->uplo() == BLAS_Cpp::upper ? Range1D(1,j) : Range1D(j,m) );
		S->gms().col(j)(rng) = t3(rng);
	}
}

// Overridden from MatrixSymNonsing

void MatrixSymNonsingularSerial::M_StMtInvMtM(
	MatrixSymOp* symwo_lhs, value_type alpha
	,const MatrixOp& mwo, BLAS_Cpp::Transp mwo_trans
	,EMatrixDummyArg dummy
	) const
{
	using DynamicCastHelperPack::dyn_cast;
	this->M_StMtInvMtM(
		&MatrixDenseSymMutableEncap(symwo_lhs)(), alpha
		,dyn_cast<const MatrixWithOpSerial>(mwo), mwo_trans
		,dummy );
}

} // end namespace AbstractLinAlgPack
