// /////////////////////////////////////////////////////////////////
// MatrixSymOpSerial.cpp
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

#include "AbstractLinAlgPack/src/serial/interfaces/MatrixSymOpSerial.hpp"
#include "AbstractLinAlgPack/src/serial/interfaces/MatrixSymOpGetGMSSymMutable.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/GenPermMatrixSlice.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/EtaVector.hpp"
#include "DenseLinAlgPack/src/DMatrixOp.hpp"
#include "DenseLinAlgPack/src/DMatrixAsTriSym.hpp"
#include "DenseLinAlgPack/src/LinAlgOpPack.hpp"
#include "DenseLinAlgPack/src/DenseLinAlgPackAssertOp.hpp"
#include "dynamic_cast_verbose.hpp"

namespace AbstractLinAlgPack {

void MatrixSymOpSerial::Mp_StPtMtP(
	DMatrixSliceSym* S, value_type a
	,EMatRhsPlaceHolder
	,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	,value_type b
	) const
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
	using BLAS_Cpp::cols;
	//
	// S = b*S
	//
	// S += a*op(P')*M*op(P)
	//
	// We will perform this operation for each column in op(P) as:
	//
	// op(P)(:,j(k)) = e(i(k)) <: R^n
	//
	// S += a*op(P')*M*[ ... e(i(1)) ... e(i(k)) ... e(i(nz)) ... ]
	//                                     j(k)
	//
	// We will perform this by column as:
	//
	// for k = 1...nz
	//    S(:,j(k)) += a*y_k
	//
	//    where:
	//       y_k = op(P')*M*e(i(k))
	//
	// Above we only need to set the portion of S(:,j(k)) for the stored part
	// of the symmetric matrix (i.e. upper part for upper and lower part for lower).
	//
	DenseLinAlgPack::MtM_assert_sizes(
		this->rows(), this->cols(), no_trans
		, P.rows(), P.cols(), P_trans );
	DenseLinAlgPack::Mp_M_assert_sizes(
		S->rows(), S->cols(), no_trans
		, cols( P.rows(), P.cols(), P_trans )
		, cols( P.rows(), P.cols(), P_trans )
		, no_trans );
	//
	const size_type
		n = this->rows(),
		m = S->rows();
	// S = b*S
	if( b != 1.0 )
		DenseLinAlgPack::Mt_S( &DenseLinAlgPack::nonconst_tri_ele(S->gms(),S->uplo()), b );
	// Set the colums of S
	DVector y_k_store(m);
	DVectorSlice y_k = y_k_store();
	for( GenPermMatrixSlice::const_iterator P_itr = P.begin(); P_itr != P.end(); ++P_itr )
	{
		const size_type
			i_k = P_trans == no_trans ? P_itr->row_i() : P_itr->col_j(),
			j_k = P_trans == no_trans ? P_itr->col_j() : P_itr->row_i();
		// e(i(k))
		EtaVector
			e_i_k(i_k,n);
		// y_k = op(P')*M*e(i(k))
		AbstractLinAlgPack::Vp_StPtMtV( &y_k, 1.0, P, trans_not(P_trans), *this, no_trans, e_i_k(), 0.0 );
		// S(:,j(k)) += a*y_k
		if( S->uplo() == BLAS_Cpp::upper )
			DenseLinAlgPack::Vp_StV( &S->gms().col(j_k)(1,j_k), a, y_k(1,j_k) );
		else
			DenseLinAlgPack::Vp_StV( &S->gms().col(j_k)(j_k,m), a, y_k(j_k,m) );
	}
}

void MatrixSymOpSerial::Mp_StMtMtM(
	DMatrixSliceSym* sym_lhs, value_type alpha
	,EMatRhsPlaceHolder dummy_place_holder
	,const MatrixOpSerial& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
	,value_type beta
	) const
{
	assert(0);	// ToDo: Give this some default implementation for 
				// this at some point in the future?
}

// Overridden from MatrixSymOp

const VectorSpace& MatrixSymOpSerial::space_rows() const
{
	return MatrixOpSerial::space_rows();
}

void MatrixSymOpSerial::Mp_StPtMtP(
	MatrixSymOp* symwo_lhs, value_type alpha
	,EMatRhsPlaceHolder dummy
	,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
	,value_type beta
	) const
{
	this->Mp_StPtMtP(
		&MatrixDenseSymMutableEncap(symwo_lhs)(), alpha, dummy
		,gpms_rhs, gpms_rhs_trans
		, beta );
}

void MatrixSymOpSerial::Mp_StMtMtM(
	MatrixSymOp* symwo_lhs, value_type alpha
	,EMatRhsPlaceHolder dummy
	,const MatrixOp& mwo_rhs, BLAS_Cpp::Transp mwo_rhs_trans
	,value_type beta
	) const
{
	using DynamicCastHelperPack::dyn_cast;
	this->Mp_StMtMtM(
		&MatrixDenseSymMutableEncap(symwo_lhs)(), alpha, dummy
		,dyn_cast<const MatrixOpSerial>(mwo_rhs), mwo_rhs_trans
		,beta );
}

}	// end namespace AbstractLinAlgPack 
