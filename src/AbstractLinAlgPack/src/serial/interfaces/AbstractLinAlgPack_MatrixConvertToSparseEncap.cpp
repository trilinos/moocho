// //////////////////////////////////////////////////////////////////////////////
// MatrixConvertToSparseEncap.cpp
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

#include "SparseLinAlgPack/include/MatrixConvertToSparseEncap.h"
#include "SparseLinAlgPack/include/MatrixExtractSparseElements.h"
#include "LinAlgPack/include/IVector.h"
#include "ThrowException.h"

namespace SparseLinAlgPack {

// Constructors/initializers

MatrixConvertToSparseEncap::MatrixConvertToSparseEncap()
	:row_rng_(Range1D::Invalid)
	,col_rng_(Range1D::Invalid)
	,mese_trans_(BLAS_Cpp::no_trans)
	,alpha_(0.0)
	,nz_full_(0)
{}

MatrixConvertToSparseEncap::MatrixConvertToSparseEncap(
	const mese_ptr_t           &mese
	,const i_vector_ptr_t      &inv_row_perm
	,const Range1D             &row_rng
	,const i_vector_ptr_t      &inv_col_perm
	,const Range1D             &col_rng
	,const BLAS_Cpp::Transp    mese_trans
	,const value_type          alpha
	)
{
	this->initialize(mese,inv_row_perm,row_rng,inv_col_perm,col_rng,mese_trans,alpha);
}

void MatrixConvertToSparseEncap::initialize(
	const mese_ptr_t           &mese
	,const i_vector_ptr_t      &inv_row_perm
	,const Range1D             &row_rng
	,const i_vector_ptr_t      &inv_col_perm
	,const Range1D             &col_rng
	,const BLAS_Cpp::Transp    mese_trans
	,const value_type          alpha
	)
{
#ifdef _DEBUG
	const char msg_head[] = "MatrixConvertToSparseEncap::initialize(...): Error!";
	THROW_EXCEPTION( mese.get() == NULL, std::logic_error, msg_head );
	const size_type mese_rows = mese->rows(), mese_cols = mese->cols();
	THROW_EXCEPTION( !(inv_row_perm.get() == NULL || inv_row_perm->size() == mese_rows), std::logic_error, msg_head );
	THROW_EXCEPTION( row_rng.ubound() > mese_rows, std::logic_error, msg_head );
	THROW_EXCEPTION( !(inv_col_perm.get() == NULL || inv_col_perm->size() == mese_cols), std::logic_error, msg_head );
	THROW_EXCEPTION( col_rng.ubound() > mese->cols(), std::logic_error, msg_head );
#endif
	mese_           = mese;
	inv_row_perm_   = inv_row_perm;
	row_rng_        = row_rng;
	inv_col_perm_   = inv_col_perm;
	col_rng_        = col_rng;
	mese_trans_     = mese_trans_;
	alpha_          = alpha;
	nz_full_        = this->num_nonzeros(EXTRACT_FULL_MATRIX,ELEMENTS_ALLOW_DUPLICATES_SUM);
}

void MatrixConvertToSparseEncap::set_uninitialized()
{
	namespace mmp   = MemMngPack;
	mese_           = mmp::null;
	inv_row_perm_   = mmp::null;
	row_rng_        = Range1D::Invalid;
	inv_col_perm_   = mmp::null;
	col_rng_        = Range1D::Invalid;
	mese_trans_     = BLAS_Cpp::no_trans;
	alpha_          = 0.0;
	nz_full_        = 0;
}

// Overridden from MatrixBase

size_type MatrixConvertToSparseEncap::rows() const
{
	return mese_trans_ == BLAS_Cpp::no_trans ? row_rng_.size() : col_rng_.size();
}

size_type MatrixConvertToSparseEncap::cols() const
{
	return mese_trans_ == BLAS_Cpp::no_trans ? col_rng_.size() : row_rng_.size();
}

size_type MatrixConvertToSparseEncap::nz() const
{
	return nz_full_;
}

// Overridden from MatrixConvertToSparse

index_type MatrixConvertToSparseEncap::num_nonzeros(
	EExtractRegion        extract_region_in
	,EElementUniqueness   element_uniqueness
	) const
{
	index_type dl = 0, du = 0;
	EExtractRegion extract_region = extract_region_in;
	if( mese_trans_ == BLAS_Cpp::trans )
		extract_region
			= ( extract_region_in == EXTRACT_FULL_MATRIX
				? EXTRACT_FULL_MATRIX
				: ( extract_region_in == EXTRACT_UPPER_TRIANGULAR
					? EXTRACT_LOWER_TRIANGULAR
					: EXTRACT_UPPER_TRIANGULAR ) );
	switch(extract_region) {
		case EXTRACT_FULL_MATRIX:
			dl = -(index_type)row_rng_.ubound() + (index_type)col_rng_.lbound();
			du = +(index_type)col_rng_.ubound() - (index_type)row_rng_.lbound();
			break;
		case EXTRACT_UPPER_TRIANGULAR:
			dl = -(index_type)row_rng_.lbound() + (index_type)col_rng_.lbound();
			du = +(index_type)col_rng_.ubound() - (index_type)row_rng_.lbound();
			break;
		case EXTRACT_LOWER_TRIANGULAR:
			dl = -(index_type)row_rng_.ubound() + (index_type)col_rng_.lbound();
			du = +(index_type)col_rng_.lbound() - (index_type)row_rng_.lbound();
			break;
		default:
			assert(0);
	}
	const index_type
		*inv_row_perm = inv_row_perm_.get() ? &(*inv_row_perm_)(1) : NULL,
		*inv_col_perm = inv_col_perm_.get() ? &(*inv_col_perm_)(1) : NULL;
	return mese_->count_nonzeros(
		element_uniqueness
		,inv_row_perm
		,inv_col_perm
		,row_rng_
		,col_rng_
		,dl
		,du
		);
}

void MatrixConvertToSparseEncap::coor_extract_nonzeros(
	EExtractRegion                extract_region
	,EElementUniqueness           element_uniqueness
	,const index_type             len_Aval
	,value_type                   Aval[]
	,const index_type             len_Aij
	,index_type                   Arow[]
	,index_type                   Acol[]
	,const index_type             row_offset
	,const index_type             col_offset
	) const
{
	index_type dl = 0, du = 0;  // This may not be correct!
	switch(extract_region) {
		case EXTRACT_FULL_MATRIX:
			dl = -(index_type)row_rng_.ubound() + (index_type)col_rng_.lbound();
			du = +(index_type)col_rng_.ubound() - (index_type)row_rng_.lbound();
			break;
		case EXTRACT_UPPER_TRIANGULAR:
			dl = -(index_type)row_rng_.lbound() + (index_type)col_rng_.lbound();
			du = +(index_type)col_rng_.ubound() - (index_type)row_rng_.lbound();
			break;
		case EXTRACT_LOWER_TRIANGULAR:
			dl = -(index_type)row_rng_.ubound() + (index_type)col_rng_.lbound();
			du = +(index_type)col_rng_.lbound() - (index_type)row_rng_.lbound();
			break;
		default:
			assert(0);
	}
	const index_type
		*inv_row_perm = inv_row_perm_.get() ? &(*inv_row_perm_)(1) : NULL,
		*inv_col_perm = inv_col_perm_.get() ? &(*inv_col_perm_)(1) : NULL;
	mese_->coor_extract_nonzeros(
		element_uniqueness
		,inv_row_perm
		,inv_col_perm
		,row_rng_
		,col_rng_
		,dl
		,du
		,alpha_
		,len_Aval
		,Aval
		,len_Aij
		,mese_trans_ == BLAS_Cpp::no_trans ? Arow : Acol
		,mese_trans_ == BLAS_Cpp::no_trans ? Acol : Arow
		,mese_trans_ == BLAS_Cpp::no_trans ? row_offset : col_offset
		,mese_trans_ == BLAS_Cpp::no_trans ? col_offset : row_offset
		);
}

}	// end namespace SparseLinAlgPack 