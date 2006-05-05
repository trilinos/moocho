// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <assert.h>

#include "AbstractLinAlgPack_MatrixConvertToSparseEncap.hpp"
#include "AbstractLinAlgPack_MatrixExtractSparseElements.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "DenseLinAlgPack_IVector.hpp"
#include "Teuchos_TestForException.hpp"

namespace AbstractLinAlgPack {

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
  ,const Range1D             &row_rng_in
  ,const i_vector_ptr_t      &inv_col_perm
  ,const Range1D             &col_rng_in
  ,const BLAS_Cpp::Transp    mese_trans
  ,const value_type          alpha
  )
{
  const size_type mese_rows = mese->rows(), mese_cols = mese->cols();
  const Range1D row_rng = RangePack::full_range(row_rng_in,1,mese_rows);
  const Range1D col_rng = RangePack::full_range(col_rng_in,1,mese_cols);
#ifdef _DEBUG
  const char msg_head[] = "MatrixConvertToSparseEncap::initialize(...): Error!";
  TEST_FOR_EXCEPTION( mese.get() == NULL, std::logic_error, msg_head );
  TEST_FOR_EXCEPTION( inv_row_perm.get() != NULL && inv_row_perm->size() != mese_rows, std::logic_error, msg_head );
  TEST_FOR_EXCEPTION( row_rng.ubound() > mese_rows, std::logic_error, msg_head );
  TEST_FOR_EXCEPTION( inv_col_perm.get() != NULL && inv_col_perm->size() != mese_cols, std::logic_error, msg_head );
  TEST_FOR_EXCEPTION( col_rng.ubound() > mese->cols(), std::logic_error, msg_head );
#endif
  mese_           = mese;
  mese_trans_     = mese_trans;
  alpha_          = alpha;
  inv_row_perm_   = inv_row_perm;
  inv_col_perm_   = inv_col_perm;
  row_rng_        = row_rng;
  col_rng_        = col_rng;
  nz_full_        = this->num_nonzeros(EXTRACT_FULL_MATRIX,ELEMENTS_ALLOW_DUPLICATES_SUM);
  space_cols_     = ( mese_trans_ == BLAS_Cpp::no_trans
            ? mese_->space_cols().sub_space(row_rng_)
            : mese_->space_rows().sub_space(col_rng_) );
  space_rows_     = ( mese_trans_ == BLAS_Cpp::no_trans
            ? mese_->space_rows().sub_space(col_rng_)
            : mese_->space_cols().sub_space(row_rng_) );
}

void MatrixConvertToSparseEncap::set_uninitialized()
{
  namespace mmp   = MemMngPack;
  mese_           = Teuchos::null;
  inv_row_perm_   = Teuchos::null;
  row_rng_        = Range1D::Invalid;
  inv_col_perm_   = Teuchos::null;
  col_rng_        = Range1D::Invalid;
  mese_trans_     = BLAS_Cpp::no_trans;
  alpha_          = 0.0;
  nz_full_        = 0;
}

// Overridden from MatrixBase

const VectorSpace& MatrixConvertToSparseEncap::space_cols() const
{
  return *space_cols_;
}

const VectorSpace& MatrixConvertToSparseEncap::space_rows() const
{
  return *space_rows_;
}

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

}	// end namespace AbstractLinAlgPack 
