// /////////////////////////////////////////////////////////////////////
// MatrixSparseCOORSerial.cpp
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

#include "SparseLinAlgPack/src/MatrixSparseCOORSerial.h"
#include "SparseLinAlgPack/src/MatrixCOORTmplItfc.h"
#include "SparseLinAlgPack/src/COOMatrixTmplOp.h"
#include "SparseLinAlgPack/src/VectorDenseEncap.h"
#include "AbstractLinAlgPack/src/AbstractLinAlgPackAssertOp.h"
#include "LinAlgPack/src/LinAlgOpPack.h"
#include "ThrowException.h"
#include "dynamic_cast_verbose.h"

namespace SparseLinAlgPack {

MatrixSparseCOORSerial::ReleaseValRowColArrays::~ReleaseValRowColArrays()
{
	if(owns_mem_) {
		if(val_)   delete [] val_;
		if(row_i_) delete [] row_i_;
		if(col_j_) delete [] col_j_;
	}
}

bool MatrixSparseCOORSerial::ReleaseValRowColArrays::resource_is_bound() const
{
	return val_ != NULL;
}

// static members

MatrixSparseCOORSerial::release_resource_ptr_t
MatrixSparseCOORSerial::release_resource_null_ = MemMngPack::null;

// Constructors / initializers

MatrixSparseCOORSerial::MatrixSparseCOORSerial()
	:rows_(0)
	,cols_(0)
	,max_nz_(0)
	,nz_(0)
	,val_(NULL)
	,row_i_(NULL)
	,col_j_(NULL)
	,self_allocate_(true)
{}

void MatrixSparseCOORSerial::set_buffers(
	size_type                      max_nz
	,value_type                    *val
	,index_type                    *row_i
	,index_type                    *col_j
	,const release_resource_ptr_t  &release_resource
	,size_type                     rows
	,size_type                     cols
	,size_type                     nz
	,bool                          check_input
	)
{
#ifdef _DEBUG
	const char msg_err[] = "MatrixSparseCOORSerial::set_buffer(...) : Error,!";
	THROW_EXCEPTION( max_nz <= 0, std::invalid_argument, msg_err );
	THROW_EXCEPTION( val == NULL || row_i == NULL || col_j == NULL, std::invalid_argument, msg_err );
	THROW_EXCEPTION( rows > 0 && cols <= 0 , std::invalid_argument, msg_err );
	THROW_EXCEPTION( rows > 0 && (nz < 0 || nz > max_nz), std::invalid_argument, msg_err );
#endif
	max_nz_           = max_nz;
	val_              = val;
	row_i_            = row_i;
	col_j_            = col_j;
	release_resource_ = release_resource;
	self_allocate_    = false;
	if(rows) {
		rows_ = rows;
		cols_ = cols;
		nz_   = nz;
		space_cols_.initialize(rows);
		space_rows_.initialize(cols);
		if( nz && check_input ) {
			assert(0); // Todo: Check that row_i[] and col_j[] are in bounds
		}
	}
	else {
		rows_ = 0;
		cols_ = 0;
		nz_   = 0;
		space_cols_.initialize(0);
		space_rows_.initialize(0);
	}
}

void MatrixSparseCOORSerial::set_uninitialized()
{
	max_nz_           = 0;
	val_              = NULL;
	row_i_            = NULL;
	col_j_            = NULL;
	release_resource_ = MemMngPack::null;
	self_allocate_    = true;
	rows_             = 0;
	cols_             = 0;
	nz_               = 0;
	space_cols_.initialize(0);
	space_rows_.initialize(0);
}

// Overridden from MatrixBase

size_type MatrixSparseCOORSerial::rows() const
{
	return rows_;
}

size_type MatrixSparseCOORSerial::cols() const
{
	return cols_;
}

size_type MatrixSparseCOORSerial::nz() const
{
	return nz_;
}

// Overridden from MatrixWithOp

const VectorSpace& MatrixSparseCOORSerial::space_cols() const
{
	return space_cols_;
}

const VectorSpace& MatrixSparseCOORSerial::space_rows() const
{
	return space_rows_;
}

MatrixWithOp& MatrixSparseCOORSerial::operator=(const MatrixWithOp& M)
{
	using DynamicCastHelperPack::dyn_cast;
	const MatrixSparseCOORSerial
		&Mc = dyn_cast<const MatrixSparseCOORSerial>(M);
	if( this == &Mc )
		return *this; // assignment to self
	// A shallow copy is fine as long as we are carefull.
	max_nz_           = Mc.max_nz_;
	val_              = Mc.val_;
	row_i_            = Mc.row_i_;
	col_j_            = Mc.col_j_;
	release_resource_ = Mc.release_resource_;
	self_allocate_    = Mc.self_allocate_;
	rows_             = Mc.rows_;
	cols_             = Mc.cols_;
	nz_               = Mc.nz_;
	space_cols_.initialize(rows_);
	space_rows_.initialize(cols_);
	return *this;
}

std::ostream& MatrixSparseCOORSerial::output(std::ostream& out) const
{
	const size_type
		rows = this->rows(),
		cols = this->cols(),
		nz   = this->nz();
	out
		<< "Sparse " << rows << " x " << cols << " matrix with "
		<< nz << " nonzero entries:\n";
	const value_type
		*itr_val      = &val_[0],
		*itr_val_end  = itr_val + nz_;
	const index_type
		*itr_ivect    = &row_i_[0],
		*itr_jvect    = &col_j_[0];
	for(; itr_val != itr_val_end; ++itr_val, ++itr_ivect, ++itr_jvect)
		out << *itr_val << ":" << *itr_ivect << ":" << *itr_jvect << " ";
	out << "\n";
	
	if( rows * cols <= 400 ) {
		out << "Converted to dense =\n";
		MatrixWithOp::output(out);
	}

	return out;
}

void MatrixSparseCOORSerial::Vp_StMtV(
	VectorWithOpMutable* y, value_type a, BLAS_Cpp::Transp M_trans
	, const VectorWithOp& x, value_type b
	) const
{
	AbstractLinAlgPack::Vp_MtV_assert_compatibility( y, *this, M_trans, x );
	VectorDenseMutableEncap   y_d(*y);
	VectorDenseEncap          x_d(x);
	LinAlgPack::Vt_S(&y_d(),b);
	Vp_StCOOMtV(
		&y_d(),a,MatrixCOORTmplItfc<value_type,index_type>(
			rows_,cols_,nz_,0,0,val_,row_i_,col_j_ )
		,M_trans, x_d()
		);
}

// Overridden from MatrixLoadSparseElements

void MatrixSparseCOORSerial::reinitialize(
	size_type                  rows
	,size_type                 cols
	,size_type                 max_nz
	,EAssumeElementUniqueness  element_uniqueness
	)
{
	namespace rcp = MemMngPack;
#ifdef _DEBUG
	const char msg_err_head[] = "MatrixSparseCOORSerial::reinitialize(...) : Error";
	THROW_EXCEPTION( max_nz <= 0, std::invalid_argument, msg_err_head<<"!" );
	THROW_EXCEPTION( rows <= 0 || cols <= 0 , std::invalid_argument, msg_err_head<<"!" );
#endif
	rows_               = rows;
	cols_               = cols;	
	element_uniqueness_ = element_uniqueness;
	if( self_allocate_ ) {
		if(max_nz_ < max_nz) {
			release_resource_ = rcp::rcp(
				new ReleaseValRowColArrays(
					val_    = new value_type[max_nz]
					,row_i_ = new index_type[max_nz]
					,col_j_ = new index_type[max_nz]
					)
				);
			max_nz_ = max_nz;
		}

	}
	else {
#ifdef _DEBUG
		THROW_EXCEPTION(
			max_nz <= max_nz_, std::invalid_argument
			,msg_err_head << "Buffers set up by client in set_buffers() only allows storage for "
			"max_nz_ = " << max_nz_ << " nonzero entries while client requests storage for "
			"max_nz = " << max_nz << " nonzero entries!" );
#endif
	}
	reload_val_only_         = false;
	reload_val_only_nz_last_ = 0;
	nz_                      = 0;
	max_nz_load_             = 0;
}

void MatrixSparseCOORSerial::reset_to_load_values()
{
#ifdef _DEBUG
	THROW_EXCEPTION(
		rows_ == 0 || cols_ == 0, std::invalid_argument
		,"MatrixSparseCOORSerial::reset_to_load_values(...) : Error, "
		"this matrix is not initialized so it can't be rest to load "
		"new values for nonzero entries." );
#endif
	reload_val_only_         = true;
	reload_val_only_nz_last_ = nz_;
	nz_                      = 0;
	max_nz_load_             = 0;
}

void MatrixSparseCOORSerial::get_load_nonzeros_buffers(
	size_type      max_nz_load
	,value_type    **val
	,index_type    **row_i
	,index_type    **col_j
	)
{
#ifdef _DEBUG
	THROW_EXCEPTION(
		max_nz_load_ != 0 , std::logic_error
		,"MatrixSparseCOORSerial::get_load_nonzeros_buffers(...) : Error, "
		"You must call commit_load_nonzeros_buffers() between calls to this method!" );
	THROW_EXCEPTION(
		max_nz_load <= 0 || max_nz_load > max_nz_ - nz_, std::invalid_argument
		,"MatrixSparseCOORSerial::get_load_nonzeros_buffers(...) : Error, "
		"The number of nonzeros to load max_nz_load = " << max_nz_load << " can not "
		"be greater than max_nz - nz = " << max_nz_ << " - " << nz_ << " = " << (max_nz_-nz_) <<
		" entries!" );
	THROW_EXCEPTION(
		reload_val_only_ && (row_i != NULL || col_j != NULL), std::invalid_argument
		,"MatrixSparseCOORSerial::get_load_nonzeros_buffers(...) : Error, "
		"reset_to_load_values() was called and therefore the structure of the matrix "
		"can not be set!" );
	THROW_EXCEPTION(
		!reload_val_only_  && (row_i == NULL || col_j == NULL), std::invalid_argument
		,"MatrixSparseCOORSerial::get_load_nonzeros_buffers(...) : Error, "
		"both *row_i and *col_j must be non-NULL since reinitialize() was called" );
#endif
	max_nz_load_ = max_nz_load;
	*val   = val_   + nz_;
	if(!reload_val_only_)
		*row_i = row_i_ + nz_;
	if(!reload_val_only_)
		*col_j = col_j_ + nz_;
}

void MatrixSparseCOORSerial::commit_load_nonzeros_buffers(
	size_type      nz_commit
	,value_type    **val
	,index_type    **row_i
	,index_type    **col_j
	)
{
#ifdef _DEBUG
	THROW_EXCEPTION(
		max_nz_load_ == 0 , std::logic_error
		,"MatrixSparseCOORSerial::commit_load_nonzeros_buffers(...) : Error, "
		"You must call get_load_nonzeros_buffers() before calling this method!" );
	THROW_EXCEPTION(
		nz_commit > max_nz_load_ , std::logic_error
		,"MatrixSparseCOORSerial::commit_load_nonzeros_buffers(...) : Error, "
		"You can not commit more nonzero entries than you requested buffer space for in "
		"get_load_nonzeros_buffers(...)!" );
	THROW_EXCEPTION(
		*val != val_ + nz_
		, std::logic_error
		,"MatrixSparseCOORSerial::commit_load_nonzeros_buffers(...) : Error, "
		"This is not the buffer I give you in get_load_nonzeros_buffers(...)!" );
	THROW_EXCEPTION(
		reload_val_only_ && (row_i != NULL || col_j != NULL), std::invalid_argument
		,"MatrixSparseCOORSerial::commit_load_nonzeros_buffers(...) : Error, "
		"reset_to_load_values() was called and therefore the structure of the matrix "
		"can not be set!" );
#endif
	nz_            += nz_commit;
	max_nz_load_    = 0;
}

void MatrixSparseCOORSerial::finish_construction( bool test_setup )
{
	THROW_EXCEPTION(
		reload_val_only_ == true && reload_val_only_nz_last_ != nz_, std::logic_error
		,"MatrixSparseCOORSerial::finish_construction() : Error, the number of nonzeros on"
		" the initial load with row and column indexes was = " << reload_val_only_nz_last_ <<
		" and does not agree with the number of nonzero values = " << nz_ << " loaded this time!" );
	if( test_setup ) {
		for( size_type k = 0; k < nz_; ++k ) {
			const index_type
				i = row_i_[k],
				j = col_j_[k];
			THROW_EXCEPTION(
				i < 1 || rows_ < i, std::logic_error
				,"MatrixSparseCOORSerial::finish_construction(true) : Error, "
				"row_i[" << k << "] = " << i << " is not in the range [1,rows] = [1,"<<rows_<<"]!" );
			THROW_EXCEPTION(
				j < 1 || cols_ < j, std::logic_error
				,"MatrixSparseCOORSerial::finish_construction(true) : Error, "
				"col_j[" << k << "] = " << j << " is not in the range [1,cols] = [1,"<<cols_<<"]!" );
		}
	}
	space_cols_.initialize(rows_);
	space_rows_.initialize(cols_);
}

// Overridden from MatrixExtractSparseElements

#ifdef _DEBUG
#define VALIDATE_ROW_COL_IN_RANGE() \
THROW_EXCEPTION( \
	i < 1 || rows_ < i, std::invalid_argument \
	,err_msg_head<<", i = inv_row_perm[(row_i["<<k<<"]=="<<*row_i<<")-1] = "<<i<<" > rows = "<<rows_ ); \
THROW_EXCEPTION( \
	j < 1 || cols_ < j, std::invalid_argument \
	,err_msg_head<<", j = inv_col_perm[(col_j["<<k<<"]=="<<*col_j<<")-1] = "<<j<<" > rows = "<<cols_ );
#else
#define VALIDATE_ROW_COL_IN_RANGE()
#endif

index_type MatrixSparseCOORSerial::count_nonzeros(
	EElementUniqueness    element_uniqueness
	,const index_type     inv_row_perm[]
	,const index_type     inv_col_perm[]
	,const Range1D        &row_rng_in
	,const Range1D        &col_rng_in
	,index_type           dl
	,index_type           du
	) const
{
#ifdef _DEBUG
	const char err_msg_head[] = "MatrixSparseCOORSerial::count_nonzeros(...): Error";
	THROW_EXCEPTION(
		element_uniqueness_ == ELEMENTS_ASSUME_DUPLICATES_SUM && element_uniqueness == ELEMENTS_FORCE_UNIQUE
		,std::logic_error
		,err_msg_head << ", the client requests a count for unique "
		"elements but this sparse matrix object is not allowed to assume this!" );
#endif
	const Range1D
		row_rng = RangePack::full_range(row_rng_in,1,rows_),
		col_rng = RangePack::full_range(col_rng_in,1,rows_),
		row_rng_full(1,rows_),
		col_rng_full(1,cols_);
	const index_type
		*row_i    = row_i_,
		*col_j    = col_j_;
	index_type
		cnt_nz    = 0,
		k         = 0;
	if( dl == -row_rng.ubound() + col_rng.lbound() && du == +col_rng.ubound() - row_rng.lbound() ) {
		// The diagonals are not limiting so we can ignore them
		if( row_rng == row_rng_full && col_rng == col_rng_full ) {
			// The row and column ranges are not limiting either
			cnt_nz = nz_; // Just return the count of all the elements!
		}
		else {
			// The row or column range is limiting
			if( inv_row_perm == NULL && inv_col_perm == NULL ) {
				// Neither the rows nor columns are permuted
				for( k = 0; k < nz_; ++row_i, ++col_j, ++k ) {
					const index_type
						i = (*row_i),
						j = (*col_j);
					VALIDATE_ROW_COL_IN_RANGE();
					cnt_nz += row_rng.in_range(i) && col_rng.in_range(j) ? 1 : 0;
				}
			}
			else if ( inv_row_perm != NULL && inv_col_perm == NULL ) {
				// Only the rows are permuted 
				for( k = 0; k < nz_; ++row_i, ++col_j, ++k ) {
					const index_type
						i = inv_row_perm[(*row_i)-1],
						j = (*col_j);
					VALIDATE_ROW_COL_IN_RANGE();
					cnt_nz += row_rng.in_range(i) && col_rng.in_range(j) ? 1 : 0;
				}
			}
			else if ( inv_row_perm == NULL && inv_col_perm != NULL ) {
				// Only the columns are permuted 
				for( k = 0; k < nz_; ++row_i, ++col_j, ++k ) {
					const index_type
						i = (*row_i),
						j = inv_col_perm[(*col_j)-1];
					VALIDATE_ROW_COL_IN_RANGE();
					cnt_nz += row_rng.in_range(i) && col_rng.in_range(j) ? 1 : 0;
				}
			}
			else {
				// Both the rows and columns are permuted!
				for( k = 0; k < nz_; ++row_i, ++col_j, ++k ) {
					const index_type
						i = inv_row_perm[(*row_i)-1],
						j = inv_col_perm[(*col_j)-1];
					VALIDATE_ROW_COL_IN_RANGE();
					cnt_nz += row_rng.in_range(i) && col_rng.in_range(j) ? 1 : 0;
				}
			}
		}
	}
	else {
		// We have to consider the diagonals dl and du
		assert(0); // ToDo: Implement!
	}
	return cnt_nz;
}

void MatrixSparseCOORSerial::coor_extract_nonzeros(
	EElementUniqueness    element_uniqueness
	,const index_type     inv_row_perm[]
	,const index_type     inv_col_perm[]
	,const Range1D        &row_rng_in
	,const Range1D        &col_rng_in
	,index_type           dl
	,index_type           du
	,value_type           alpha
	,const index_type     len_Aval
	,value_type           Aval[]
	,const index_type     len_Aij
	,index_type           Arow[]
	,index_type           Acol[]
	,const index_type     row_offset
	,const index_type     col_offset
	) const
{
#ifdef _DEBUG
	const char err_msg_head[] = "MatrixSparseCOORSerial::count_nonzeros(...): Error";
	THROW_EXCEPTION(
		element_uniqueness_ == ELEMENTS_ASSUME_DUPLICATES_SUM && element_uniqueness == ELEMENTS_FORCE_UNIQUE
		,std::logic_error
		,err_msg_head << ", the client requests extraction of unique "
		"elements but this sparse matrix object can not guarantee this!" );
#endif
	const Range1D
		row_rng = RangePack::full_range(row_rng_in,1,rows_),
		col_rng = RangePack::full_range(col_rng_in,1,rows_),
		row_rng_full(1,rows_),
		col_rng_full(1,cols_);
	value_type
		*val      = val_;
	const index_type
		*row_i    = row_i_,
		*col_j    = col_j_;
	index_type
		cnt_nz    = 0,
		k         = 0;
	if( dl == -row_rng.ubound() + col_rng.lbound() && du == +col_rng.ubound() - row_rng.lbound() ) {
		// The diagonals are not limiting so we can ignore them
		if( row_rng == row_rng_full && col_rng == col_rng_full ) {
			// The row and column ranges are not limiting either
			if( inv_row_perm == NULL && inv_col_perm == NULL ) {
				// We are just extracting the whole, unpermuted matrix
				for( k = 0; k < nz_; ++val, ++row_i, ++col_j, ++k ) {
					++cnt_nz;
					if( len_Aval )
						*Aval++ = *val;           // ToDo: Split into different loops (no inner if())
					if( len_Aij ) {
						*Arow++ = *row_i + row_offset;
						*Acol++ = *col_j + col_offset;
					}
				}
			}
			else {
				assert(0); // ToDo: Implement!
			}
		}
		else {
			// The row or column range is limiting
			if( inv_row_perm == NULL && inv_col_perm == NULL ) {
				// There are no permutations to consider
				for( k = 0; k < nz_; ++val, ++row_i, ++col_j, ++k ) {
					const index_type
						i = (*row_i),
						j = (*col_j);
					VALIDATE_ROW_COL_IN_RANGE();
					if( row_rng.in_range(i) && col_rng.in_range(j) ) {
						++cnt_nz;
						if( len_Aval )
							*Aval++ = *val;           // ToDo: Split into different loops (no inner if())
						if( len_Aij ) {
							*Arow++ = i + row_offset;
							*Acol++ = j + col_offset;
						}
					}
				}
			}
			else if( inv_row_perm != NULL && inv_col_perm == NULL ) {
				// We must consider only row permutations
				for( k = 0; k < nz_; ++val, ++row_i, ++col_j, ++k ) {
					const index_type
						i = inv_row_perm[(*row_i)-1],
						j = (*col_j);
					VALIDATE_ROW_COL_IN_RANGE();
					if( row_rng.in_range(i) && col_rng.in_range(j) ) {
						++cnt_nz;
						if( len_Aval )
							*Aval++ = *val;           // ToDo: Split into different loops (no inner if())
						if( len_Aij ) {
							*Arow++ = i + row_offset;
							*Acol++ = j + col_offset;
						}
					}
				}
			}
			else if( inv_row_perm == NULL && inv_col_perm != NULL ) {
				// We must consider only column permutations
				for( k = 0; k < nz_; ++val, ++row_i, ++col_j, ++k ) {
					const index_type
						i = (*row_i),
						j = inv_col_perm[(*col_j)-1];
					VALIDATE_ROW_COL_IN_RANGE();
					if( row_rng.in_range(i) && col_rng.in_range(j) ) {
						++cnt_nz;
						if( len_Aval )
							*Aval++ = *val;           // ToDo: Split into different loops (no inner if())
						if( len_Aij ) {
							*Arow++ = i + row_offset;
							*Acol++ = j + col_offset;
						}
					}
				}
			}
			else {
				// We must consider row and column permutations
				for( k = 0; k < nz_; ++val, ++row_i, ++col_j, ++k ) {
					const index_type
						i = inv_row_perm[(*row_i)-1],
						j = inv_col_perm[(*col_j)-1];
					VALIDATE_ROW_COL_IN_RANGE();
					if( row_rng.in_range(i) && col_rng.in_range(j) ) {
						++cnt_nz;
						if( len_Aval )
							*Aval++ = *val;           // ToDo: Split into different loops (no inner if())
						if( len_Aij ) {
							*Arow++ = i + row_offset;
							*Acol++ = j + col_offset;
						}
					}
				}
			}
		}
	}
	else {
		// We have to consider the diagonals dl and du
		assert(0); // ToDo: Implement!
	}
	assert( len_Aval == 0 || len_Aval == cnt_nz );
	assert( len_Aij == 0  || len_Aij  == cnt_nz );
}

// private

void MatrixSparseCOORSerial::make_storage_unique()
{
	if( release_resource_.count() > 1 ) {
		assert(0); // ToDo: Allocate new storage and copy this memory.
		self_allocate_ = true;
	}
}

} // end namespace SparseLinAlgPack
