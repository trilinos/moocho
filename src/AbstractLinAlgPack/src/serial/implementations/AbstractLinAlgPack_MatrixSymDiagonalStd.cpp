// ////////////////////////////////////////////////////////////////
// MatrixSymDiagonalStd.cpp

#include "SparseLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/GenMatrixOp.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"

namespace SparseLinAlgPack {

MatrixSymDiagonalStd::MatrixSymDiagonalStd()
{}

VectorSlice MatrixSymDiagonalStd::diag()
{
	return diag_();
}

const VectorSlice MatrixSymDiagonalStd::diag() const
{
	return diag_();
}

// Overridden from MatrixSymInitDiagonal

void MatrixSymDiagonalStd::init_identity( size_type n, value_type alpha = 1.0 )
{
	diag_.resize(n);
	if(n)
		diag_ = alpha;
}

void MatrixSymDiagonalStd::init_diagonal( const VectorSlice& diag )
{
	diag_ = diag;
}

// Overridden from Matrix

size_type MatrixSymDiagonalStd::rows() const
{
	return diag_.size();
}

// Overridden from MatrixWithOp

void MatrixSymDiagonalStd::Mp_StM(
	GenMatrixSlice* gms_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const
{
	using LinAlgPack::Vp_StV;
	// Assert that the dimensions match up.
	LinAlgPack::Mp_M_assert_sizes( gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
		, rows(), cols(), trans_rhs );
	// Add to the diagonal
	Vp_StV( &gms_lhs->diag(), alpha, diag_ );
}

void MatrixSymDiagonalStd::Vp_StMtV(
	VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2, value_type beta) const
{
	const VectorSlice diag = this->diag();
	size_type n = diag.size();

	// y = b*y + a * op(A) * x
	//
	// A is symmetric and diagonal A = diag(diag) so:
	//
	// y(j) = b*y(j) + a * diag(j) * x(j), for j = 1...n
	
	LinAlgPack::Vp_MtV_assert_sizes(
		vs_lhs->size(), n, n, trans_rhs1, vs_rhs2.size() );

	if( vs_rhs2.stride() == 1 && vs_lhs->stride() == 1 ) {
		// Optimized implementation
		const value_type
			*d_itr      = diag.raw_ptr(),
			*x_itr      = vs_rhs2.raw_ptr();
		value_type
			*y_itr      = vs_lhs->raw_ptr(),
			*y_end      = y_itr + vs_lhs->size();

		if( beta == 0.0 ) {
			while( y_itr != y_end )
				*y_itr++ = alpha * (*d_itr++) * (*x_itr++);
		}
		else if( beta == 1.0 ) {
			while( y_itr != y_end )
				*y_itr++ += alpha * (*d_itr++) * (*x_itr++);
		}
		else {
			for( ; y_itr != y_end; ++y_itr )
				*y_itr = beta * (*y_itr) + alpha * (*d_itr++) * (*x_itr++);
		}
	}
	else {
		// Generic implementation
		VectorSlice::const_iterator
			d_itr = diag.begin(),
			x_itr = vs_rhs2.begin();
		VectorSlice::iterator
			y_itr = vs_lhs->begin(),
			y_end = vs_lhs->end();
		for( ; y_itr != y_end; ++y_itr, ++d_itr, ++x_itr ) {
#ifdef LINALGPACK_CHECK_RANGE
			assert( d_itr < diag.end() );
			assert( x_itr < vs_rhs2.end() );
			assert( y_itr < vs_lhs->end() );
#endif
			*y_itr = beta * (*y_itr) + alpha * (*d_itr) * (*x_itr);
		}
	}
}

void MatrixSymDiagonalStd::Vp_StMtV(
	VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2, value_type beta) const
{
	const VectorSlice diag = this->diag();
	size_type n = diag.size();

	// y = b*y + a * op(A) * x
	//
	// A is symmetric and diagonal A = diag(diag) so:
	//
	// y(j) = b*y(j) + a * diag(j) * x(j), for j = 1...n
	//
	// x is sparse so take account of this.
	
	LinAlgPack::Vp_MtV_assert_sizes( vs_lhs->size()
		, n, n, trans_rhs1, sv_rhs2.size() );

	for(   SpVectorSlice::const_iterator x_itr = sv_rhs2.begin()
		 ; x_itr != sv_rhs2.end()
		 ; ++x_itr )
	{
		(*vs_lhs)(x_itr->indice() + sv_rhs2.offset())
			= beta * (*vs_lhs)(x_itr->indice() + sv_rhs2.offset())
				+ alpha * diag(x_itr->indice() + sv_rhs2.offset()) * x_itr->value();
			// Note: The indice x(i) invocations are ranged check
			// if this is compiled into the code.
	}
}

// Overridden from MatrixWithOpFactorized

void MatrixSymDiagonalStd::V_InvMtV(
	VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
	, const VectorSlice& vs_rhs2) const
{
	const VectorSlice diag = this->diag();
	size_type n = diag.size();

	// y = inv(op(A)) * x
	//
	// A is symmetric and diagonal (A = diag(diag)) so:
	//
	// y(j) = x(j) / diag(j), for j = 1...n

	LinAlgPack::Vp_MtV_assert_sizes( vs_lhs->size()
		, n, n, trans_rhs1, vs_rhs2.size() );
	
	if( vs_rhs2.stride() == 1 && vs_lhs->stride() == 1 ) {
		// Optimized implementation
		const value_type
			*d_itr      = diag.raw_ptr(),
			*x_itr      = vs_rhs2.raw_ptr();
		value_type
			*y_itr      = vs_lhs->raw_ptr(),
			*y_end      = y_itr + vs_lhs->size();
		while( y_itr != y_end )
			*y_itr++ = (*x_itr++) / (*d_itr++);
	}
	else {
		// Generic implementation
		VectorSlice::const_iterator
			d_itr = diag.begin(),
			x_itr = vs_rhs2.begin();
		VectorSlice::iterator
			y_itr = vs_lhs->begin(),
			y_end = vs_lhs->end();
		for( ; y_itr != y_end; ++y_itr, ++d_itr, ++x_itr ) {
			assert( d_itr < diag.end() );
			assert( x_itr < vs_rhs2.end() );
			assert( y_itr < vs_lhs->end() );
			*y_itr = (*x_itr)/(*d_itr);
		}
	}
}

void MatrixSymDiagonalStd::V_InvMtV(
	VectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
	, const SpVectorSlice& sv_rhs2) const
{
	const VectorSlice diag = this->diag();
	size_type n = diag.size();

	// y = inv(op(A)) * x
	//
	// A is symmetric and diagonal A = diag(diag) so:
	//
	// y(j) = x(j) / diag(j), for j = 1...n
	//
	// x is sparse so take account of this.
	
	LinAlgPack::Vp_MtV_assert_sizes( vs_lhs->size()
		, n, n, trans_rhs1, sv_rhs2.size() );

	for(   SpVectorSlice::const_iterator x_itr = sv_rhs2.begin()
		 ; x_itr != sv_rhs2.end()
		 ; ++x_itr )
	{
		(*vs_lhs)(x_itr->indice() + sv_rhs2.offset())
			= x_itr->value() / diag(x_itr->indice() + sv_rhs2.offset());
			// Note: The indice x(i) invocations are ranged check
			// if this is compiled into the code.
	}
}

} // end namespace SparseLinAlgPack