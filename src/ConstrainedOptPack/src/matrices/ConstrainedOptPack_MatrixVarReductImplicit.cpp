// //////////////////////////////////////////////////////////////////////////////////
// MatrixVarReductImplicit.cpp
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

#include "ConstrainedOptPack_MatrixVarReductImplicit.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSlice.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_AssertOp.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TestForException.hpp"

namespace {

//
// Implicit matrix-vector multiplication:
//
// y = b*y + a*op(inv(C)*N)*x
//
template<class V>
void imp_Vp_StMtV_implicit(
	AbstractLinAlgPack::VectorMutable               *y
	,AbstractLinAlgPack::value_type                       a
	,const AbstractLinAlgPack::MatrixOpNonsing    &C
	,const AbstractLinAlgPack::MatrixOp               &N
	,BLAS_Cpp::Transp                                     D_trans
	,const V                                              &x
	,DenseLinAlgPack::value_type                               b
	)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	namespace alap = AbstractLinAlgPack;
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

	const DenseLinAlgPack::size_type
		r   = C.rows(),
		dof = N.cols();

	// ToDo: Pass in workspace vectors to save some allocation time!

	if( D_trans == no_trans ) {
		// y = b*y
		alap::Vt_S(y,b);
		//
		// y += -a * inv(C) * ( N * x )
		//
		alap::VectorSpace::vec_mut_ptr_t
			t1 = N.space_cols().create_member(),
			t2 = C.space_rows().create_member();
		// t1 = N*x
		LinAlgOpPack::V_MtV( t1.get(), N, no_trans, x );
		// t2 = inv(C) * t1
		alap::V_InvMtV( t2.get(), C, no_trans, *t1 );
		// y += a*t2
		alap::Vp_StV( y, -a, *t2 );
	}
	else {
		//
		// y = b*y - a * N' * ( inv(C') * x )
		//
		alap::VectorSpace::vec_mut_ptr_t
			t1 = C.space_cols().create_member();
		// t1 = inv(C')*x
		alap::V_InvMtV( t1.get(), C, trans, x );
		// y = b*y - a*N'*t1
	    alap::Vp_StMtV( y, -a,  N, trans, *t1, b );
	}
}

/*

//
// Generate a row of inv(C)*N if not already computed.
//
void imp_compute_InvCtN_row(
	DenseLinAlgPack::size_type                                                 r
	,const ConstrainedOptPack::DecompositionSystemVarReduct      &decomp_sys
	,DenseLinAlgPack::size_type                                                j
	,DenseLinAlgPack::DVectorSlice                                              *e_j  // Set to all zeros on input and output!
	,ConstrainedOptPack::MatrixVarReductImplicit::InvCtN_rows_t  *InvCtN_rows
	)
{
	typedef  DenseLinAlgPack::value_type  value_type;
	using DenseLinAlgPack::DVectorSlice;
	if( (*InvCtN_rows)[j-1] == NULL ) {
		// Generate row j of inv(C)*N
		value_type *vec = (*InvCtN_rows)[j-1] = new value_type[r]; // ToDo: We may want to allocate more vectors at once!
		DVectorSlice row_j(vec,r);
		// row_j = N'*inv(C')*e_j
		(*e_j)(j) = 1.0;
		imp_Vp_StMtV_implicit( &row_j, 1.0, decomp_sys, BLAS_Cpp::trans, *e_j, 0.0 );
		(*e_j)(j) = 0.0;
	}
}

//
// Perform the matrix-vector multiplication:
// 
// y = b*y -a * op(P) * [inv(C) * N] * x
//
// by generating rows [inv(C)*N](j,:) for each nonzero entry op(P)(i,j).
//
template<class V>
void imp_Vp_StPtMtV_by_row(
	DenseLinAlgPack::DVectorSlice                                               *y
	,DenseLinAlgPack::value_type                                               a
	,const ConstrainedOptPack::GenPermMatrixSlice                &P
	,BLAS_Cpp::Transp                                                     P_trans
	,const ConstrainedOptPack::DecompositionSystemVarReduct      &decomp_sys
	,const V                                                              &x
	,DenseLinAlgPack::value_type                                               b
	,ConstrainedOptPack::MatrixVarReductImplicit::InvCtN_rows_t *InvCtN_rows
	)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::trans_not;
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	using DenseLinAlgPack::dot;
	using DenseLinAlgPack::DVectorSlice;
	using AbstractLinAlgPack::dot;
	using AbstractLinAlgPack::Vp_StMtV;
	using AbstractLinAlgPack::GenPermMatrixSlice;
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
	
	const DenseLinAlgPack::size_type
		D_rows = decomp_sys.C().rows(),
		D_cols = decomp_sys.N().cols();
	// y = b*y
	if(b==0.0)       *y = 0.0;
	else if(b!=1.0)  DenseLinAlgPack::Vt_S(y,b);
	// Compute t = N'*inv(C')*e(j) then y(i) += -a*t'*x where op(P)(i,j) = 1.0
	Workspace<DenseLinAlgPack::value_type>   e_j_ws(wss,D_rows);
	DVectorSlice                              e_j(&e_j_ws[0],e_j_ws.size());
	e_j = 0.0;
	for( GenPermMatrixSlice::const_iterator itr = P.begin(); itr != P.end(); ++itr ) {
		const DenseLinAlgPack::size_type
			i = P_trans == no_trans ? itr->row_i() : itr->col_j(),
			j = P_trans == no_trans ? itr->col_j() : itr->row_i();
		// t = op(M') * e(j)
		imp_compute_InvCtN_row(D_rows,decomp_sys,j,&e_j,InvCtN_rows);
		DVectorSlice t((*InvCtN_rows)[j-1],D_cols);
		// y(i) += -a*t'*x
		(*y)(i) += (-a) * dot(t,x);
	}
}

*/

} // end namespace

namespace ConstrainedOptPack {

void MatrixVarReductImplicit::initialize(
	const mat_nonsing_ptr_t          &C
	,const mat_ptr_t                 &N
	,const mat_ptr_t                 &D_direct
	)
{
	namespace rcp = MemMngPack;
	// Validate the inputs
	TEST_FOR_EXCEPTION(
		C.get() == NULL, std::invalid_argument
		,"MatrixVarReductImplicit::initialize(...): Error, "
		"C.get() must not be NULL" );
	TEST_FOR_EXCEPTION(
		N.get() == NULL, std::invalid_argument
		,"MatrixVarReductImplicit::initialize(...): Error, "
		"N.get() must not be NULL" );
	if( D_direct.get() ) {
		const bool is_compatible_cols = D_direct->space_cols().is_compatible(C->space_cols());
		TEST_FOR_EXCEPTION(
			!is_compatible_cols, VectorSpace::IncompatibleVectorSpaces
			,"MatrixVarReductImplicit::initialize(...): Error, "
			"D_direct->space_cols() is not compatible with C->space_cols()" );
		const bool is_compatible_rows = D_direct->space_rows().is_compatible(N->space_rows());
		TEST_FOR_EXCEPTION(
			!is_compatible_rows, VectorSpace::IncompatibleVectorSpaces
			,"MatrixVarReductImplicit::initialize(...): Error, "
			"D_direct->space_rows() is not compatible with N->space_rows()" );
	}
	// Set the members
	C_        = C;
	N_        = N;
	D_direct_ = D_direct;
	if(!InvCtN_rows_set_list_.empty()) { // Free previously allocated vectors
		for( InvCtN_rows_set_list_t::iterator itr = InvCtN_rows_set_list_.begin();
			 itr != InvCtN_rows_set_list_.end(); ++itr )
        {
			InvCtN_rows_[*itr] = Teuchos::null;
		}
		InvCtN_rows_set_list_.clear();
	}
}

void MatrixVarReductImplicit::set_uninitialized()
{
	namespace rcp = MemMngPack;
	C_        = Teuchos::null;
	N_        = Teuchos::null;
	D_direct_ = Teuchos::null;
}

// Overridden from MatrixBase

size_type MatrixVarReductImplicit::rows() const
{
	return C_.get() ? C_->rows() : 0;
}

size_type MatrixVarReductImplicit::cols() const
{
	return N_.get() ? N_->cols() : 0;
}

// Overridden from MatrixOp

const VectorSpace& MatrixVarReductImplicit::space_cols() const
{
	assert_initialized();
	return C_->space_cols();
}

const VectorSpace& MatrixVarReductImplicit::space_rows() const
{
	assert_initialized();
	return N_->space_rows();
}

MatrixOp& MatrixVarReductImplicit::operator=(const MatrixOp& M)
{
	assert_initialized();
	assert(0); // ToDo: Finish!
	return *this;
}

std::ostream& MatrixVarReductImplicit::output(std::ostream& o) const
{
	o << "This is a " << this->rows() << " x " << this->cols()
	  << " variable reduction matrix D = -inv(C)*N where C and N are:\n"
	  << "C =\n" << *C_
	  << "N =\n" << *N_;
	return o;
}

void MatrixVarReductImplicit::Vp_StMtV(
	VectorMutable* y, value_type a
	,BLAS_Cpp::Transp D_trans
	,const Vector& x, value_type b
	) const
{
	assert_initialized();
	AbstractLinAlgPack::Vp_MtV_assert_compatibility(y,*this,D_trans,x);
	imp_Vp_StMtV_implicit( y, a, *C_, *N_, D_trans, x, b );
}

void MatrixVarReductImplicit::Vp_StMtV(
	VectorMutable* y, value_type a
	,BLAS_Cpp::Transp D_trans
	,const SpVectorSlice& x, value_type b
	) const
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

	assert_initialized();
	AbstractLinAlgPack::Vp_MtV_assert_compatibility(y,*this,D_trans,x);
	imp_Vp_StMtV_implicit( y, a, *C_, *N_, D_trans, x, b );
/*
	const size_type
		D_rows = this->rows(), D_cols = this->cols(),
		opD_rows = rows(D_rows,D_cols,D_trans), opD_cols = cols(D_rows,D_cols,D_trans);
	DenseLinAlgPack::Vp_MtV_assert_sizes(y->size(),D_rows,D_cols,D_trans,x.size());
	if( use_dense_mat_vec_ && D_dense_.rows() > 0 ) {
		LinAlgOpPack::Vp_StMtV( y, a, D_dense_, D_trans, x, b );
	}
	else {
		if( x.nz() == x.size() ) {  // This is B.S.  Should use MatrixWithOpFactorized object for C!
			DVectorSlice dx = AbstractLinAlgPack::dense_view(x);
			imp_Vp_StMtV_implicit( y, -a, *decomp_sys_, D_trans, dx, b );
		}
		else if( D_trans == BLAS_Cpp::trans && x.nz() < D_cols ) {
			//
			// y = b*y + (-a)*[N'*inv(C')]*x
			//
			// We can do something crafty here.  We can generate columns of N'*inv(C')
			// and then perform y += -a*[N'*inv(C')](:,j)*x(j) for nonzero x(j)
			//
			Workspace<DenseLinAlgPack::value_type>   e_j_ws(wss,D_rows);
			DVectorSlice                              e_j(&e_j_ws[0],e_j_ws.size());
			e_j = 0.0;
			// y = b*y
			if(b==0.0)       *y = 0.0;
			else if(b!=1.0)  Vt_S(y,b);
			// y += -a*[N'*inv(C')](:,j)*x(j), for j <: { j | x(j) != 0.0 }
			const SpVectorSlice::difference_type o = x.offset();
			for( SpVectorSlice::const_iterator itr = x.begin(); itr != x.end(); ++itr ) {
				const size_type j = itr->indice() + o;
				imp_compute_InvCtN_row(D_rows,*decomp_sys_,j,&e_j,&InvCtN_rows_);
				DenseLinAlgPack::Vp_StV( y, -a * itr->value(), DVectorSlice(InvCtN_rows_[j-1],D_cols) );
			}
		}
		else {   // This is B.S.  Should use MatrixWithOpFactorized object for C!
			DVector dx;
			LinAlgOpPack::assign( &dx, x );
			imp_Vp_StMtV_implicit( y, -a, *decomp_sys_, D_trans, dx(), b );
		}
	}
*/
	// ToDo: Consider implementing the above specialized implementation!
}

void MatrixVarReductImplicit::Vp_StPtMtV(
	VectorMutable* y, value_type a
	,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	,BLAS_Cpp::Transp D_trans
	,const Vector& x, value_type b
	) const
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	assert_initialized();
/*
	const size_type
		D_rows = this->rows(), D_cols = this->cols(),
		opD_rows = rows(D_rows,D_cols,D_trans), opD_cols = cols(D_rows,D_cols,D_trans);
	DenseLinAlgPack::Vp_MtV_assert_sizes(y->size(),P.rows(),P.cols(),P_trans,opD_rows);
	DenseLinAlgPack::Vp_MtV_assert_sizes(cols(P.rows(),P.cols(),P_trans),D_rows,D_cols,D_trans,x.size());
	if( D_dense_.rows() > 0 ) {
		AbstractLinAlgPack::dense_Vp_StPtMtV(y,a,P,P_trans,D_dense_,D_trans,x,b);
	}
	else if( P.nz() > D_cols || D_trans == BLAS_Cpp::trans ) {
		// Just use the default implementation
		MatrixOp::Vp_StPtMtV(y,a,P,P_trans,D_trans,x,b);
	}
	else {
		imp_Vp_StPtMtV_by_row(y,a,P,P_trans,*decomp_sys_,x,b,&InvCtN_rows_);
	}
*/
	MatrixOp::Vp_StPtMtV(y,a,P,P_trans,D_trans,x,b); // ToDo:Update specialized implementation above!
}

void MatrixVarReductImplicit::Vp_StPtMtV(
	VectorMutable* y, value_type a
	,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	,BLAS_Cpp::Transp D_trans
	,const SpVectorSlice& x, value_type b
	) const
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	assert_initialized();
/*
	const size_type
		D_rows = this->rows(), D_cols = this->cols(),
		opD_rows = rows(D_rows,D_cols,D_trans), opD_cols = cols(D_rows,D_cols,D_trans);
	DenseLinAlgPack::Vp_MtV_assert_sizes(y->size(),P.rows(),P.cols(),P_trans,opD_rows);
	DenseLinAlgPack::Vp_MtV_assert_sizes(cols(P.rows(),P.cols(),P_trans),D_rows,D_cols,D_trans,x.size());
	if( D_dense_.rows() > 0 ) {
		AbstractLinAlgPack::dense_Vp_StPtMtV(y,a,P,P_trans,D_dense_,D_trans,x,b);
	}
	else if( P.nz() > D_cols || D_trans == BLAS_Cpp::trans ) {
		// Just use the default implementation
		MatrixOp::Vp_StPtMtV(y,a,P,P_trans,D_trans,x,b);
	}
	else {
		imp_Vp_StPtMtV_by_row(y,a,P,P_trans,*decomp_sys_,x,b,&InvCtN_rows_);
	}
*/
	MatrixOp::Vp_StPtMtV(y,a,P,P_trans,D_trans,x,b); // ToDo:Update specialized implementation above!
}

// Private member functions

void MatrixVarReductImplicit::assert_initialized() const
{
	TEST_FOR_EXCEPTION(
		C_.get() == NULL, std::logic_error
		,"MatrixVarReductImplicit::assert_initialized(): Error, "
		"initialize(...) has not been called yet!" );
}

}	// end namespace ConstrainedOptPack 
