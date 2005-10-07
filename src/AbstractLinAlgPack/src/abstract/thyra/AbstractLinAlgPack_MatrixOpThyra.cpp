// ///////////////////////////////////////////////////////////
// MatrixOpThyra.cpp
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

#include "MatrixOpThyra.hpp"
#include "VectorMutableThyra.hpp"
#include "ThyraAccessors.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace AbstractLinAlgPack {

// Constructors / Initializers

MatrixOpThyra::MatrixOpThyra()
{}

MatrixOpThyra::MatrixOpThyra(
	const Teuchos::RefCountPtr<const Thyra::LinearOpBase<value_type> >   &thyra_linear_op
	,BLAS_Cpp::Transp                                                    thyra_linear_op_trans
	)
{
	this->initialize(thyra_linear_op,thyra_linear_op_trans);
}

void MatrixOpThyra::initialize(
	const Teuchos::RefCountPtr<const Thyra::LinearOpBase<value_type> >   &thyra_linear_op
	,BLAS_Cpp::Transp                                                    thyra_linear_op_trans
	)
{
	namespace mmp = MemMngPack;
	TEST_FOR_EXCEPTION(
		thyra_linear_op.get()==NULL, std::invalid_argument
		,"MatrixOpThyra::initialize(thyra_linear_op): Error!"
		);
	const bool adjointSupported = ( ::Thyra::opSupported(*thyra_linear_op,Thyra::NOTRANS) && ::Thyra::opSupported(*thyra_linear_op,Thyra::TRANS) );
	TEST_FOR_EXCEPTION(
		!adjointSupported, std::invalid_argument
		,"MatrixOpThyra::initialize(thyra_linear_op): Error, the operator opSupported(thyra_linear_op,transp) must return true "
		"for both values of transp==NOTRANS and transp=TRANS!"
		);
	thyra_linear_op_       = thyra_linear_op;
	thyra_linear_op_trans_ = thyra_linear_op_trans;
	space_cols_.initialize(thyra_linear_op_->range());
	space_rows_.initialize(thyra_linear_op_->domain());
}

Teuchos::RefCountPtr<const Thyra::LinearOpBase<value_type> > 
MatrixOpThyra::set_uninitialized()
{
	Teuchos::RefCountPtr<const Thyra::LinearOpBase<value_type> > tmp_thyra_linear_op = thyra_linear_op_;
	thyra_linear_op_       = Teuchos::null;
	thyra_linear_op_trans_ = BLAS_Cpp::no_trans;
	space_cols_.set_uninitialized();
	space_rows_.set_uninitialized();
	return tmp_thyra_linear_op;
}

// Overridden from MatrixBase

const VectorSpace&
MatrixOpThyra::space_cols() const
{
	return space_cols_;
}

const VectorSpace&
MatrixOpThyra::space_rows() const
{
	return space_rows_;
}

// Overridden from MatrixOp

MatrixOp::mat_mut_ptr_t
MatrixOpThyra::clone()
{
	return Teuchos::rcp(new MatrixOpThyra(thyra_linear_op_->clone()));
}

MatrixOp& MatrixOpThyra::operator=(const MatrixOp& mwo_rhs)
{
	using Teuchos::dyn_cast;
	const MatrixOpThyra &mwo_rhs_thyra = dyn_cast<const MatrixOpThyra>(mwo_rhs);
	thyra_linear_op_       = mwo_rhs_thyra.thyra_linear_op_;  // ToDo: Clone this!
	thyra_linear_op_trans_ = mwo_rhs_thyra.thyra_linear_op_trans_;
	space_cols_           = mwo_rhs_thyra.space_cols_;  // ToDo: Clone this!
	space_rows_           = mwo_rhs_thyra.space_rows_;  // ToDo: Clone this!
}

void MatrixOpThyra::Vp_StMtV(
	VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
	,const Vector& v_rhs2, value_type beta
	) const
{
	using BLAS_Cpp::trans_trans;
	// Get Thyra views of the vectors
	Teuchos::RefCountPtr<const Thyra::VectorBase<value_type> > thyra_vec_rhs2;
	get_thyra_vector( BLAS_Cpp::no_trans==trans_rhs1 ? space_rows_ : space_cols_, v_rhs2, &thyra_vec_rhs2 );
	Teuchos::RefCountPtr<Thyra::VectorBase<value_type> > thyra_vec_lhs;
	get_thyra_vector( BLAS_Cpp::no_trans==trans_rhs1 ? space_cols_ : space_rows_, v_lhs, &thyra_vec_lhs );
	// Perform the multiplication
  ::Thyra::apply(
    *thyra_linear_op_
		,trans_trans(trans_rhs1,thyra_linear_op_trans())==BLAS_Cpp::no_trans ? Thyra::NOTRANS : Thyra::TRANS  // M_trans
		,*thyra_vec_rhs2                                                                                      // x
		,&*thyra_vec_lhs                                                                                      // y
		,alpha                                                                                                // alpha
		,beta                                                                                                // beta
		);
	// Free/commit Thyra vector views
	free_thyra_vector( BLAS_Cpp::no_trans==trans_rhs1 ? space_rows_ : space_cols_, v_rhs2, &thyra_vec_rhs2 );
	commit_thyra_vector( BLAS_Cpp::no_trans==trans_rhs1 ? space_cols_ : space_rows_, v_lhs, &thyra_vec_lhs );
}

bool MatrixOpThyra::Mp_StMtM(
	MatrixOp* mwo_lhs, value_type alpha
	,BLAS_Cpp::Transp trans_rhs1
	,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	,value_type beta
	) const
{
	return false; // ToDo: Specialize!
}

} // end namespace AbstractLinAlgPack