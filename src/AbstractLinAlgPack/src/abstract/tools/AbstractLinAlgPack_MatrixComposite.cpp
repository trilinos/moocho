// /////////////////////////////////////////////////////////////////////
// MatrixCompositeStd.cpp
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

#include "AbstractLinAlgPack/include/MatrixCompositeStd.h"
#include "AbstractLinAlgPack/include/SpVectorClass.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutableCompositeStd.h"
#include "AbstractLinAlgPack/include/LinAlgPackAssertOp.h"
//#include "AbstractLinAlgPack/include/GenPermMatrixSliceOp.h"
#include "WorkspacePack.h"
#include "Range1D.h"
#include "ThrowException.h"
#include "profile_hack.h"

namespace {

// Get an element from a vector (two versions)

inline
AbstractLinAlgPack::value_type
get_element( const AbstractLinAlgPack::VectorWithOp& v, AbstractLinAlgPack::index_type i )
{
 	return v.get_ele(i);
}

inline
AbstractLinAlgPack::value_type
get_element( const AbstractLinAlgPack::SpVectorSlice& v, AbstractLinAlgPack::index_type i )
{
	const AbstractLinAlgPack::SpVectorSlice::element_type
		*ele = v.lookup_element(i);
	return ele != NULL ? ele->value() : 0.0;
}

// Get a view of a vector (two versions)

inline
ReferenceCountingPack::ref_count_ptr<const AbstractLinAlgPack::VectorWithOp>
get_view(
	const AbstractLinAlgPack::VectorWithOp& v
	,AbstractLinAlgPack::index_type l
	,AbstractLinAlgPack::index_type u
	)
{
 	return v.sub_view(l,u);
}

inline
ReferenceCountingPack::ref_count_ptr<const AbstractLinAlgPack::SpVectorSlice>
get_view(
	const AbstractLinAlgPack::SpVectorSlice& v
	,AbstractLinAlgPack::index_type l
	,AbstractLinAlgPack::index_type u
	)
{
	return ReferenceCountingPack::ref_count_ptr<const AbstractLinAlgPack::SpVectorSlice>(
		new AbstractLinAlgPack::SpVectorSlice( v(l,u) ) );
}

// Perform a matrix vector multiplication

template<class V>
void Vp_StMtV_imp(
	AbstractLinAlgPack::VectorWithOpMutable* y, AbstractLinAlgPack::value_type a
	,AbstractLinAlgPack::size_type M_rows, AbstractLinAlgPack::size_type M_cols
	,const AbstractLinAlgPack::MatrixCompositeStd::matrix_list_t& mat_list
	,const AbstractLinAlgPack::MatrixCompositeStd::vector_list_t& vec_list
	,BLAS_Cpp::Transp M_trans
	,const V& x, AbstractLinAlgPack::value_type b
	)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	using AbstractLinAlgPack::dot;  // We should not have to do this but some compiles &%!#$
	typedef AbstractLinAlgPack::MatrixCompositeStd::matrix_list_t  mat_list_t;
	typedef AbstractLinAlgPack::MatrixCompositeStd::vector_list_t  vec_list_t;

	AbstractLinAlgPack::Vp_MtV_assert_sizes( y->dim(), M_rows, M_cols, M_trans, x.dim() );

    AbstractLinAlgPack::Vt_S( y, b );  // Will take care of b == 0.0

	if( vec_list.size() ) {
		for( vec_list_t::const_iterator itr = vec_list.begin(); itr != vec_list.end(); ++itr ) {
			const BLAS_Cpp::Transp
				op_v_trans = ( M_trans == itr->v_trans_ ? no_trans : trans );
			const AbstractLinAlgPack::index_type
				r = ( M_trans == no_trans ? itr->r_l_ : itr->c_l_ ),
				c = ( M_trans == no_trans ? itr->c_l_ : itr->r_l_ );
			if( itr->rng_G_.full_range() ) { // op(v)
				if( op_v_trans == no_trans ) {
					//
					//         [ y1 ]         [              ]   [ x1 ]
					// r:r+n-1 [ y2 ] +=  a * [   beta * v   ] * [ x2 ] c:c 
					//         [ y3 ]         [              ]   [ x3 ]
					//                            \______/        
					//                              c:c
					// =>
					//
					// y(r:r+n-1) += (a * beta * x(c)) * v
					// 
					AbstractLinAlgPack::Vp_StV( y->sub_view(r,r+itr->v_->dim()-1).get(), a * itr->beta_ * get_element(x,c), *itr->v_ );
				}
				else {
					//
					//     [ y1 ]        [                ]   [ x1 ]
					// r:r [ y2 ] += a * [     beta*v'    ] * [ x2 ] c:c+n-1
					//     [ y3 ]        [                ]   [ x3 ]
					//                         \_____/  
					//                         c:c+n-1
					// =>
					//
					// y(r) += a * beta * v'*x(c,c+n-1)
					//
//					y->set_ele( r, y->get_ele(r) + a * itr->beta_ * dot( *itr->v_, *get_view(x,c,c+itr->v_->dim()-1) ) );
					assert(0); // ToDo: Implement the above method in VectorStdOps for VectorWithOp,SpVectorSlice!
				}
			}
			else { // op(op(G)*v) or op(v(rng_G))
				assert(0); // ToDo: Implement when needed!
			}
		}
	}
	if( mat_list.size() ) {
		for( mat_list_t::const_iterator itr = mat_list.begin(); itr != mat_list.end(); ++itr ) {
			const AbstractLinAlgPack::index_type
				rl = rows(itr->r_l_,itr->c_l_,M_trans),
				ru = rows(itr->r_u_,itr->c_u_,M_trans),
				cl = cols(itr->r_l_,itr->c_l_,M_trans),
				cu = cols(itr->r_u_,itr->c_u_,M_trans);
			const BLAS_Cpp::Transp
				op_P_trans = ( M_trans == itr->P_trans_ ? no_trans : trans ),
				op_A_trans = ( M_trans == itr->A_trans_ ? no_trans : trans ),
				op_Q_trans = ( M_trans == itr->Q_trans_ ? no_trans : trans );
			if( itr->rng_P_.full_range() && itr->rng_Q_.full_range() ) { // op(A)
				//
				//       [ y1 ]        [                        ]   [ x1 ]
				// rl:ru [ y2 ] += a * [    alpha * op(op(A))   ] * [ x2 ] cl:cu
				//       [ y3 ]        [                        ]   [ x3 ]
				//                          \_______________/
				//                               cl:cu
				// =>
				//
				// y(rl:ru) += a * alpha * op(op(A)) * x(cl:cu)
				//
				AbstractLinAlgPack::Vp_StMtV( y->sub_view(rl,ru).get(), a * itr->alpha_, *itr->A_, op_A_trans, *get_view(x,cl,cu) );
			}
			else {
				if( itr->A_ == NULL ) { // op(P)
					assert( itr->P_.get() && !itr->P_->is_identity() );
					//
					//       [ y1 ]        [                        ]   [ x1 ]
					// rl:ru [ y2 ] += a * [    alpha * op(op(P))   ] * [ x2 ] cl:cu
					//       [ y3 ]        [                        ]   [ x3 ]
					//                          \_______________/
					//                               cl:cu
					// =>
					//
					// y(rl:ru) += a * alpha * op(op(P)) * x(cl:cu)
					//
// 					AbstractLinAlgPack::Vp_StMtV( y->sub_view(rl,ru).get(), a * itr->alpha_, itr->P_, op_P_trans, *get_view(x,cl,cu) );
					assert(0); // ToDo: Implement the above method properly!
				}
				else { // op(P)*op(A)*op(Q)  [or some simplification]
					//
					//       [ y1 ]        [                                ]   [ x1 ]
					// rl:ru [ y2 ] += a * [    alpha * op(P)*op(A)*op(Q)   ] * [ x2 ] cl:cu
					//       [ y3 ]        [                                ]   [ x3 ]
					//                          \______________________/
					//                                    cl:cu
					// =>
					//
					// y(rl:ru) += a * alpha * op(P)*op(A)*op(Q) * x(cl:cu)
					//
					if( !itr->rng_P_.full_range() && !itr->rng_Q_.full_range() ) { // op(A)(rng_Q,rng_Q)
						AbstractLinAlgPack::Vp_StMtV(
							y->sub_view(rl,ru).get(), a * itr->alpha_
							,*itr->A_->sub_view(
								itr->A_trans_==no_trans  ? itr->rng_P_ : itr->rng_Q_
								,itr->A_trans_==no_trans ? itr->rng_Q_ : itr->rng_P_ )
							,op_A_trans
							,*get_view(x,cl,cu)
							);
					}
					else {
						assert(0); // ToDo: Implement when needed!
					}
				}
			}
		}
	}
}

// Define a comparison object for comparing SubVectorEntry or SubMatrixEntry
// objects for storting by row or by column
template<class T>
class CompSubEntry {
public:
	enum EByRowCol { BY_ROW, BY_COL };
	CompSubEntry(enum EByRowCol by_row_col)
		: by_row_col_(by_row_col)
		{}
	bool operator()( const T& e1, const T& e2 )
		{
			return
				( by_row_col_ == BY_ROW
				  ? e1.r_l_ < e2.r_l_
				  : e1.c_l_ < e2.c_l_
					);
		}
private:
	EByRowCol  by_row_col_;
	CompSubEntry(); // not defined and not to be called
}; // end class CompSubVectorEntry


} // end namespace

namespace AbstractLinAlgPack {

MatrixCompositeStd::MatrixCompositeStd( size_type rows, size_type cols )
{
	reinitalize(rows,cols);
}

void MatrixCompositeStd::reinitalize( size_type rows, size_type cols )
{
	fully_constructed_ = true;
	rows_ = rows;
	cols_ = cols;
	if(matrix_list_.size()) {
		matrix_list_.erase(matrix_list_.begin(),matrix_list_.end());
	}
	if(vector_list_.size()) {
		vector_list_.erase(vector_list_.begin(),vector_list_.end());
	}
	space_rows_  = NULL;
	space_cols_  = NULL;
}

void MatrixCompositeStd::add_vector(
	size_type                      row_offset
	,size_type                     col_offset
	,value_type                    beta
	,const GenPermMatrixSlice      *G
	,const release_resource_ptr_t  &G_release
	,BLAS_Cpp::Transp              G_trans
	,const VectorWithOp            *v
	,const release_resource_ptr_t  &v_release
	,BLAS_Cpp::Transp              v_trans
	)
{
	fully_constructed_ = false;
	assert(0); // ToDo: Finish!
}

void MatrixCompositeStd::add_vector(
	size_type                      row_offset
	,size_type                     col_offset
	,value_type                    beta
	,const Range1D                 &rng_G
	,const VectorWithOp            *v
	,const release_resource_ptr_t  &v_release
	,BLAS_Cpp::Transp              v_trans
	)
{
	fully_constructed_ = false;
	assert(0); // ToDo: Finish!
}

void MatrixCompositeStd::add_vector(
	size_type                      row_offset
	,size_type                     col_offset
	,value_type                    beta
	,const VectorWithOp            *v
	,const release_resource_ptr_t  &v_release
	,BLAS_Cpp::Transp              v_trans
	)
{
	assert( beta != 0.0 );
	assert( v != NULL );
	fully_constructed_ = false;
	if( v_trans == BLAS_Cpp::no_trans ) {
		assert( row_offset + v->dim() <= rows_ );
		assert( col_offset + 1 <= cols_ );
	}
	else {
		assert( row_offset + 1 <= rows_ );
		assert( col_offset + v->dim() <= cols_ );
	}

	vector_list_.push_back(
		SubVectorEntry(
			row_offset+1,col_offset+1,beta,Range1D()
			,NULL,NULL,BLAS_Cpp::no_trans
			,v,v_release,v_trans ) );
}

void MatrixCompositeStd::remove_vector( vector_list_t::iterator itr )
{
	fully_constructed_ = false;
	vector_list_.erase(itr);
}

void MatrixCompositeStd::add_matrix(
	size_type                      row_offset
	,size_type                     col_offset
	,value_type                    alpha
	,const GenPermMatrixSlice      *P
	,const release_resource_ptr_t  &P_release
	,BLAS_Cpp::Transp              P_trans
	,const MatrixWithOp            *A
	,const release_resource_ptr_t  &A_release
	,BLAS_Cpp::Transp              A_trans
	,const GenPermMatrixSlice      *Q
	,const release_resource_ptr_t  &Q_release
	,BLAS_Cpp::Transp              Q_trans
	)
{
	fully_constructed_ = false;
	assert(0); // ToDo: Finish!
}

void MatrixCompositeStd::add_matrix(
	size_type                      row_offset
	,size_type                     col_offset
	,value_type                    alpha
	,const Range1D                 &rng_P
	,const MatrixWithOp            *A
	,const release_resource_ptr_t  &A_release
	,BLAS_Cpp::Transp              A_trans
	,const GenPermMatrixSlice      *Q
	,const release_resource_ptr_t  &Q_release
	,BLAS_Cpp::Transp              Q_trans
	)
{
	fully_constructed_ = false;
	assert(0); // ToDo: Finish!
}

void MatrixCompositeStd::add_matrix(
	size_type                      row_offset
	,size_type                     col_offset
	,value_type                    alpha
	,const GenPermMatrixSlice      *P
	,const release_resource_ptr_t  &P_release
	,BLAS_Cpp::Transp              P_trans
	,const MatrixWithOp            *A
	,const release_resource_ptr_t  &A_release
	,BLAS_Cpp::Transp              A_trans
	,const Range1D                 &rng_Q
	)
{
	fully_constructed_ = false;
	assert(0); // ToDo: Finish!
}

void MatrixCompositeStd::add_matrix(
	size_type                      row_offset
	,size_type                     col_offset
	,value_type                    alpha
	,const Range1D                 &rng_P_in
	,const MatrixWithOp            *A
	,const release_resource_ptr_t  &A_release
	,BLAS_Cpp::Transp              A_trans
	,const Range1D                 &rng_Q_in
	)
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	using RangePack::full_range;

	THROW_EXCEPTION(
		alpha == 0.0, std::invalid_argument
		,"MatrixCompositeStd::add_matrix(...) : Error!" );
	THROW_EXCEPTION(
		A == NULL, std::invalid_argument
		,"MatrixCompositeStd::add_matrix(...) : Error!" );

	const size_type
		A_rows   = A->rows(),
		A_cols   = A->cols(),
		opA_rows = rows(A_rows,A_cols,A_trans),
		opA_cols = cols(A_rows,A_cols,A_trans);
	const Range1D
		rng_P = full_range(rng_P_in,1,opA_rows),
		rng_Q = full_range(rng_Q_in,1,opA_cols);
	const size_type
		opPopAopQ_rows = rng_P.size(),
		opPopAopQ_cols = rng_Q.size();

	THROW_EXCEPTION(
		row_offset + opPopAopQ_rows > rows_, std::invalid_argument
		,"MatrixCompositeStd::add_matrix(...) : Error!" );
	THROW_EXCEPTION(
		col_offset + opPopAopQ_cols > cols_, std::invalid_argument
		,"MatrixCompositeStd::add_matrix(...) : Error!" );

	fully_constructed_ = false;
	
	matrix_list_.push_back(
		SubMatrixEntry(
			row_offset+1,row_offset+opPopAopQ_rows
			,col_offset+1,col_offset+opPopAopQ_rows
			,alpha
			,rng_P
			,NULL,NULL,BLAS_Cpp::no_trans
			,A,A_release,A_trans
			,rng_Q
			,NULL,NULL,BLAS_Cpp::no_trans
			)
		);

	// ToDo: We could create a sub-view of the matrix using rng_P and rng_Q
	// and then store the sub-view in the data structure?  However, this
	// would confuse an external client who wanted to iterate through
	// the matrix list using the public iterators.
}

void MatrixCompositeStd::add_matrix(
	size_type                      row_offset
	,size_type                     col_offset
	,value_type                    alpha
	,const MatrixWithOp            *A
	,const release_resource_ptr_t  &A_release
	,BLAS_Cpp::Transp              A_trans
	)
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;

	THROW_EXCEPTION(
		alpha == 0.0, std::invalid_argument
		,"MatrixCompositeStd::add_matrix(...) : Error!" );
	THROW_EXCEPTION(
		A == NULL, std::invalid_argument
		,"MatrixCompositeStd::add_matrix(...) : Error!" );

	const size_type
		A_rows   = A->rows(),
		A_cols   = A->cols(),
		opA_rows = rows(A_rows,A_cols,A_trans),
		opA_cols = cols(A_rows,A_cols,A_trans);

	THROW_EXCEPTION(
		row_offset + opA_rows > rows_, std::invalid_argument
		,"MatrixCompositeStd::add_matrix(...) : Error!" );
	THROW_EXCEPTION(
		col_offset + opA_cols > cols_, std::invalid_argument
		,"MatrixCompositeStd::add_matrix(...) : Error!" );

	fully_constructed_ = false;

	matrix_list_.push_back(
		SubMatrixEntry(
			row_offset+1,row_offset+opA_rows
			,col_offset+1,col_offset+opA_cols
			,alpha
			,Range1D()
			,NULL,NULL,BLAS_Cpp::no_trans
			,A,A_release,A_trans
			,Range1D()
			,NULL,NULL,BLAS_Cpp::no_trans
			)
		);
}

void MatrixCompositeStd::add_matrix(
	size_type                      row_offset
	,size_type                     col_offset
	,value_type                    alpha
	,const GenPermMatrixSlice      *P
	,const release_resource_ptr_t  &P_release
	,BLAS_Cpp::Transp              P_trans
	)
{
	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;

	assert( alpha != 0.0 );
	assert( P != NULL );

	fully_constructed_ = false;

	const size_type
		P_rows   = P->rows(),
		P_cols   = P->cols(),
		opP_rows = rows(P_rows,P_cols,P_trans),
		opP_cols = cols(P_rows,P_cols,P_trans);

	assert( row_offset + opP_rows <= rows_ );
	assert( col_offset + opP_cols <= cols_ );

	matrix_list_.push_back(
		SubMatrixEntry(
			row_offset+1,row_offset+opP_rows,col_offset+1,col_offset+opP_cols,alpha
			,Range1D::Invalid
			,new GenPermMatrixSlice(*P),P_release,P_trans
			,NULL,NULL,BLAS_Cpp::no_trans
			,Range1D()
			,NULL,NULL,BLAS_Cpp::no_trans
			)
		);
}

void MatrixCompositeStd::remove_matrix( matrix_list_t::iterator itr )
{
	fully_constructed_ = false;
	matrix_list_.erase(itr);
}

void MatrixCompositeStd::finish_construction(
	const VectorSpace::space_ptr_t&  space_rows
	,const VectorSpace::space_ptr_t& space_cols
	)
{
	THROW_EXCEPTION(
		!space_rows.get(), std::invalid_argument
		,"MatrixCompositeStd::finish_construction(...): Error, space_rows.get() can not be NULL" );
	THROW_EXCEPTION(
		!space_cols.get(), std::invalid_argument
		,"MatrixCompositeStd::finish_construction(...): Error, space_cols.get() can not be NULL" );
	THROW_EXCEPTION(
		space_rows->dim() != cols_, std::invalid_argument
		,"MatrixCompositeStd::finish_construction(...): Error, space_rows->dim() = " << space_rows->dim()
		<< " != cols = " << cols_ << " where cols was passed to this->reinitialize(...)" );
	THROW_EXCEPTION(
		space_cols->dim() != rows_, std::invalid_argument
		,"MatrixCompositeStd::finish_construction(...): Error, space_colss->dim() = " << space_cols->dim()
		<< " != rows = " << rows_ << " where cols was passed to this->reinitialize(...)" );

	space_rows_ = space_rows;
	space_cols_ = space_cols;
	fully_constructed_ = true;
}

// Member access

MatrixCompositeStd::vector_list_t::iterator
MatrixCompositeStd::vectors_begin()
{
	return vector_list_.begin();
}

MatrixCompositeStd::vector_list_t::iterator
MatrixCompositeStd::vectors_end()
{
	return vector_list_.end();
}

MatrixCompositeStd::vector_list_t::const_iterator
MatrixCompositeStd::vectors_begin() const
{
	return vector_list_.begin();
}

MatrixCompositeStd::vector_list_t::const_iterator
MatrixCompositeStd::vectors_end() const
{
	return vector_list_.end();
}

MatrixCompositeStd::matrix_list_t::iterator
MatrixCompositeStd::matrices_begin()
{
	return matrix_list_.begin();
}

MatrixCompositeStd::matrix_list_t::iterator
MatrixCompositeStd::matrices_end()
{
	return matrix_list_.end();
}

MatrixCompositeStd::matrix_list_t::const_iterator
MatrixCompositeStd::matrices_begin() const
{
	return matrix_list_.begin();
}

MatrixCompositeStd::matrix_list_t::const_iterator
MatrixCompositeStd::matrices_end() const
{
	return matrix_list_.end();
}

// Overridden from Matrix

size_type MatrixCompositeStd::rows() const
{
	return fully_constructed_ ? rows_ : 0;
}

size_type MatrixCompositeStd::cols() const
{
	return fully_constructed_ ? cols_ : 0;
}

// Overridden from MatrixWithOp

const VectorSpace& MatrixCompositeStd::space_rows() const
{
	assert_fully_constructed();
	return *space_rows_;
}

const VectorSpace& MatrixCompositeStd::space_cols() const
{
	assert_fully_constructed();
	return *space_cols_;
}

MatrixWithOp::mat_ptr_t
MatrixCompositeStd::sub_view(const Range1D& row_rng, const Range1D& col_rng) const
{
	assert_fully_constructed();
	// For this initial implementation we will just look through the list of matrices
	// and find an exact match.  If we don't find it we will just return the result
	// from the default implementation of this method.

	// ToDo: Implement!

	// Just return the default sub-view
	return MatrixWithOp::sub_view(row_rng,col_rng);
}

void MatrixCompositeStd::Vp_StMtV(
	VectorWithOpMutable* y, value_type a, BLAS_Cpp::Transp M_trans
	, const VectorWithOp& x, value_type b
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixCompositeStd::Vp_StMtV(...VectorSlice...)" );
#endif
	assert_fully_constructed();
	Vp_StMtV_imp(y,a,rows_,cols_,matrix_list_,vector_list_,M_trans,x,b);
}

void MatrixCompositeStd::Vp_StMtV(
	VectorWithOpMutable* y, value_type a, BLAS_Cpp::Transp M_trans
	, const SpVectorSlice& x, value_type b
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixCompositeStd::Vp_StMtV(...SpVectorSlice...)" );
#endif
	assert_fully_constructed();
	Vp_StMtV_imp(y,a,rows_,cols_,matrix_list_,vector_list_,M_trans,x,b);
}

void MatrixCompositeStd::Vp_StPtMtV(
	VectorWithOpMutable* y, value_type a
	, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp M_trans
	, const VectorWithOp& x, value_type b
	) const
{
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixCompositeStd::Vp_StPtMtV(...VectorSlice...)" );
#endif
	assert_fully_constructed();
	MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,M_trans,x,b); // ToDo: Implement when needed!
}

void MatrixCompositeStd::Vp_StPtMtV(
	VectorWithOpMutable* y, value_type a
	, const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
	, BLAS_Cpp::Transp M_trans
	, const SpVectorSlice& x, value_type b
	) const
{	
#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "MatrixCompositeStd::Vp_StPtMtV(...SpVectorSlice...)" );
#endif
	assert_fully_constructed();
	MatrixWithOp::Vp_StPtMtV(y,a,P,P_trans,M_trans,x,b); // ToDo: Implement when needed!
}

// private

void MatrixCompositeStd::assert_fully_constructed() const
{
	const bool fully_constructed = fully_constructed_;
	THROW_EXCEPTION(
		!fully_constructed, std::logic_error
		,"MatrixCompositeStd::assert_fully_constructed() : Error, not fully constructed!");
}

} // end namespace AbstractLinAlgPack