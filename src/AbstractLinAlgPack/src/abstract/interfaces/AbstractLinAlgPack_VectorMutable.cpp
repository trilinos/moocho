// //////////////////////////////////////////////////////////////////////
// VectorMutable.cpp
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

#include "VectorMutable.hpp"
#include "VectorMutableSubView.hpp"
#include "VectorSpace.hpp"
#include "RTOpStdOpsLib/src/RTOp_TOp_assign_scalar.h"
#include "RTOpStdOpsLib/src/RTOp_TOp_assign_vectors.h"
#include "RTOpStdOpsLib/src/RTOp_TOp_axpy.h"
#include "RTOpStdOpsLib/src/RTOp_TOp_set_ele.h"
#include "RTOpStdOpsLib/src/RTOp_TOp_set_sub_vector.h"
#include "RTOpPack/src/RTOpCppC.hpp"
#include "Range1D.hpp"
#include "ThrowException.hpp"

namespace {

// vector scalar assignment operator
static RTOpPack::RTOpC          assign_scalar_op;
// vector assignment operator
static RTOpPack::RTOpC          assign_vec_op;
// set element operator
static RTOpPack::RTOpC          set_ele_op;
// set sub-vector operator
static RTOpPack::RTOpC          set_sub_vector_op;
// axpy operator
static RTOpPack::RTOpC          axpy_op;

// Simple class for an object that will initialize the RTOp_Server for get_ele operator.
class init_rtop_server_t {
public:
	init_rtop_server_t() {
		// DVector scalar assignment operator
		if(0>RTOp_TOp_assign_scalar_construct( 0.0, &assign_scalar_op.op() ))
			assert(0);
		if(0>RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_assign_scalar_name
			   ,&RTOp_TOp_assign_scalar_vtbl
			   ))
			assert(0);
		// DVector assignment operator
		if(0>RTOp_TOp_assign_vectors_construct( &assign_vec_op.op() ))
			assert(0);
		if(0>RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_assign_vectors_name
			   ,&RTOp_TOp_assign_vectors_vtbl
			   ))
			assert(0);
		// Set element operator
		if(0>RTOp_TOp_set_ele_construct( 0, 0.0, &set_ele_op.op() ))
			assert(0);
		if(0>RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_set_ele_name
			   ,&RTOp_TOp_set_ele_vtbl
			   ))
			assert(0);	
		// Set sub-vector operator
		RTOp_SubVector sub_vec;
		RTOp_sub_vector_null(&sub_vec);
		if(0>RTOp_TOp_set_sub_vector_construct( &sub_vec, &set_sub_vector_op.op() ))
			assert(0);
		if(0>RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_set_sub_vector_name
			   ,&RTOp_TOp_set_sub_vector_vtbl
			   ))
			assert(0);	
		// axpy operator
		if(0>RTOp_TOp_axpy_construct( 0.0, &axpy_op.op() ))
			assert(0);
		if(0>RTOp_Server_add_op_name_vtbl(
			   RTOp_TOp_axpy_name
			   ,&RTOp_TOp_axpy_vtbl
			   ))
			assert(0);	
	}
}; 

// When the program starts, this object will be created and the RTOp_Server object will
// be initialized before main() gets underway!
init_rtop_server_t  init_rtop_server;

} // end namespace

namespace AbstractLinAlgPack {

// VectorMutable

VectorMutable& VectorMutable::operator=(value_type alpha)
{
	if(0!=RTOp_TOp_assign_scalar_set_alpha( alpha, &assign_scalar_op.op() ))
		assert(0);
	this->apply_transformation(assign_scalar_op,0,NULL,0,NULL,RTOp_REDUCT_OBJ_NULL);
	return *this;
}

VectorMutable& VectorMutable::operator=(const Vector& vec)
{
	if( dynamic_cast<const void*>(&vec) == dynamic_cast<const void*>(this) )
		return *this; // Assignment to self!
	const int num_vecs = 1;
	const Vector*
		vec_args[1] = { &vec };
	this->apply_transformation(assign_vec_op,num_vecs,vec_args,0,NULL,RTOp_REDUCT_OBJ_NULL);
	return *this;
}

VectorMutable& VectorMutable::operator=(const VectorMutable& vec)
{
	return this->operator=(static_cast<const Vector&>(vec));
}

void VectorMutable::set_ele( index_type i, value_type alpha )
{
	if(0!=RTOp_TOp_set_ele_set_i_alpha( i, alpha, &set_ele_op.op() ))
		assert(0);
	this->apply_transformation(
		set_ele_op,0,NULL,0,NULL,RTOp_REDUCT_OBJ_NULL
		,i,1,i-1 // first_ele, sub_dim, global_offset
		);
}

VectorMutable::vec_mut_ptr_t
VectorMutable::sub_view( const Range1D& rng_in )
{
	namespace rcp = MemMngPack;
	const index_type dim = this->dim();
	const Range1D    rng = rng_in.full_range() ? Range1D(1,dim) : rng_in;
#ifdef _DEBUG
	THROW_EXCEPTION(
		rng.ubound() > dim, std::out_of_range
		,"VectorMutable::sub_view(rng): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1,"<<dim<<"]" );
#endif	
	if( rng.lbound() == 1 && rng.ubound() == dim )
		return rcp::rcp( this, false );
	return rcp::rcp(
		new VectorMutableSubView(
			rcp::rcp( this, false )
			,rng ) );
}

void VectorMutable::zero()
{
	this->operator=(0.0);
}

void VectorMutable::axpy( value_type alpha, const Vector& x )
{
	if( 0!=RTOp_TOp_axpy_set_alpha( alpha, &axpy_op.op() ) )
		assert(0);
	const int num_vecs = 1;
	const Vector*
		vec_args[1] = { dynamic_cast<const Vector*>(&x) };
	if( vec_args[0] == NULL )
		throw VectorSpace::IncompatibleVectorSpaces(
			"Vector::axpy(alpha,x): Error, x is not of type Vector!" );
	this->apply_transformation(axpy_op,num_vecs,vec_args,0,NULL,RTOp_REDUCT_OBJ_NULL);
}

void VectorMutable::get_sub_vector(
	const Range1D& rng, RTOp_MutableSubVector* sub_vec )
{
	// Here we get a copy of the data for the sub-vector that the client will
	// modify.  We must later commit these changes to the actual vector
	// when the client calls commit_sub_vector(...).
	// Note, this implementation is very dependent on the behavior of the default
	// implementation of Vector::get_sub_vector(...) and
	// Vector::set_sub_vector(...)!
	RTOp_SubVector _sub_vec;
	RTOp_sub_vector_null( &_sub_vec );
	Vector::get_sub_vector(
		rng
		,DENSE
		,&_sub_vec
		);
	RTOp_mutable_sub_vector(
		_sub_vec.global_offset
		,_sub_vec.sub_dim
		,const_cast<value_type*>(_sub_vec.values)
		,_sub_vec.values_stride
		,sub_vec
		);
}

void VectorMutable::commit_sub_vector( RTOp_MutableSubVector* sub_vec )
{
	RTOp_SubVector _sub_vec;
	RTOp_sub_vector_dense(
		sub_vec->global_offset
		,sub_vec->sub_dim
		,sub_vec->values
		,sub_vec->values_stride
		,&_sub_vec
		);
	VectorMutable::set_sub_vector( _sub_vec ); // Commit the changes!
	Vector::free_sub_vector( &_sub_vec );      // Free the memory!
}

void VectorMutable::set_sub_vector( const RTOp_SubVector& sub_vec )
{
	if(0!=RTOp_TOp_set_sub_vector_set_sub_vec( &sub_vec, &set_sub_vector_op.op() ))
		assert(0);
	this->apply_transformation(
		set_sub_vector_op,0,NULL,0,NULL,RTOp_REDUCT_OBJ_NULL
		,sub_vec.global_offset+1,sub_vec.sub_dim,sub_vec.global_offset // first_ele, sub_dim, global_offset
		);
}

void VectorMutable::Vp_StMtV(
	value_type                       alpha
	,const GenPermMatrixSlice        &P
	,BLAS_Cpp::Transp                P_trans
	,const Vector              &x
	,value_type                      beta
	)
{
	assert(0); // ToDo: Implement using reduction or transformation operators
}

// Overridden from Vector

Vector::vec_ptr_t
VectorMutable::sub_view( const Range1D& rng ) const
{
	namespace rcp = MemMngPack;
	return const_cast<VectorMutable*>(this)->sub_view(rng);
}

// protected

void VectorMutable::finalize_apply_transformation(
	const size_t num_targ_vecs, VectorMutable** targ_vecs
	)
{
	this->has_changed();
	for( int k = 0; k < num_targ_vecs; ++k )
		targ_vecs[k]->has_changed();
}

} // end namespace AbstractLinAlgPack
