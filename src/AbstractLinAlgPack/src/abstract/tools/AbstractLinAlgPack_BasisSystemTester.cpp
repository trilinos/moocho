// ///////////////////////////////////////////////////////////
// BasisSystemTester.cpp
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

#include <math.h>

#include <ostream>

#include "AbstractLinAlgPack/include/BasisSystemTester.h"
#include "AbstractLinAlgPack/include/BasisSystem.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/include/MatrixCompositeStd.h"
#include "AbstractLinAlgPack/include/assert_print_nan_inf.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"

namespace AbstractLinAlgPack {

BasisSystemTester::BasisSystemTester(
	EPrintTestLevel  print_tests
	,bool            dump_all
	,bool            throw_exception
	,size_type       num_random_tests
	,value_type      warning_tol
	,value_type      error_tol
	)
	:print_tests_(print_tests)
	,dump_all_(dump_all)
	,throw_exception_(throw_exception)
	,num_random_tests_(num_random_tests)
	,warning_tol_(warning_tol)
	,error_tol_(error_tol)
{}

bool BasisSystemTester::test_basis_system(
	const BasisSystem               &bs
	,const MatrixWithOp             *Gc
	,const MatrixWithOp             *Gh
	,const MatrixWithOpNonsingular  *C
	,const MatrixWithOp             *N_in
	,const MatrixWithOp             *D
	,const MatrixWithOp             *GcUP
	,const MatrixWithOp             *GhUP
	,std::ostream                   *out
	)
{
	namespace rcp = MemMngPack;
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using AbstractLinAlgPack::sum;
	using AbstractLinAlgPack::dot;
	using AbstractLinAlgPack::Vp_StV;
	using AbstractLinAlgPack::assert_print_nan_inf;
	using AbstractLinAlgPack::random_vector;
	using LinAlgOpPack::V_StMtV;
	using LinAlgOpPack::V_MtV;
	using LinAlgOpPack::V_StV;
	using LinAlgOpPack::V_VpV;
	using LinAlgOpPack::Vp_V;
	
	// ToDo: Check the preconditions
	
	bool success = true, result, lresult, llresult;
	const value_type
		rand_y_l  = -1.0,
		rand_y_u  = 1.0,
		small_num = ::pow(std::numeric_limits<value_type>::epsilon(),0.25),
		alpha     = 2.0;
	
	EPrintTestLevel
		print_tests = ( this->print_tests() == PRINT_NOT_SELECTED ? PRINT_NONE : this->print_tests() );
	
	// Print the input?
	if( out && print_tests != PRINT_NONE ) {
		if( print_tests >= PRINT_BASIC )
			*out << "\n*************************************************"
				 << "\n*** BasisSystemTester::test_basis_system(...) ***"
				 << "\n*************************************************\n";
		if(dump_all()) {
			if(Gc)
				*out << "\nGc =\n"    << *Gc;
			if(Gh)
				*out << "\nGh =\n"    << *Gh;
			if(C)
				*out << "\nC =\n"     << *C;
			if(N_in)
				*out << "\nN =\n"     << *N_in;
			if(D)
				*out << "\nD =\n"     << *D;
			if(GcUP)
				*out << "\nGcUP =\n"  << *GcUP;
			if(GhUP)
				*out << "\nGhUP =\n"  << *GhUP;
		}
	}

	//
	// Check the dimensions of everything
	//

	const Range1D
		var_dep          = bs.var_dep(),
		var_indep        = bs.var_indep(),
		equ_decomp       = bs.equ_decomp(),
		equ_undecomp     = bs.equ_undecomp(),
		inequ_decomp     = bs.inequ_decomp(),
		inequ_undecomp   = bs.inequ_undecomp();

	if( out && print_tests >= PRINT_MORE ) {
		*out
			<< "\nbs.var_dep()        = ["<<var_dep.lbound()<<","<<var_dep.ubound()<<"]"
			<< "\nbs.var_indep( )     = ["<<var_indep.lbound()<<","<<var_indep.ubound()<<"]"
			<< "\nbs.equ_decomp()     = ["<<equ_decomp.lbound()<<","<<equ_decomp.ubound()<<"]"
			<< "\nbs.equ_undecomp()   = ["<<equ_undecomp.lbound()<<","<<equ_undecomp.ubound()<<"]"
			<< "\nbs.inequ_decomp()   = ["<<inequ_decomp.lbound()<<","<<inequ_decomp.ubound()<<"]"
			<< "\nbs.inequ_undecomp() = ["<<inequ_undecomp.lbound()<<","<<inequ_undecomp.ubound()<<"]"
			<< std::endl;
	}

	if( out && print_tests >= PRINT_BASIC )
		*out << "\n1) Check the partitioning ranges ...";
	lresult = true;

	if( out && print_tests >= PRINT_MORE )
		*out << "\n\n1.a) check: var_dep.size() != equ_decomp.size() + inequ_decomp.size() : ";
	result = var_dep.size() == equ_decomp.size() + inequ_decomp.size();
	if(out && print_tests >= PRINT_MORE)
		*out << ( result ? "passed" : "failed" );
	if(!result) lresult = false;

	if(Gc) {
		if( out && print_tests >= PRINT_MORE )
			*out << "\n1.b) check: var_dep.size() + var_indep.size() == Gc->rows() : ";
		result = var_dep.size() + var_indep.size() == Gc->rows();
		if(out && print_tests >= PRINT_MORE)
			*out << ( result ? "passed" : "failed" );
		if(!result) lresult = false;
	}	

	if(Gh) {
		if( out && print_tests >= PRINT_MORE )
			*out << "\n1.c) check: var_dep.size() + var_indep.size() == Gh->rows() : ";
		result = var_dep.size() + var_indep.size() == Gh->rows();
		if(out && print_tests >= PRINT_MORE)
			*out << ( result ? "passed" : "failed" );
		if(!result) lresult = false;
	}

	if(Gc) {
		if( out && print_tests >= PRINT_MORE )
			*out << "\n1.d) check: equ_decomp.size() + equ_undecomp.size() == Gc->cols() : ";
		result = equ_decomp.size() + equ_undecomp.size() == Gc->cols();
		if(out && print_tests >= PRINT_MORE)
			*out << ( result ? "passed" : "failed" );
		if(!result) lresult = false;
	}	

	if(Gh) {
		if( out && print_tests >= PRINT_MORE )
			*out << "\n1.e) check: inequ_decomp.size() + inequ_undecomp.size() == Gh->cols() : ";
		result = inequ_decomp.size() + inequ_undecomp.size() == Gh->cols();
		if(out && print_tests >= PRINT_MORE)
			*out << ( result ? "passed" : "failed" );
		if(!result) lresult = false;
	}	

	if(out && print_tests >= PRINT_MORE)
		*out << std::endl;

	if(!lresult) success = false;
	if( out && print_tests == PRINT_BASIC )
		*out << " : " << ( lresult ? "passed" : "failed" );

	// Create the N matrix if not input
	rcp::ref_count_ptr<const AbstractLinAlgPack::MatrixWithOp>
		N = rcp::rcp(N_in,false);
	if( (Gc || Gh) && C && N.get() == NULL ) {
		if(out && print_tests >= PRINT_BASIC)
			*out
				<< "\nCreating the matrix N since it was not input by the client ...";
		if(out && print_tests >= PRINT_MORE)
			*out
				<< std::endl;
		rcp::ref_count_ptr<AbstractLinAlgPack::MatrixCompositeStd>
			N_comp = rcp::rcp(new AbstractLinAlgPack::MatrixCompositeStd(var_dep.size(),var_indep.size()));
		if( equ_decomp.size() )
			N_comp->add_matrix( 0, 0, 1.0, equ_decomp, Gc, rcp::null, BLAS_Cpp::trans, var_indep );
		if( Gh && inequ_decomp.size() )
			N_comp->add_matrix( equ_decomp.size(), 0, 1.0, inequ_decomp, Gh, rcp::null, BLAS_Cpp::trans, var_indep );
		N_comp->finish_construction(
			Gc->space_rows().sub_space(equ_decomp)->clone()
			,Gc->space_cols().sub_space(var_indep)->clone()
			);
		if( out && dump_all() )
			*out << "\nN =\n" << *N_comp;
		N = N_comp;
	}

	// Create the other auxillary matrix objects
	if( equ_undecomp.size() ) {
		assert(0); // ToDo: Create matrix objects for Gc(var_dep,equ_undecomp) and Gc(var_indep,equ_undecomp)
	}
	if( inequ_undecomp.size() ) {
		assert(0); // ToDo: Create matrix objects for Gh(var_dep,inequ_undecomp) and Gh(var_indep,inequ_undecomp)
	}

	//
	// Perform the tests
	//

	if( C && N.get() ) {

		if(out && print_tests >= PRINT_BASIC)
			*out
				<< "\n2) Check the compatibility of the vector spaces for C, N, D, Gc and Gh ...";
		lresult = true;

		if(out && print_tests >= PRINT_MORE)
			*out
				<< "\n2.a) Check consistency of the vector spaces for:"
				<< "\n    C.space_cols() == N.space_cols()";
		llresult = true;
		if(out && print_tests >= PRINT_ALL)
 			*out << "\n\n2.a.1) C->space_cols().is_compatible(N->space_cols()) == true : ";
		result = C->space_cols().is_compatible(N->space_cols());	
		if(out && print_tests >= PRINT_ALL)
		 	*out << ( result ? "passed" : "failed" )
				 << std::endl;
		if(!result) llresult = false;
		if(!llresult) lresult = false;
		if( out && print_tests == PRINT_MORE )
		*out << " : " << ( llresult ? "passed" : "failed" );

		if(D) {
			if(out && print_tests >= PRINT_MORE)
				*out
					<< "\n2.b) Check consistency of the vector spaces for:"
					<< "\n    D.space_cols() == C.space_cols() and D.space_rows() == N.space_rows()";
			llresult = true;
			if(out && print_tests >= PRINT_ALL)
				*out << "\n2.b.1) D->space_cols().is_compatible(C->space_cols()) == true : ";
			result = D->space_cols().is_compatible(C->space_cols());
			if(out && print_tests >= PRINT_ALL)
				*out << ( result ? "passed" : "failed" );
			if(!result) llresult = false;
			if(out && print_tests >= PRINT_ALL)
				*out << "\n2.b.2) D->space_rows().is_compatible(N->space_rows()) == true : ";
			result = D->space_rows().is_compatible(N->space_rows());
			if(out && print_tests >= PRINT_ALL)
				*out << ( result ? "passed" : "failed" )
					 << std::endl;
			if(!result) llresult = false;
			if(!llresult) lresult = false;
			if( out && print_tests == PRINT_MORE )
				*out << " : " << ( llresult ? "passed" : "failed" );
		}

		if(out && print_tests >= PRINT_MORE)
			*out
				<< "\n2.c) Check consistency of the vector spaces for:"
				<< "\n    [ Gc'(equ_decomp,  var_dep) ]"
				<< "\n    [ Gh'(inequ_decomp,var_dep) ] == C";
		llresult = true;
		if( equ_decomp.size() ) {
			if(out && print_tests >= PRINT_ALL)
				*out << "\n2.c.1) Gc->space_rows().sub_space(equ_decomp)->is_compatible(*C->space_cols().sub_space(equ_decomp)) == true : ";
			result = Gc->space_rows().sub_space(equ_decomp)->is_compatible(*C->space_cols().sub_space(equ_decomp));
			if(out && print_tests >= PRINT_ALL)
				*out << ( result ? "passed" : "failed" );
			if(!result) llresult = false;
			if(out && print_tests >= PRINT_ALL)
				*out << "\n2.c.2) Gc->space_cols().sub_space(var_dep)->is_compatible(C->space_rows()) == true : ";
			result = Gc->space_cols().sub_space(var_dep)->is_compatible(C->space_rows());
			if(out && print_tests >= PRINT_ALL)
				*out << ( result ? "passed" : "failed" );
			if(!result) llresult = false;
		}
		if( inequ_decomp.size() ) {
			if(out && print_tests >= PRINT_ALL)
				*out << "\n2.c.3) Gh->space_rows().sub_space(inequ_decomp)->is_compatible(*C->space_cols().sub_space(inequ_decomp)) == true : ";
			result = Gh->space_rows().sub_space(inequ_decomp)->is_compatible(*C->space_cols().sub_space(inequ_decomp));
			if(out && print_tests >= PRINT_ALL)
				*out << ( result ? "passed" : "failed" );
			if(!result) llresult = false;
			if(out && print_tests >= PRINT_ALL)
				*out << "\n2.c.4) Gh->space_cols().sub_space(var_dep)->is_compatible(C->space_rows()) == true : ";
			result = Gh->space_cols().sub_space(var_dep)->is_compatible(C->space_rows());
			if(out && print_tests >= PRINT_ALL)
				*out << ( result ? "passed" : "failed" );
			if(!result) llresult = false;
		}
		if(out && print_tests >= PRINT_ALL)
			*out << std::endl;
		if(!llresult) lresult = false;
		if( out && print_tests == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" );

		if(out && print_tests >= PRINT_MORE)
			*out
				<< "\n2.d) Check consistency of the vector spaces for:"
				<< "\n    [ Gc'(equ_decomp,  var_indep) ]"
				<< "\n    [ Gh'(inequ_decomp,var_indep) ] == N";
		llresult = true;
		if( equ_decomp.size() ) {
			if(out && print_tests >= PRINT_ALL)
				*out << "\n2.d.1) Gc->space_rows().sub_space(equ_decomp)->is_compatible(*N->space_cols().sub_space(equ_decomp)) == true : ";
			result = Gc->space_rows().sub_space(equ_decomp)->is_compatible(*N->space_cols().sub_space(equ_decomp));
			if(out && print_tests >= PRINT_ALL)
				*out << ( result ? "passed" : "failed" );
			if(!result) llresult = false;
			if(out && print_tests >= PRINT_ALL)
				*out << "\n2.d.2) Gc->space_cols().sub_space(var_indep)->is_compatible(N->space_rows()) == true : ";
			result = Gc->space_cols().sub_space(var_indep)->is_compatible(N->space_rows());
			if(out && print_tests >= PRINT_ALL)
				*out << ( result ? "passed" : "failed" );
			if(!result) llresult = false;
		}
		if( inequ_decomp.size() ) {
			if(out && print_tests >= PRINT_ALL)
				*out << "\n2.d.3) Gh->space_rows().sub_space(inequ_decomp)->is_compatible(*N->space_cols().sub_space(inequ_decomp)) == true : ";
			result = Gh->space_rows().sub_space(inequ_decomp)->is_compatible(*N->space_cols().sub_space(inequ_decomp));
			if(out && print_tests >= PRINT_ALL)
				*out << ( result ? "passed" : "failed" );
			if(!result) llresult = false;
			if(out && print_tests >= PRINT_ALL)
				*out << "\n2.d.4) Gh->space_cols().sub_space(var_indep)->is_compatible(N->space_rows()) == true : ";
			result = Gh->space_cols().sub_space(var_indep)->is_compatible(N->space_rows());
			if(out && print_tests >= PRINT_ALL)
				*out << ( result ? "passed" : "failed" );
			if(!result) llresult = false;
		}
		if(out && print_tests >= PRINT_ALL)
			*out << std::endl;
		if(!llresult) lresult = false;
		if( out && print_tests == PRINT_MORE )
			*out << " : " << ( llresult ? "passed" : "failed" )
				 << std::endl;

		if(!lresult) success = false;
		if( out && print_tests == PRINT_BASIC )
			*out << " : " << ( lresult ? "passed" : "failed" );

		if(out && print_tests >= PRINT_BASIC)
			*out
				<< "\n3) Check the compatibility of the matrices C, N, D, Gc and Gh numerically ...";

		if(out && print_tests >= PRINT_MORE)
			*out
				<< std::endl
				<< "\n3.a) Check consistency of:"
				<< "\n                [ Gc'(equ_decomp,  var_dep) Gc'(equ_decomp,  var_indep) ]"
				<< "\n    op ( alpha* [ Gh'(inequ_decomp,var_dep) Gh'(inequ_decomp,var_indep) ] ) * v"
				<< "\n         \\______________________________________________________________/"
				<< "\n                                        A"
				<< "\n    ==  op( alpha*[ C  N ] ) * v"
				<< "\n            \\____________/"
				<< "\n                   B"
				<< "\nfor random vectors v ...";
		{

			VectorSpace::vec_mut_ptr_t
				Gc_v_x    = ( Gc ? Gc->space_cols().create_member() : rcp::null ),
				Gc_v_c    = ( Gc ? Gc->space_rows().create_member() : rcp::null ),
				Gh_v_x    = ( Gh ? Gh->space_cols().create_member() : rcp::null ),
				Gh_v_h    = ( Gh ? Gh->space_rows().create_member() : rcp::null ),
				C_v_xD    = C->space_rows().create_member(),
				C_v_chD   = C->space_cols().create_member(),
				N_v_xI    = N->space_rows().create_member(),
				N_v_chD   = N->space_cols().create_member(),
				v_x       = ( Gc ? Gc->space_cols().create_member() : Gh->space_cols().create_member() ),
				v_x_tmp   = v_x->space().create_member(),
				v_chD     = C_v_xD->space().create_member(),
				v_chD_tmp = v_chD->space().create_member();

			if(out && print_tests >= PRINT_MORE)
				*out << "\n\n3.a.1) Testing non-transposed A*v == B*v ...";
			if(out && print_tests > PRINT_MORE)
				*out << std::endl;
			llresult = true;
 			{for( int k = 1; k <= num_random_tests(); ++k ) {
				random_vector( rand_y_l, rand_y_u, v_x.get() );
			 	if(out && print_tests >= PRINT_ALL) {
					*out
						<< "\n3.a.1."<<k<<") random vector " << k << " ( ||v_x||_1 / n = " << (v_x->norm_1() / v_x->dim()) << " )\n";
					if(dump_all() && print_tests >= PRINT_ALL)
						*out << "\nv_x =\n" << *v_x;
				}
				if(Gc && equ_decomp.size()) {
					V_StMtV( Gc_v_c.get(), alpha, *Gc, trans, *v_x );
					*v_chD_tmp->sub_view(equ_decomp)
						= *Gc_v_c->sub_view(equ_decomp);
				}
				if(Gh && inequ_decomp.size()) {
					V_StMtV( Gh_v_h.get(), alpha, *Gh, trans, *v_x );
					*v_chD_tmp->sub_view(equ_decomp.size() + inequ_decomp)
						= *Gh_v_h->sub_view(inequ_decomp);
				}
				V_StMtV( C_v_chD.get(), alpha, *C, no_trans, *v_x->sub_view(var_dep) );
				V_StMtV( N_v_chD.get(), alpha, *N, no_trans, *v_x->sub_view(var_indep) );
				V_VpV( v_chD.get(), *C_v_chD, *N_v_chD );
				const value_type
					sum_Bv  = sum(*v_chD),
					sum_Av  = sum(*v_chD_tmp);
				assert_print_nan_inf(sum_Bv, "sum(B*v_x)",true,out);
				assert_print_nan_inf(sum_Av, "sum(A*v_x)",true,out);
				const value_type
					calc_err = ::fabs( ( sum_Av - sum_Bv )
									   /( ::fabs(sum_Av) + ::fabs(sum_Bv) + small_num ) );
				if(out && print_tests >= PRINT_ALL)
					*out
						<< "\nrel_err(sum(A*v_x),sum(B*v_x)) = "
						<< "rel_err(" << sum_Av << "," << sum_Bv << ") = "
						<< calc_err << std::endl;
				if( calc_err >= warning_tol() ) {
					if(out && print_tests >= PRINT_ALL)
						*out
							<< std::endl
							<< ( calc_err >= error_tol() ? "Error" : "Warning" )
							<< ", rel_err(sum(A*v_x),sum(B*v_x)) = "
							<< "rel_err(" << sum_Av << "," << sum_Bv << ") = "
							<< calc_err
							<< " exceeded "
							<< ( calc_err >= error_tol() ? "error_tol" : "warning_tol" )
							<< " = "
							<< ( calc_err >= error_tol() ? error_tol() : warning_tol() )
							<< std::endl;
					if(calc_err >= error_tol()) {
						if(dump_all() && print_tests >= PRINT_ALL) {
							*out << "\nalpha = " << alpha << std::endl;
							*out << "\nv_x =\n"              << *v_x;
							if(Gc_v_c.get() && equ_decomp.size())
								*out << "\nalpha*Gc*v_x =\n" << *Gc_v_c;
							if(Gh_v_h.get() && inequ_decomp.size())
								*out << "\nalpha*Gh*v_x =\n" << *Gh_v_h;
							*out << "A*v =\n"                << *v_chD_tmp;
							*out << "\nalpha*C*v_x =\n"      << *C_v_chD;
							*out << "\nalpha*N*v_x =\n"      << *N_v_chD;
							*out << "\nB*v =\n"              << *v_chD;
						}
						llresult = false;
					}
				}
			}}
			if(!llresult) lresult = false;
			if( out && print_tests == PRINT_MORE )
				*out << " : " << ( llresult ? "passed" : "failed" )
					 << std::endl;

			if(out && print_tests >= PRINT_MORE)
				*out << "\n3.a.2) Testing transposed A'*v == B'*v ...";
			if(out && print_tests > PRINT_MORE)
				*out << std::endl;
			llresult = true;
			{for( int k = 1; k <= num_random_tests(); ++k ) {
				random_vector( rand_y_l, rand_y_u, v_chD.get() );
			 	if(out && print_tests >= PRINT_ALL) {
					*out
						<< "\n3.a.2."<<k<<") random vector " << k << " ( ||v_chD||_1 / n = " << (v_chD->norm_1() / v_chD->dim()) << " )\n";
					if(dump_all() && print_tests >= PRINT_ALL)
						*out << "\nv_chD =\n" << *v_chD;
				}
				*v_x_tmp = 0.0;
				if(Gc && equ_decomp.size()) {
					*Gc_v_c->sub_view(equ_decomp) = *v_chD->sub_view(equ_decomp);
					if(equ_undecomp.size())
						*Gc_v_c->sub_view(equ_undecomp) = 0.0;
					V_StMtV( Gc_v_x.get(), alpha, *Gc, no_trans, *Gc_v_c );
					Vp_V( v_x_tmp.get(), *Gc_v_x );
				}
				if(Gh && inequ_decomp.size()) {
					*Gh_v_h->sub_view(inequ_decomp) = *v_chD->sub_view(equ_decomp.size()+inequ_decomp);
					if(inequ_undecomp.size())
						*Gh_v_h->sub_view(inequ_undecomp) = 0.0;
					V_StMtV( Gh_v_x.get(), alpha, *Gh, no_trans, *Gh_v_h );
					Vp_V( v_x_tmp.get(), *Gh_v_x );
				}
				V_StMtV( C_v_xD.get(), alpha, *C, trans, *v_chD );
				*v_x->sub_view(var_dep) = *C_v_xD;
				V_StMtV( N_v_xI.get(), alpha, *N, trans, *v_chD );
				*v_x->sub_view(var_indep) = *N_v_xI;
				const value_type
					sum_BTv  = sum(*v_x),
					sum_ATv  = sum(*v_x_tmp);
				assert_print_nan_inf(sum_BTv, "sum(B'*v_chD)",true,out);
				assert_print_nan_inf(sum_ATv, "sum(A'*v_chD)",true,out);
				const value_type
					calc_err = ::fabs( ( sum_ATv - sum_BTv )
									   /( ::fabs(sum_ATv) + ::fabs(sum_BTv) + small_num ) );
				if(out && print_tests >= PRINT_ALL)
					*out
						<< "\nrel_err(sum(A'*v_chD),sum(B'*v_chD)) = "
						<< "rel_err(" << sum_ATv << "," << sum_BTv << ") = "
						<< calc_err << std::endl;
				if( calc_err >= warning_tol() ) {
					if(out && print_tests >= PRINT_ALL)
						*out
							<< std::endl
							<< ( calc_err >= error_tol() ? "Error" : "Warning" )
							<< ", rel_err(sum(A'*v_chD),sum(B'*v_chD)) = "
							<< "rel_err(" << sum_ATv << "," << sum_BTv << ") = "
							<< calc_err << std::endl
							<< " exceeded "
							<< ( calc_err >= error_tol() ? "error_tol" : "warning_tol" )
							<< " = "
							<< ( calc_err >= error_tol() ? error_tol() : warning_tol() )
							<< std::endl;
					if(calc_err >= error_tol()) {
						if(dump_all() && print_tests >= PRINT_ALL) {
							*out << "\nalpha = " << alpha << std::endl;
							*out << "\nv_chD =\n"            << *v_chD;
							if(Gc_v_x.get() && equ_decomp.size()) {
								*out << "\nGc_v_c =\n"       << *Gc_v_c;
								*out << "\nalpha*Gc'*[v_chD(equ_decomp); 0] =\n"
									 << *Gc_v_x;
							}
							if(Gh_v_x.get() && inequ_decomp.size()) {
								*out << "\nGh_v_h =\n"       << *Gh_v_h;
								*out << "\nalpha*Gh'*[v_chD(equ_decomp.size()+inequ_decomp); 0] =\n"
									 << *Gh_v_x;
							}
							*out << "A'*v =\n"               << *v_x_tmp;
							*out << "\nalpha*C*v_chD =\n"    << *C_v_xD;
							*out << "\nalpha*N*v_chD =\n"    << *N_v_xI;
							*out << "\nB'*v =\n"             << *v_x;
						}
						llresult = false;
					}
				}
			}}
			if(!llresult) lresult = false;
			if( out && print_tests == PRINT_MORE )
				*out << " : " << ( llresult ? "passed" : "failed" )
					 << std::endl;
		}
		
		if(out && print_tests >= PRINT_MORE)
			*out
				<< "\n3.b) Check consistency of:"
				<< "\n    alpha*op(C)*(op(inv(C)) * v) == alpha*v"
				<< "\nfor random vectors v ...";
		{
			VectorSpace::vec_mut_ptr_t
				v_xD      = C->space_rows().create_member(),
				v_xD_tmp  = C->space_rows().create_member(),
				v_chD     = C->space_cols().create_member(),
				v_chD_tmp = C->space_cols().create_member();

			if(out && print_tests >= PRINT_MORE)
				*out << "\n\n3.b.1) Testing non-transposed: alpha*C*(inv(C)*v) == alpha*v ...";
			if(out && print_tests > PRINT_MORE)
				*out << std::endl;
			llresult = true;
 			{for( int k = 1; k <= num_random_tests(); ++k ) {
				random_vector( rand_y_l, rand_y_u, v_chD.get() );
			 	if(out && print_tests >= PRINT_ALL) {
					*out
						<< "\n\n3.b.1."<<k<<") random vector " << k << " ( ||v_chD||_1 / n = " << (v_chD->norm_1() / v_chD->dim()) << " )\n";
					if(dump_all() && print_tests >= PRINT_ALL)
						*out << "\nv_chD =\n" << *v_chD;
				}
				V_InvMtV( v_xD_tmp.get(), *C, no_trans, *v_chD );
				V_StMtV( v_chD_tmp.get(), alpha, *C, no_trans, *v_xD_tmp );
				const value_type
					sum_aCICv  =         sum(*v_chD_tmp),
					sum_av     = alpha * sum(*v_chD);
				assert_print_nan_inf(sum_aCICv, "sum(alpha*C*(inv(C)*v)",true,out);
				assert_print_nan_inf(sum_av, "sum(alpha*v)",true,out);
				const value_type
					calc_err = ::fabs( ( sum_aCICv - sum_av )
									   /( ::fabs(sum_aCICv) + ::fabs(sum_av) + small_num ) );
				if(out && print_tests >= PRINT_ALL)
					*out
						<< "\nrel_err(sum(alpha*C*(inv(C)*v),sum(alpha*v)) = "
						<< "rel_err(" << sum_aCICv << "," << sum_av << ") = "
						<< calc_err << std::endl;
				if( calc_err >= warning_tol() ) {
					if(out && print_tests >= PRINT_ALL)
						*out
							<< std::endl
							<< ( calc_err >= error_tol() ? "Error" : "Warning" )
							<< ", rel_err(sum(alpha*C*(inv(C)*v)),sum(alpha*v)) = "
							<< "rel_err(" << sum_aCICv << "," << sum_av << ") = "
							<< calc_err
							<< " exceeded "
							<< ( calc_err >= error_tol() ? "error_tol" : "warning_tol" )
							<< " = "
							<< ( calc_err >= error_tol() ? error_tol() : warning_tol() )
							<< std::endl;
					if(calc_err >= error_tol()) {
						if(dump_all() && print_tests >= PRINT_ALL) {
							*out << "\nalpha = " << alpha << std::endl;
							*out << "\nv_chD =\n"                << *v_chD;
							*out << "\ninv(C)*v_chD =\n"         << *v_xD_tmp;
							*out << "\nalpha*C*inv(C)*v_chD =\n" << *v_chD_tmp;
						}
						llresult = false;
					}
				}
			}}
			if(!llresult) lresult = false;
			if( out && print_tests == PRINT_MORE )
				*out << " : " << ( llresult ? "passed" : "failed" )
					 << std::endl;

			if(out && print_tests >= PRINT_MORE)
				*out << "\n3.b.2) Testing transposed: alpha*C'*(inv(C')*v) == alpha*v ...";
			if(out && print_tests > PRINT_MORE)
				*out << std::endl;
			llresult = true;
 			{for( int k = 1; k <= num_random_tests(); ++k ) {
				random_vector( rand_y_l, rand_y_u, v_xD.get() );
			 	if(out && print_tests >= PRINT_ALL) {
					*out
						<< "\n3.b.2."<<k<<") random vector " << k << " ( ||v_xD||_1 / n = " << (v_xD->norm_1() / v_xD->dim()) << " )\n";
					if(dump_all() && print_tests >= PRINT_ALL)
						*out << "\nv_xD =\n" << *v_xD;
				}
				V_InvMtV( v_chD_tmp.get(), *C, trans, *v_xD );
				V_StMtV( v_xD_tmp.get(), alpha, *C, trans, *v_chD_tmp );
				const value_type
					sum_aCICv  =         sum(*v_xD_tmp),
					sum_av     = alpha * sum(*v_xD);
				assert_print_nan_inf(sum_aCICv, "sum(alpha*C'*(inv(C')*v)",true,out);
				assert_print_nan_inf(sum_av, "sum(alpha*v)",true,out);
				const value_type
					calc_err = ::fabs( ( sum_aCICv - sum_av )
									   /( ::fabs(sum_aCICv) + ::fabs(sum_av) + small_num ) );
				if(out && print_tests >= PRINT_ALL)
					*out
						<< "\nrel_err(sum(alpha*C'*(inv(C')*v)),sum(alpha*v)) = "
						<< "rel_err(" << sum_aCICv << "," << sum_av << ") = "
						<< calc_err << std::endl;
				if( calc_err >= warning_tol() ) {
					if(out && print_tests >= PRINT_ALL)
						*out
							<< std::endl
							<< ( calc_err >= error_tol() ? "Error" : "Warning" )
							<< ", rel_err(sum(alpha*C'*(inv(C')*v)),sum(alpha*v)) = "
							<< "rel_err(" << sum_aCICv << "," << sum_av << ") = "
							<< calc_err
							<< " exceeded "
							<< ( calc_err >= error_tol() ? "error_tol" : "warning_tol" )
							<< " = "
							<< ( calc_err >= error_tol() ? error_tol() : warning_tol() )
							<< std::endl;
					if(calc_err >= error_tol()) {
						if(dump_all() && print_tests >= PRINT_ALL) {
							*out << "\nalpha = " << alpha << std::endl;
							*out << "\nv_xD =\n"                   << *v_xD;
							*out << "\ninv(C')*v_chD =\n"          << *v_chD_tmp;
							*out << "\nalpha*C'*inv(C')*v_chD =\n" << *v_xD_tmp;
						}
						llresult = false;
					}
				}
			}}
			if(!llresult) lresult = false;
			if( out && print_tests == PRINT_MORE )
				*out << " : " << ( llresult ? "passed" : "failed" )
					 << std::endl;
		}
		
		if(D) {
			if(out && print_tests >= PRINT_MORE)
				*out
					<< "\n3.c) Check consistency of:"
					<< "\n    alpha * op(-inv(C) * N) * v == alpha * op(D) * v"
					<< "\nfor random vectors v ...";
			
		{
				VectorSpace::vec_mut_ptr_t
					v_xD      = C->space_rows().create_member(),
					v_xI      = N->space_rows().create_member(),
					v_xD_tmp  = C->space_rows().create_member(),
					v_xI_tmp  = N->space_rows().create_member(),
					v_chD_tmp = C->space_cols().create_member();

				if(out && print_tests >= PRINT_MORE)
					*out << "\n\n3.b.1) Testing non-transposed: inv(C)*(-alpha*N*v) == alpha*D*v ...";
				if(out && print_tests > PRINT_MORE)
					*out << std::endl;
				llresult = true;
 				{for( int k = 1; k <= num_random_tests(); ++k ) {
					random_vector( rand_y_l, rand_y_u, v_xI.get() );
				 	if(out && print_tests >= PRINT_ALL) {
						*out
							<< "\n\n3.b.1."<<k<<") random vector " << k << " ( ||v_xI||_1 / n = " << (v_xI->norm_1() / v_xI->dim()) << " )\n";
						if(dump_all() && print_tests >= PRINT_ALL)
							*out << "\nv_xI =\n" << *v_xI;
					}
					V_StMtV( v_chD_tmp.get(), -alpha, *N, no_trans, *v_xI );
					V_InvMtV( v_xD_tmp.get(), *C, no_trans, *v_chD_tmp );
					V_StMtV( v_xD.get(), alpha, *D, no_trans, *v_xI );
					const value_type
						sum_ICaNv  = sum(*v_xD_tmp),
						sum_aDv    = sum(*v_xD);
					assert_print_nan_inf(sum_ICaNv, "sum(inv(C)*(-alpha*N*v))",true,out);
					assert_print_nan_inf(sum_aDv, "sum(alpha*D*v)",true,out);
					const value_type
						calc_err = ::fabs( ( sum_ICaNv - sum_aDv )
										   /( ::fabs(sum_ICaNv) + ::fabs(sum_aDv) + small_num ) );
					if(out && print_tests >= PRINT_ALL)
						*out
							<< "\nrel_err(sum(inv(C)*(-alpha*N*v)),sum(alpha*D*v)) = "
							<< "rel_err(" << sum_ICaNv << "," << sum_aDv << ") = "
							<< calc_err << std::endl;
					if( calc_err >= warning_tol() ) {
						if(out && print_tests >= PRINT_ALL)
							*out
								<< std::endl
								<< ( calc_err >= error_tol() ? "Error" : "Warning" )
								<< ", rel_err(sum(inv(C)*(-alpha*N*v))),sum(alpha*D*v)) = "
								<< "rel_err(" << sum_ICaNv << "," << sum_aDv << ") = "
								<< calc_err
								<< " exceeded "
								<< ( calc_err >= error_tol() ? "error_tol" : "warning_tol" )
								<< " = "
								<< ( calc_err >= error_tol() ? error_tol() : warning_tol() )
								<< std::endl;
						if(calc_err >= error_tol()) {
							if(dump_all() && print_tests >= PRINT_ALL) {
								*out << "\nalpha = " << alpha << std::endl;
								*out << "\nv_xI =\n"                   << *v_xI;
								*out << "\n-alpha*N*v_xI =\n"          << *v_chD_tmp;
								*out << "\ninv(C)*(-alpha*N*v_xI) =\n" << *v_xD_tmp;
								*out << "\nalpha*D*v_xI =\n"           << *v_xD;
							}
							llresult = false;
						}
					}
				}}
				if(!llresult) lresult = false;
				if( out && print_tests == PRINT_MORE )
					*out << " : " << ( llresult ? "passed" : "failed" )
						 << std::endl;

				if(out && print_tests >= PRINT_MORE)
					*out << "\n3.b.1) Testing transposed: -alpha*N'*(inv(C')*v) == alpha*D'*v ...";
				if(out && print_tests > PRINT_MORE)
					*out << std::endl;
				llresult = true;
 				{for( int k = 1; k <= num_random_tests(); ++k ) {
					random_vector( rand_y_l, rand_y_u, v_xD.get() );
				 	if(out && print_tests >= PRINT_ALL) {
						*out
							<< "\n\n3.b.1."<<k<<") random vector " << k << " ( ||v_xD||_1 / n = " << (v_xD->norm_1() / v_xD->dim()) << " )\n";
						if(dump_all() && print_tests >= PRINT_ALL)
							*out << "\nv_xD =\n" << *v_xD;
					}
					V_InvMtV( v_chD_tmp.get(), *C, trans, *v_xD );
					V_StMtV( v_xI_tmp.get(), -alpha, *N, trans, *v_chD_tmp );
					V_StMtV( v_xI.get(), alpha, *D, trans, *v_xD );
					const value_type
						sum_aNTICTv  = sum(*v_xI_tmp),
						sum_aDTv     = sum(*v_xI);
					assert_print_nan_inf(sum_aNTICTv, "sum(-alpha*N'*(inv(C')*v))",true,out);
					assert_print_nan_inf(sum_aDTv, "sum(alpha*D'*v)",true,out);
					const value_type
						calc_err = ::fabs( ( sum_aNTICTv - sum_aDTv )
										   /( ::fabs(sum_aNTICTv) + ::fabs(sum_aDTv) + small_num ) );
					if(out && print_tests >= PRINT_ALL)
						*out
							<< "\nrel_err(sum(-alpha*N'*(inv(C')*v)),sum(alpha*D'*v)) = "
							<< "rel_err(" << sum_aNTICTv << "," << sum_aDTv << ") = "
							<< calc_err << std::endl;
					if( calc_err >= warning_tol() ) {
						if(out && print_tests >= PRINT_ALL)
							*out
								<< std::endl
								<< ( calc_err >= error_tol() ? "Error" : "Warning" )
								<< ", rel_err(sum(-alpha*N'*(inv(C')*v))),sum(alpha*D'*v)) = "
								<< "rel_err(" << sum_aNTICTv << "," << sum_aDTv << ") = "
								<< calc_err
								<< " exceeded "
								<< ( calc_err >= error_tol() ? "error_tol" : "warning_tol" )
								<< " = "
								<< ( calc_err >= error_tol() ? error_tol() : warning_tol() )
								<< std::endl;
						if(calc_err >= error_tol()) {
							if(dump_all() && print_tests >= PRINT_ALL) {
								*out << "\nalpha = " << alpha << std::endl;
								*out << "\nv_xD =\n"                      << *v_xD;
								*out << "\ninv(C')**v_xD =\n"             << *v_chD_tmp;
								*out << "\n-alpha*N'*(inv(C')**v_xD) =\n" << *v_xI_tmp;
								*out << "\nalpha*D*'v_xD =\n"             << *v_xI;
							}
							llresult = false;
						}
					}
				}}
				if(!llresult) lresult = false;
				if( out && print_tests == PRINT_MORE )
					*out << " : " << ( llresult ? "passed" : "failed" )
						 << std::endl;

			}
		}
		
		if( GcUP ) {
				assert(0); // ToDo: Validate GcUP and the related matrices
		}

		if( GhUP ) {
				assert(0); // ToDo: Validate GhUP and the related matrices
		}

		if(!lresult) success = false;
		if( out && print_tests == PRINT_BASIC )
			*out << " : " << ( lresult ? "passed" : "failed" );
	}

	if(out && print_tests != PRINT_NONE ) {
		if(success)
			*out << "\nCongradulations! The BasisSystem object and its associated matrix objects seem to check out!\n";
		else
			*out << "\nOops! At last one of the tests did not check out!\n";
		if(print_tests >= PRINT_BASIC)
			*out << "\nEnd BasisSystemTester::test_basis_system(...)\n";
	}
	
	return success;
}

} // end namespace AbstractLinAlgPack
