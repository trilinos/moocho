// ///////////////////////////////////////////////////////////
// MatrixWithOpNonsingularTester.cpp
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

#include "AbstractLinAlgPack/include/MatrixWithOpNonsingularTester.h"
#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"
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

MatrixWithOpNonsingularTester::MatrixWithOpNonsingularTester(
	ETestLevel       test_level
	,EPrintTestLevel print_tests
	,bool            dump_all
	,bool            throw_exception
	,size_type       num_random_tests
	,value_type      warning_tol
	,value_type      error_tol
	)
	:test_level_(test_level)
	,print_tests_(print_tests)
	,dump_all_(dump_all)
	,throw_exception_(throw_exception)
	,num_random_tests_(num_random_tests)
	,warning_tol_(warning_tol)
	,error_tol_(error_tol)
{}

bool MatrixWithOpNonsingularTester::test_matrix(
	const MatrixWithOpNonsingular   &M
	,const char                     M_name[]
	,std::ostream                   *out
	)
{
	namespace rcp = ReferenceCountingPack;
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using BLAS_Cpp::left;
	using BLAS_Cpp::right;
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
	
	bool success = true, result, lresult;
	const value_type
		rand_y_l  = -1.0,
		rand_y_u  = 1.0,
		small_num = ::pow(std::numeric_limits<value_type>::epsilon(),0.25),
		alpha     = 2.0;
	
	//
	// Perform the tests
	//

	if(out && print_tests() >= PRINT_BASIC)
		*out
			<< "\nCheck: alpha*op(op("<<M_name<<")*op("<<M_name<<"))*v == alpha*v ...";
	if(out && print_tests() > PRINT_BASIC)
		*out << std::endl;

	VectorSpace::vec_mut_ptr_t
		v_c1 = M.space_cols().create_member(),
		v_c2 = M.space_cols().create_member(),
		v_r1 = M.space_rows().create_member(),
		v_r2 = M.space_rows().create_member();

	// Side of the matrix inverse	
	const BLAS_Cpp::Side    a_side[2]  = { BLAS_Cpp::left,     BLAS_Cpp::right };
	// If the matrices are transposed or not
	const BLAS_Cpp::Transp  a_trans[2] = { BLAS_Cpp::no_trans, BLAS_Cpp::trans };

	for( int side_i = 0; side_i < 2; ++side_i ) {
		for( int trans_i = 0; trans_i < 2; ++trans_i ) {
			const BLAS_Cpp::Side    t_side  = a_side[side_i];
			const BLAS_Cpp::Transp  t_trans = a_trans[trans_i];
			if(out && print_tests() >= PRINT_MORE)
				*out
					<< "\n" << side_i+1 << "." << trans_i+1 << ") "
					<< "Check: (t2 = "<<(t_side==left?"inv(":"alpha * ")<< M_name<<(t_trans==trans?"\'":"")<<(t_side==left?")":"")
					<< " * (t1 = "<<(t_side==right?"inv(":"alpha * ")<<M_name<<(t_trans==trans?"\'":"")<<(t_side==right?")":"")
					<< " * v)) == alpha * v ...";
			if(out && print_tests() > PRINT_MORE)
				*out << std::endl;
			result = true;
			VectorWithOpMutable
				*v  = NULL,
				*t1 = NULL,
				*t2 = NULL;
			if( (t_side == left && t_trans == no_trans) || (t_side == right && t_trans == trans) ) {
				// (inv(R)*R*v or R'*inv(R')*v
				v  = v_r1.get();
				t1 = v_c1.get();
				t2 = v_r2.get();
			}
			else {
				// (inv(R')*R'*v or R*inv(R)*v
				v  = v_c1.get();
				t1 = v_r1.get();
				t2 = v_c2.get();
			}
			for( int k = 1; k <= num_random_tests(); ++k ) {
				lresult = true;
				random_vector( rand_y_l, rand_y_u, v );
					if(out && print_tests() >= PRINT_ALL) {
					*out
						<< "\n"<<side_i+1<<"."<<trans_i+1<<"."<<k<<") random vector " << k
						<< " ( ||v||_1 / n = " << (v->norm_1() / v->dim()) << " )\n";
					if(dump_all() && print_tests() >= PRINT_ALL)
						*out << "\nv =\n" << *v;
				}
				// t1
				if( t_side == right ) {
					// t1 = inv(op(M))*v
					V_InvMtV( t1, M, t_trans, *v );
				}
				else {
					// t1 = alpha*op(M)*v
					V_StMtV( t1, alpha, M, t_trans, *v );
				}
				// t2
				if( t_side == left ) {
					// t2 = inv(op(M))*t1
					V_InvMtV( t2, M, t_trans, *t1 );
				}
				else {
					// t2 = alpha*op(M)*t1
					V_StMtV( t2, alpha, M, t_trans, *t1 );
				}
				const value_type
					sum_t2  = sum(*t2),
					sum_av  = alpha*sum(*v);
				assert_print_nan_inf(sum_t2, "sum(t2)",true,out);
				assert_print_nan_inf(sum_av, "sum(alpha*t1)",true,out);
				const value_type
					calc_err = ::fabs( ( sum_av - sum_t2 )
									   /( ::fabs(sum_av) + ::fabs(sum_t2) + small_num ) );
				if(out && print_tests() >= PRINT_ALL)
					*out
						<< "\nrel_err(sum(alpha*v),sum(t2)) = "
						<< "rel_err(" << sum_av << "," << sum_t2 << ") = "
						<< calc_err << std::endl;
				if( calc_err >= warning_tol() ) {
					if(out && print_tests() >= PRINT_ALL)
						*out
							<< std::endl
							<< ( calc_err >= error_tol() ? "Error" : "Warning" )
							<< ", rel_err(sum(alpha*v),sum(t2)) = "
							<< "rel_err(" << sum_av << "," << sum_t2 << ") = "
							<< calc_err
							<< " exceeded "
							<< ( calc_err >= error_tol() ? "error_tol" : "warning_tol" )
							<< " = "
							<< ( calc_err >= error_tol() ? error_tol() : warning_tol() )
							<< std::endl;
					if(calc_err >= error_tol()) {
						if(dump_all() && print_tests() >= PRINT_ALL) {
							*out << "\nalpha = " << alpha << std::endl;
							*out << "\nv =\n"    << *v;
							*out << "\nt1 =\n"   << *t2;
							*out << "\nt2 =\n"   << *t2;
						}
						lresult = false;
					}
				}
				if(!lresult) result = false;
			}
			if(!result) success = false;
			if( out && print_tests() == PRINT_MORE )
				*out << " : " << ( result ? "passed" : "failed" )
					 << std::endl;
		}
	}

	if( out && print_tests() == PRINT_BASIC )
		*out << " : " << ( success ? "passed" : "failed" );

	return success;
}

} // end namespace AbstractLinAlgPack
