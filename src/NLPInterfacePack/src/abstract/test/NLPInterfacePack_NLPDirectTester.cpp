// ///////////////////////////////////////////////////////////
// NLPDirectTester.cpp
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
#include <math.h>

#include <ostream>
#include <iomanip>
#include <sstream>

#include "NLPDirectTester.hpp"
#include "NLPInterfacePack/src/abstract/interfaces/NLPDirect.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorMutable.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorOut.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorSpace.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorStdOps.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/LinAlgOpPack.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/assert_print_nan_inf.hpp"
#include "Range1D.hpp"
#include "update_success.hpp"
#include "ThrowException.hpp"

namespace {
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
} // end namespace

namespace NLPInterfacePack {

NLPDirectTester::NLPDirectTester(
	const calc_fd_prod_ptr_t  &calc_fd_prod
	,ETestingMethod           Gf_testing_method
	,ETestingMethod           Gc_testing_method
	,value_type               Gf_warning_tol
	,value_type               Gf_error_tol
	,value_type               Gc_warning_tol
	,value_type               Gc_error_tol
	,value_type               Gh_warning_tol
	,value_type               Gh_error_tol
	,size_type                num_fd_directions
	)
	:calc_fd_prod_(calc_fd_prod)
	,Gf_testing_method_(Gf_testing_method)
	,Gc_testing_method_(Gc_testing_method)
	,Gf_warning_tol_(Gf_warning_tol)
	,Gf_error_tol_(Gf_error_tol)
	,Gc_warning_tol_(Gc_warning_tol)
	,Gc_error_tol_(Gc_error_tol)
	,Gh_warning_tol_(Gh_warning_tol)
	,Gh_error_tol_(Gh_error_tol)
	,num_fd_directions_(num_fd_directions)
{}

bool NLPDirectTester::finite_diff_check(
	NLPDirect     *nlp
	,const Vector     &xo
	,const Vector     *xl
	,const Vector     *xu
	,const Vector     *c
	,const Vector     *h
	,const Vector     *Gf
	,const Vector     *py
	,const Vector     *rGf
	,const MatrixOp     *GcU
	,const MatrixOp     *Gh
	,const MatrixOp     *D
	,const MatrixOp     *Uz
	,const MatrixOp     *Vz
	,bool                   print_all_warnings
	,std::ostream           *out
	) const
{

	using std::setw;
	using std::endl;
	using std::right;

	using AbstractLinAlgPack::sum;
	using AbstractLinAlgPack::dot;
	using AbstractLinAlgPack::Vp_StV;
	using AbstractLinAlgPack::random_vector;
	using AbstractLinAlgPack::assert_print_nan_inf;
	using LinAlgOpPack::V_StV;
	using LinAlgOpPack::V_StMtV;
	using LinAlgOpPack::Vp_MtV;
	using LinAlgOpPack::M_StM;
	using LinAlgOpPack::M_StMtM;

	typedef VectorSpace::vec_mut_ptr_t  vec_mut_ptr_t;

//	using AbstractLinAlgPack::TestingPack::CompareDenseVectors;
//	using AbstractLinAlgPack::TestingPack::CompareDenseSparseMatrices;

	using TestingHelperPack::update_success;

	bool success = true, result, preformed_fd;
	if(out) {
		*out << std::boolalpha
			 << std::endl
			 << "*********************************************************\n"
			 << "*** NLPDirectTester::finite_diff_check(...) ***\n"
			 << "*********************************************************\n";
	}

	const size_type
		n  = nlp->n(),
		m  = nlp->m(),
		mI = nlp->mI();
	const Range1D
		var_dep      = nlp->var_dep(),
		var_indep    = nlp->var_indep(),
		con_decomp   = nlp->con_decomp(),
		con_undecomp = nlp->con_undecomp();
	NLP::vec_space_ptr_t
		space_x = nlp->space_x(),
		space_c = nlp->space_c(),
		space_h = nlp->space_h();

	// //////////////////////////////////////////////
	// Validate the input

	THROW_EXCEPTION(
		py && !c, std::invalid_argument
		,"NLPDirectTester::finite_diff_check(...) : "
		"Error, if py!=NULL then c!=NULL must also be true!" );

	const CalcFiniteDiffProd
		&fd_deriv_prod = this->calc_fd_prod();

	const value_type
		rand_v_l = -1.0, rand_v_u = 1.0,
		small_num = ::sqrt(std::numeric_limits<value_type>::epsilon());

	try {

	// ///////////////////////////////////////////////
	// (1) Check Gf

	if( Gf ) {
		switch( Gf_testing_method() ) {
			case FD_COMPUTE_ALL: {
				// Compute FDGf outright
				assert(0); // ToDo: update above!
				break;
			}
			case FD_DIRECTIONAL: {
				// Compute FDGF'*y using random y's
				if(out)
					*out
						<< "\nComparing products Gf'*y with finite difference values FDGf'*y for "
						<< "random y's ...\n";
				vec_mut_ptr_t y = space_x->create_member();
				value_type max_warning_viol = 0.0;
				int num_warning_viol = 0;
				for( int direc_i = 1; direc_i <= num_fd_directions(); ++direc_i ) {
					random_vector( rand_v_l, rand_v_u, y.get() );
					value_type
						Gf_y = dot( *Gf, *y ),
						FDGf_y;
					preformed_fd = fd_deriv_prod.calc_deriv_product(
						xo,xl,xu
						,*y,NULL,NULL,NULL,true,nlp
						,&FDGf_y,NULL,NULL,out
						);
					if( !preformed_fd )
						goto FD_NOT_PREFORMED;
					assert_print_nan_inf(FDGf_y, "FDGf'*y",true,out);
					const value_type
						calc_err = ::fabs( ( Gf_y - FDGf_y )/( ::fabs(Gf_y) + ::fabs(FDGf_y) + small_num ) );
					if( calc_err >= Gf_warning_tol() ) {
						max_warning_viol = my_max( max_warning_viol, calc_err );
						++num_warning_viol;
					}
					if( calc_err >= Gf_error_tol() ) {
						if(out)
							*out
								<< "\nError, rel_err(Gf'*y,FDGf'*y) = "
								<< "rel_err(" << Gf_y << "," << FDGf_y << ") = "
								<< calc_err << endl
								<< "exceeded Gf_error_tol = " << Gf_error_tol() << endl
								<< "Stoping the tests!\n";
						return false;
					}
				}
				if(out && num_warning_viol)
					*out
						<< "\nThere were " << num_warning_viol << " warning tolerance "
						<< "violations out of num_fd_directions = " << num_fd_directions()
						<< " computations of FDGf'*y\n"
						<< "and the maximum violation was " << max_warning_viol
						<< " > Gf_waring_tol = " << Gf_warning_tol() << endl;
				break;
			}
			default:
				assert(0); // Invalid value
		}
	}

	// /////////////////////////////////////////////
	// (2) Check py = -inv(C)*c
	//
	// We want to check; 
	// 
	//  FDC * (inv(C)*c) \approx c
  	//       \_________/
  	//         -py
  	//
  	// We can compute this as:
  	//           
  	// FDC * py = [ FDC, FDN ] * [ -py ; 0 ]
  	//            \__________/
  	//                FDA'
  	// 
  	// t1 =  [ -py ; 0 ]
  	// 
  	// t2 = FDA'*t1
  	// 
  	// Compare t2 \approx c
  	// 
	if(py) {
		if(out)
			*out
				<< "\nComparing c with finite difference product FDA'*[ -py; 0 ] = -FDC*py ...\n";
	  	// t1 =  [ -py ; 0 ]
		VectorSpace::vec_mut_ptr_t
			t1 = space_x->create_member();
		V_StV( t1->sub_view(var_dep).get(), -1.0, *py );
		*t1->sub_view(var_indep) = 0.0;
		// t2 = FDA'*t1
		VectorSpace::vec_mut_ptr_t
			t2 = nlp->space_c()->create_member();
		preformed_fd = fd_deriv_prod.calc_deriv_product(
			xo,xl,xu
			,*t1,NULL,NULL,NULL,true,nlp,NULL,t2.get(),NULL,out
			);
		if( !preformed_fd )
			goto FD_NOT_PREFORMED;
		const value_type
			sum_c  = sum(*c),
			sum_t2 = sum(*t2);
		assert_print_nan_inf(sum_t2, "sum(-FDC*py)",true,out);
		const value_type
			calc_err = ::fabs( ( sum_c - sum_t2 )/( ::fabs(sum_c) + ::fabs(sum_t2) + small_num ) );
		if( calc_err >= Gc_error_tol() ) {
			if(out)
				*out
					<< "\nError, rel_err(sum(c),sum(-FDC*py)) = "
					<< "rel_err(" << sum_c << "," << sum_t2 << ") = "
					<< calc_err << endl
					<< "exceeded Gc_error_tol = " << Gc_error_tol() << endl
					<< "Stoping the tests!\n";
			if(print_all_warnings)
				*out << "\nt1 = [ -py; 0 ] =\n" << *t1
					 << "\nt2 = FDA'*t1 = -FDC*py =\n"   << *t2;
			update_success( false, &success );
			return false;
		}
		if( calc_err >= Gc_warning_tol() ) {
			if(out)
				*out
					<< "\nWarning, rel_err(sum(c),sum(-FDC*py)) = "
					<< "rel_err(" << sum_c << "," << sum_t2 << ") = "
					<< calc_err << endl
					<< "exceeded Gc_warning_tol = " << Gc_warning_tol() << endl;
		}
	}

	// /////////////////////////////////////////////
	// (3) Check D = -inv(C)*N

	if(D) {
		switch( Gc_testing_method() ) {
			case FD_COMPUTE_ALL: {
				//
				// Compute FDN outright and check
				// -FDC * D \aprox FDN
				// 
				// We want to compute:
				// 
				// FDC * -D = [ FDC, FDN ] * [ -D; 0 ]
			  	//            \__________/
  				//                FDA'
				// 
				// To compute the above we perform:
				// 
				// T = FDA' * [ -D; 0 ] (one column at a time)
				// 
				// Compare T \approx FDN
				//
/*
				// FDN
				DMatrix FDN(m,n-m);
				fd_deriv_computer.calc_deriv( xo, xl, xu
					, Range1D(m+1,n), nlp, NULL
					, &FDN() ,BLAS_Cpp::trans, out );

				// T = FDA' * [ -D; 0 ] (one column at a time)
				DMatrix T(m,n-m);
				DVector t(n);
				t(m+1,n) = 0.0;
				for( int s = 1; s <= n-m; ++s ) {
					// t = [ -D(:,s); 0 ]
					V_StV( &t(1,m), -1.0, D->col(s) );
					// T(:,s) =  FDA' * t
					fd_deriv_prod.calc_deriv_product(
						xo,xl,xu,t(),nlp,NULL,&T.col(s),out);
				}				

				// Compare T \approx FDN
				if(out)
					*out
						<< "\nChecking the computed D = -inv(C)*N\n"
						<< "where D(i,j) = (-FDC*D)(i,j), dM(i,j) = FDN(i,j) ...\n";
				result = comp_M.comp(
					T(), FDN, BLAS_Cpp::no_trans
					, CompareDenseSparseMatrices::FULL_MATRIX
					, CompareDenseSparseMatrices::REL_ERR_BY_COL
					, Gc_warning_tol(), Gc_error_tol()
					, print_all_warnings, out );
				update_success( result, &success );
				if(!result) return false;
*/
				assert(0); // Todo: Implement above!
				break;
			}
			case FD_DIRECTIONAL: {
				//
				// Compute -FDC * D * v \aprox FDN * v
				// for random v's
				//
				// We will compute this as:
				// 
				// t1 = [ 0; y ] <: R^(n)
				// 
				// t2 = FDA' * t1  (  FDN * y ) <: R^(m)
				//
				// t1 = [ -D * y ; 0 ]  <: R^(n)
				// 
				// t3 = FDA' * t1  ( -FDC * D * y ) <: R^(m)
				// 
				// Compare t2 \approx t3
				//
				if(out)
					*out
						<< "\nComparing finite difference products -FDC*D*y with FDN*y for "
							"random y's ...\n";
				VectorSpace::vec_mut_ptr_t
					y  = space_x->sub_space(var_indep)->create_member(),
					t1 = space_x->create_member(),
					t2 = space_c->create_member(),
					t3 = space_c->create_member();
				value_type max_warning_viol = 0.0;
				int num_warning_viol = 0;
				for( int direc_i = 1; direc_i <= num_fd_directions(); ++direc_i ) {
					random_vector( rand_v_l, rand_v_u, y.get() );
					// t1 = [ 0; y ] <: R^(n)
					*t1->sub_view(var_dep)   = 0.0;
					*t1->sub_view(var_indep) = *y;
					// t2 = FDA' * t1  (  FDN * y ) <: R^(m)
					preformed_fd = fd_deriv_prod.calc_deriv_product(
						xo,xl,xu
						,*t1,NULL,NULL,NULL,true,nlp,NULL,t2.get(),NULL,out
						);
					if( !preformed_fd )
						goto FD_NOT_PREFORMED;
					// t1 = [ -D * y ; 0 ]  <: R^(n)
					V_StMtV( t1->sub_view(var_dep).get(), -1.0, *D, BLAS_Cpp::no_trans, *y );
					*t1->sub_view(var_indep) = 0.0;
					// t3 = FDA' * t1  ( -FDC * D * y ) <: R^(m)
					preformed_fd = fd_deriv_prod.calc_deriv_product(
						xo,xl,xu
						,*t1,NULL,NULL,NULL,true,nlp,NULL,t3.get(),NULL,out
						);
					// Compare t2 \approx t3
					const value_type
						sum_t2 = sum(*t2),
						sum_t3 = sum(*t3);
					const value_type
						calc_err = ::fabs( ( sum_t2 - sum_t3 )/( ::fabs(sum_t2) + ::fabs(sum_t3) + small_num ) );
					if( calc_err >= Gc_warning_tol() ) {
						max_warning_viol = my_max( max_warning_viol, calc_err );
						++num_warning_viol;
					}
					if( calc_err >= Gc_error_tol() ) {
						if(out)
							*out
								<< "\nError, rel_err(sum(-FDC*D*y),sum(FDN*y)) = "
								<< "rel_err(" << sum_t3 << "," << sum_t2 << ") = "
								<< calc_err << endl
								<< "exceeded Gc_error_tol = " << Gc_error_tol() << endl
								<< "Stoping the tests!\n";
						if(print_all_warnings)
							*out << "\ny =\n" << *y
								 << "\nt1 = [ -D*y; 0 ] =\n" << *t1
								 << "\nt2 =  FDA' * [ 0; y ] = FDN * y =\n" << *t2
								 << "\nt3 =  FDA' * t1 = -FDC * D * y =\n" << *t3;
						update_success( false, &success );
						return false;
					}
				}
				if(out && num_warning_viol)
					*out
						<< "\nThere were " << num_warning_viol << " warning tolerance "
						<< "violations out of num_fd_directions = " << num_fd_directions()
						<< " computations of sum(FDC*D*y) and sum(FDN*y)\n"
						<< "and the maximum relative iolation was " << max_warning_viol
						<< " > Gc_waring_tol = " << Gc_warning_tol() << endl;
				break;
			}
			default:
				assert(0);
		}
	}

	// ///////////////////////////////////////////////
	// (4) Check rGf = Gf(var_indep) + D'*Gf(var_dep)

	if(rGf) {
		if( Gf && D ) {
			if(out)
				*out
					<< "\nComparing rGf_tmp = Gf(var_indep) - D'*Gf(var_dep) with rGf ...\n";
			VectorSpace::vec_mut_ptr_t
				rGf_tmp = space_x->sub_space(var_indep)->create_member();
			*rGf_tmp = *Gf->sub_view(var_indep);
			Vp_MtV( rGf_tmp.get(), *D, BLAS_Cpp::trans, *Gf->sub_view(var_dep) );
			const value_type
				sum_rGf_tmp  = sum(*rGf_tmp),
				sum_rGf      = sum(*rGf);
			const value_type
				calc_err = ::fabs( ( sum_rGf_tmp - sum_rGf )/( ::fabs(sum_rGf_tmp) + ::fabs(sum_rGf) + small_num ) );
			if( calc_err >= Gc_error_tol() ) {
				if(out)
					*out
						<< "\nError, rel_err(sum(rGf_tmp),sum(rGf)) = "
						<< "rel_err(" << sum_rGf_tmp << "," << sum_rGf << ") = "
						<< calc_err << endl
						<< "exceeded Gc_error_tol = " << Gc_error_tol() << endl
						<< "Stoping the tests!\n";
				if(print_all_warnings)
					*out << "\nrGf_tmp =\n" << *rGf_tmp
						 << "\nrGf =\n"   << *rGf;
				update_success( false, &success );
				return false;
			}
			if( calc_err >= Gc_warning_tol() ) {
				if(out)
					*out
						<< "\nWarning, rel_err(sum(rGf_tmp),sum(rGf)) = "
						<< "rel_err(" << sum_rGf_tmp << "," << sum_rGf << ") = "
						<< calc_err << endl
						<< "exceeded Gc_warning_tol = " << Gc_warning_tol() << endl;
			}
		}
		else {
			assert(0); //  ToDo: Must validate rGf without Gf and D!
		}
	}

	// ///////////////////////////////////////////////////
	// (5) Check GcU, and/or Uz (for undecomposed equalities)

	if(GcU || Uz) {
		assert(0); // ToDo: Implement!
	}
	
	// ///////////////////////////////////////////////////
	// (6) Check Gh, and/or Vz (for general inequalities)

	if(Gh || Vz) {
		assert(0); // ToDo: Implement!
	}

FD_NOT_PREFORMED:

	if(!preformed_fd) {
		if(out)
			*out
				<< "\nError, the finite difference computation was not preformed due to cramped bounds\n"
				<< "Finite difference test failed!\n" << endl;
		return false;
	}

	} // end try
	catch( const AbstractLinAlgPack::NaNInfException& except ) {
		if(out)
			*out
				<< "Error, found a NaN or Inf.  Stoping tests\n";
		success = false;
	}

	if( out && success )
		*out
			<< "\nCongradulations, all the finite difference errors where within the\n"
				"specified error tolerances!\n";

	return success;

}

}	// end namespace NLPInterfacePack
