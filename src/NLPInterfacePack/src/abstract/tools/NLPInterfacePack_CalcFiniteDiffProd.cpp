// //////////////////////////////////////////////////////////////////////////////
// CalcFiniteDiffProd.cpp
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

//#include <iostream> // debugging only!
#include <typeinfo>
#include <iomanip>
#include <sstream>
#include <limits>

#include "NLPInterfacePack/src/abstract/tools/CalcFiniteDiffProd.hpp"
#include "NLPInterfacePack/src/abstract/interfaces/NLP.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorSpace.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorMutable.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/LinAlgOpPack.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/VectorAuxiliaryOps.hpp"
#include "Range1D.hpp"
#include "Teuchos_TestForException.hpp"

namespace NLPInterfacePack {

CalcFiniteDiffProd::CalcFiniteDiffProd(
	EFDMethodOrder              fd_method_order
	,EFDStepSelect              fd_step_select
	,value_type                 fd_step_size
	,value_type                 fd_step_size_min
	,value_type                 fd_step_size_f
	,value_type                 fd_step_size_c
	)
	:fd_method_order_(fd_method_order)
	,fd_step_select_(fd_step_select)
	,fd_step_size_(fd_step_size)
	,fd_step_size_min_(fd_step_size_min)
	,fd_step_size_f_(fd_step_size_f)
	,fd_step_size_c_(fd_step_size_c)
{}

bool CalcFiniteDiffProd::calc_deriv_product(
	const Vector       &xo
	,const Vector      *xl
	,const Vector      *xu
	,const Vector      &v
	,const value_type  *fo
	,const Vector      *co
	,bool              check_nan_inf
	,NLP               *nlp
	,value_type        *Gf_prod
	,VectorMutable     *Gc_prod
	,std::ostream      *out
	) const
{

	using std::setw;
	using std::endl;
	using std::right;

	using BLAS_Cpp::rows;
	using BLAS_Cpp::cols;
	
	namespace rcp = MemMngPack;
	typedef VectorSpace::vec_mut_ptr_t  vec_mut_ptr_t;
	using AbstractLinAlgPack::Vt_S;
	using AbstractLinAlgPack::Vp_StV;
	using AbstractLinAlgPack::max_near_feas_step;
	using AbstractLinAlgPack::assert_print_nan_inf;
	using LinAlgOpPack::V_StV;

	//
	// The gradient of the contraints is defined as the matrix
	// Gc as:
	//
	// Gc = [ Gc1, Gc2, ..., Gcm ]
	//
	//		[	dc1/dx(1)		dc2/dx(1)	...		dcm/dx(1)	]
	//		[	dc1/dx(2)		dc2/dx(2)	...		dcm/dx(2)	]
	// Gc= 	[	.				.			...		.			]
	//		[	dc1/dx(n)		dc2/dx(n)	...		dcm/dx(n)	]
	//
	//		[	(dc/dx(1))'	]
	//		[	(dc/dx(2))'	]
	// Gc =	[	.			]
	//		[	(dc/dx(n))'	]
	//
	// The gradient of the objective function is defined as the
	// vector Gf as:
	//
	//		[	(df/dx(1))'	]
	//		[	(df/dx(2))'	]
	// Gf =	[	.			]
	//		[	(df/dx(n))'	]
	//
	// To illustrate the theory behind this implementation consider
	// the generic multi-variable function g(x) <: R^n -> R.  No let's
	// consider we have the base point xo and the vector v to
	// perturb g(x) along.  First form the function g(xo+a*v) and then
	// let's compute dg/da at a = 0:
	// 
	// (1) d(g(xo+a*v))/d(a) at a = 0
	//         = sum( dg/dx(i) * dx(i)/da, i = 1...n)
	//         = sum( dg/dx(i) * v(i), i = 1...n)
	//         = Gf'*v
	//
	// Now we can approximate (1) using central differences as:
	// 
	// (2) d(g(xo+a*v))/d(a) at a = 0
	//          \approx ( g(xo+h*v) - g(xo+h*v) ) / (2*h)
	//
	// If we equate (1) and (2) we have the approximation:
	// 
	// (3) Gg' * v \approx ( g(xo+h*v) - g(xo+h*v) ) / (2*h)
	// 
	// It is clear how this applies to computing Gf'*v and Gc'*v.
	// 

	const size_type
		n = nlp->n(),
		m = nlp->m();

	const value_type
		max_bnd_viol = nlp->max_var_bounds_viol();

	// /////////////////////////////////////////
	// Validate the input

	TEST_FOR_EXCEPTION(
		m==0 && Gc_prod, std::invalid_argument
		,"CalcFiniteDiffProd::calc_deriv(...) : "
		"Error, if nlp->m() == 0, then Gc_prod must equal NULL" );
	TEST_FOR_EXCEPTION(
		Gc_prod && !Gc_prod->space().is_compatible(*nlp->space_c())
		,std::invalid_argument
		,"CalcFiniteDiffProd::calc_deriv(...) : "
		"Error, Gc_prod (type \' "<<typeid(*Gc_prod).name()<<"\' "
		"is not compatible with the NLP" );
	TEST_FOR_EXCEPTION(
		(xl && !xu) || (!xl && xu), std::invalid_argument
		,"CalcFiniteDiffProd::calc_deriv(...) : "
		"Error, both xl = "<<xl<<" and xu = "<<xu
		<<" must be NULL or not NULL" );

	assert_print_nan_inf(xo,"xo",true,out); 

	// ////////////////////////
	// Find the step size

	//
	// Get defaults for the optimal step sizes
	//

	const value_type
		sqrt_epsilon = ::pow(std::numeric_limits<value_type>::epsilon(),1.0/2.0),
		u_optimal_1  = sqrt_epsilon,
		u_optimal_2  = ::pow(sqrt_epsilon,1.0/2.0),
		u_optimal_4  = ::pow(sqrt_epsilon,1.0/4.0),
		xo_norm_inf  = xo.norm_inf();

	value_type
		uh_opt = 0.0;
	switch(this->fd_method_order()) {
		case FD_ORDER_ONE:
			uh_opt = u_optimal_1 * ( fd_step_select() == FD_STEP_ABSOLUTE ? 1.0 : xo_norm_inf + 1.0 );
			break;
		case FD_ORDER_TWO:
		case FD_ORDER_TWO_CENTRAL:
		case FD_ORDER_TWO_AUTO:
			uh_opt = u_optimal_2 * ( fd_step_select() == FD_STEP_ABSOLUTE ? 1.0 : xo_norm_inf + 1.0 );
			break;
		case FD_ORDER_FOUR:
		case FD_ORDER_FOUR_CENTRAL:
		case FD_ORDER_FOUR_AUTO:
			uh_opt = u_optimal_4 * ( fd_step_select() == FD_STEP_ABSOLUTE ? 1.0 : xo_norm_inf + 1.0 );
			break;
		default:
			assert(0); // Should not get here!
	}

	//
	// Set the step sizes used.
	//

	value_type
		uh      = this->fd_step_size(),
		uh_f    = this->fd_step_size_f(),
		uh_c    = this->fd_step_size_c(),
		uh_min  = this->fd_step_size_min();

	// uh
	if( uh < 0 )
		uh = uh_opt;
	else if(fd_step_select() == FD_STEP_RELATIVE)
		uh *= (xo_norm_inf + 1.0);
	// uh_f
	if( uh_f < 0 )
		uh_f = uh;
	else if(fd_step_select() == FD_STEP_RELATIVE)
		uh_f *= (xo_norm_inf + 1.0);
	// uh_c
	if( uh_c < 0 )
		uh_c = uh;
	else if(fd_step_select() == FD_STEP_RELATIVE)
		uh_c *= (xo_norm_inf + 1.0);
	// uh_h

	//
 	// Determine the maximum step size that can be used and
	// still stay in the relaxed bounds.
	//
	// ToDo: Consider cramped bounds, one sided differences!
	//

	value_type max_u_feas = std::numeric_limits<value_type>::max();
	if( xl ) {
		std::pair<value_type,value_type>
			u_pn
			= max_near_feas_step(
				xo
				,v
				,*xl
				,*xu
				,max_bnd_viol
				);
		if( u_pn.first < -u_pn.second )
			max_u_feas = u_pn.first;
		else
			max_u_feas = u_pn.second;
		const value_type abs_max_u_feas = ::fabs(max_u_feas);
		if( abs_max_u_feas < uh ) {
			if( abs_max_u_feas < uh_min ) {
				if(out)
					*out
						<< "\nCalcFiniteDiffProd::calc_deriv_product(...) : "
						<< "Warning, the size of the maximum finite difference step length\n"
						<< "that does not violate the relaxed variable bounds uh = "
						<< max_u_feas << " is less than the mimimum allowable step length\n"
						<< "uh_min = " << uh_min << " and the finite difference "
						<< "derivatives are not computed!\n";
				return false;
			}
			if(out)
				*out
					<< "\nCalcFiniteDiffProd::calc_deriv_product(...) : "
					<< "Warning, the size of the maximum finite difference step length\n"
					<< "that does not violate the relaxed variable bounds uh = "
					<< max_u_feas << " is less than the desired step length\n"
					<< "uh = " << uh << " and the finite difference "
					<< "derivatives may be much less accurate!\n";
		}
	}

	//
	// Set the actual method being used
	//
	// ToDo: Consider cramped bounds and method order.
	//
	
	EFDMethodOrder  fd_method_order = this->fd_method_order();
	switch(fd_method_order) {
		case FD_ORDER_TWO_AUTO:
			fd_method_order = FD_ORDER_TWO_CENTRAL;
			break;
		case FD_ORDER_FOUR_AUTO:
			fd_method_order = FD_ORDER_FOUR_CENTRAL;
			break;
	}

	// Compute the actual individual step size so as to stay in bounds
	const value_type
		abs_max_u_feas = ::fabs(max_u_feas);
	value_type
	   num_u_i = 0;
	switch(fd_method_order) {
		case FD_ORDER_ONE:
			num_u_i = 1.0;
			break;
		case FD_ORDER_TWO:
			num_u_i = 2.0;
			break;
		case FD_ORDER_TWO_CENTRAL:
			num_u_i = 1.0;
			break;
		case FD_ORDER_FOUR:
			num_u_i = 4.0;
			break;
		case FD_ORDER_FOUR_CENTRAL:
			num_u_i = 2.0;
			break;
		default:
			assert(0); // Should not get here!
	}

	uh   = ( abs_max_u_feas/num_u_i < uh   ? max_u_feas/num_u_i : uh   ); // This can be a negative number!
	uh_f = ( abs_max_u_feas/num_u_i < uh_f ? max_u_feas/num_u_i : uh_f ); //""
	uh_c = ( abs_max_u_feas/num_u_i < uh_c ? max_u_feas/num_u_i : uh_c ); //""

	if( uh_min < 0 ) {
		uh_min = uh / 100.0;
	}

//	std::cerr
//		<< "uh_opt = " << uh_opt << ", uh = " << uh << ", uh_f = " << uh_f
//		<< ", uh_c = " << uh_c << ", uh_h = " << uh_h << std::endl;

	//
	// Remember some stuff
	//
	
	value_type        *f_saved = NULL;
	VectorMutable     *c_saved = NULL;

	f_saved = nlp->get_f();
	if(m)  c_saved = nlp->get_c();

	int p_saved;
	if(out)
		p_saved = out->precision();

	// ///////////////////////////////////////////////
	// Compute the directional derivatives

	try {

	value_type
		f;
	vec_mut_ptr_t
		x = nlp->space_x()->create_member();
	vec_mut_ptr_t
		c = m  && Gc_prod ? nlp->space_c()->create_member() : Teuchos::null;
	
	// Set the quanitities used to compute with

	nlp->set_f(&f);
	if(m)  nlp->set_c( c.get() );

	const int dbl_p = 15;
	if(out)
		*out << std::setprecision(dbl_p);

	//
	// Compute the weighted sum of the terms
	//

	int          num_evals  = 0;
	value_type   dwgt       = 0.0;
	switch(fd_method_order) {
		case FD_ORDER_ONE: // may only need one eval if f(xo) etc is passed in
			num_evals = 2;
			dwgt      = 1.0;
			break;
		case FD_ORDER_TWO: // may only need two evals if c(xo) etc is passed in
			num_evals = 3;
			dwgt      = 2.0;
			break;
		case FD_ORDER_TWO_CENTRAL:
			num_evals = 2;
			dwgt      = 2.0;
			break;
		case FD_ORDER_FOUR:
			num_evals = 5;
			dwgt      = 12.0;
			break;
		case FD_ORDER_FOUR_CENTRAL:
			num_evals = 5;
			dwgt      = 12.0;
			break;
		default:
			assert(0); // Should not get here!
	}
	for( int eval_i = 1; eval_i <= num_evals; ++eval_i ) {
		// Set the step constant and the weighting constant
		value_type
			uh_i   = 0.0,
			wgt_i  = 0.0;
		switch(fd_method_order) {
			case FD_ORDER_ONE: {
				switch(eval_i) {
					case 1:
						uh_i  = +0.0;
						wgt_i = -1.0;
						break;
					case 2:
						uh_i  = +1.0;
						wgt_i = +1.0;
						break;
				}
				break;
			}
			case FD_ORDER_TWO: {
				switch(eval_i) {
					case 1:
						uh_i  = +0.0;
						wgt_i = -3.0;
						break;
					case 2:
						uh_i  = +1.0;
						wgt_i = +4.0;
						break;
					case 3:
						uh_i  = +2.0;
						wgt_i = -1.0;
						break;
				}
				break;
			}
			case FD_ORDER_TWO_CENTRAL: {
				switch(eval_i) {
					case 1:
						uh_i  = -1.0;
						wgt_i = -1.0;
						break;
					case 2:
						uh_i  = +1.0;
						wgt_i = +1.0;
						break;
				}
				break;
			}
			case FD_ORDER_FOUR: {
				switch(eval_i) {
					case 1:
						uh_i  = +0.0;
						wgt_i = -25.0;
						break;
					case 2:
						uh_i  = +1.0;
						wgt_i = +48.0;
						break;
					case 3:
						uh_i  = +2.0;
						wgt_i = -36.0;
						break;
					case 4:
						uh_i  = +3.0;
						wgt_i = +16.0;
						break;
					case 5:
						uh_i  = +4.0;
						wgt_i = -3.0;
						break;
				}
				break;
			}
			case FD_ORDER_FOUR_CENTRAL: {
				switch(eval_i) {
					case 1:
						uh_i  = -2.0;
						wgt_i = +1.0;
						break;
					case 2:
						uh_i  = -1.0;
						wgt_i = -8.0;
						break;
					case 3:
						uh_i  = +1.0;
						wgt_i = +8.0;
						break;
					case 4:
						uh_i  = +2.0;
						wgt_i = -1.0;
						break;
				}
				break;
			}
		}

		// Compute the weighted term and add it to the sum
		bool new_point = true;
		if(Gc_prod) {
			if( co && uh_i == 0.0 ) {
				*c = *co;
			}
			else {
				if( new_point || uh_c != uh ) {
					*x = xo; Vp_StV( x.get(), uh_i * uh_c, v ); // x = xo + uh_i*uh_c*v
				}
				nlp->calc_c(*x,new_point);
			}
			new_point = false;
			if(check_nan_inf)
				assert_print_nan_inf(*c,"c(xo+u*v)",true,out);
			if(eval_i == 1)
				V_StV( Gc_prod, wgt_i, *c );
			else
				Vp_StV( Gc_prod, wgt_i, *c );
		}
		
		if(Gf_prod) {
			if( fo && uh_i == 0.0 ) {
				f = *fo;
			}
			else {
				if( new_point || uh_f != uh ) {
					*x = xo; Vp_StV( x.get(), uh_i * uh_f, v ); // x = xo + uh_i*uh_f*v
				}
				nlp->calc_f(*x,new_point);
			}
			new_point = false;
			if(check_nan_inf)
				assert_print_nan_inf(f,"f(xo+u*v)",true,out);
			if(eval_i == 1)
				*Gf_prod = wgt_i * f;
			else
				*Gf_prod += wgt_i * f;
		}

	}

	//
	// Multiply by the scaling factor!
	//

	if(Gc_prod) {
		Vt_S( Gc_prod, 1.0 / (dwgt * uh_c) );
	}
		
	if(Gf_prod) {
		*Gf_prod *= ( 1.0 / (dwgt * uh_f) );
	}

	}	// end try
	catch(...) {
		nlp->set_f( f_saved );
		if(m)  nlp->set_c( c_saved );
		if(out)
			*out << std::setprecision(p_saved);
		throw;
	}
	
	nlp->set_f( f_saved );
	if(m)  nlp->set_c( c_saved );
	if(out)
		*out << std::setprecision(p_saved);
	
	return true;
}

}	// end namespace NLPInterfacePack
