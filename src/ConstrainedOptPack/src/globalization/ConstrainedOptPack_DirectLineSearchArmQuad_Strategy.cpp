// //////////////////////////////////////////////////////////////////////////////////
// DirectLineSearchArmQuad_Strategy.cpp
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

#include <ostream>
#include <iomanip>
#include <sstream>

#include "ConstrainedOptimizationPack/src/DirectLineSearchArmQuad_Strategy.h"
#include "ConstrainedOptimizationPack/src/MeritFuncCalc1D.h"
#include "RTOpPack/src/check_nan_inf.h"

namespace ConstrainedOptimizationPack {
inline value_type min(value_type v1, value_type v2) {
	return (v1 < v2) ? v1 : v2;
}
inline value_type max(value_type v1, value_type v2) {
	return (v1 > v2) ? v1 : v2;
}
}

ConstrainedOptimizationPack::DirectLineSearchArmQuad_Strategy::DirectLineSearchArmQuad_Strategy(
		int max_iter, value_type eta, value_type min_frac, value_type max_frac )
	: max_iter_(max_iter), eta_(eta), min_frac_(min_frac), max_frac_(max_frac)
{}

void ConstrainedOptimizationPack::DirectLineSearchArmQuad_Strategy::set_max_iter(int max_iter)
{
	max_iter_ = max_iter;
}

int ConstrainedOptimizationPack::DirectLineSearchArmQuad_Strategy::max_iter() const
{
	return max_iter_;
}

int ConstrainedOptimizationPack::DirectLineSearchArmQuad_Strategy::num_iterations() const
{
	return num_iter_;
}

bool ConstrainedOptimizationPack::DirectLineSearchArmQuad_Strategy::do_line_search(
	  const MeritFuncCalc1D& phi, value_type phi_k
	, value_type* alpha_k, value_type* phi_kp1
	, std::ostream* out )
{
	using std::setw;
	using std::endl;

	if(*alpha_k < 0.0) {
		throw std::logic_error(	"DirectLineSearchArmQuad_Strategy::do_line_search(): "
			"alpha_k can't start out less than 0.0"	);
	}

	validate_parameters();
	
	int w = 20;
	int prec = 8;
	int prec_saved;
	if(out) {
		prec_saved = out->precision();
		*out	<< std::setprecision(prec)
				<< "\nStarting Armijo Quadratic interpolation linesearch ...\n";
	}

	// Loop initialization (technically the first iteration)

	const value_type Dphi_k = phi.deriv();

	// output header
	if(out)
		*out	<< "\nDphi_k = "	<< Dphi_k
				<< "\nphi_k = "		<< phi_k << "\n\n"
				<< setw(5)			<< "itr"
				<< setw(w)			<< "alpha_k"
				<< setw(w)			<< "phi_kp1"
				<< setw(w)			<< "phi_kp1-frac_phi\n"
				<< setw(5)			<< "----"
				<< setw(w)			<< "----------------"
				<< setw(w)			<< "----------------"
				<< setw(w)			<< "----------------\n";

	// Check that this is a decent direction
	if(Dphi_k >= 0) throw DirectLineSearch_Strategy::NotDescentDirection(
		"DirectLineSearchArmQuad_Strategy::do_line_search(): "
		"The given d_k is not a descent direction for the given "
		"phi (phi.deriv() is >= 0)"			);

	// keep memory of the best value
	value_type	best_alpha = *alpha_k, best_phi = *phi_kp1;

	// Perform linesearch.
	bool success = false;
	for( num_iter_ = 0; num_iter_ < max_iter(); ++num_iter_ ) {

		// Print out this iteration.

		value_type frac_phi = phi_k + eta() * (*alpha_k) * Dphi_k;
		if(out)
			*out	<< setw(5)			<< num_iter_
					<< setw(w)			<< *alpha_k
					<< setw(w)			<< *phi_kp1
					<< setw(w)			<< ((*phi_kp1)-frac_phi)	<< endl;
		
		// Check that this is a valid number.
		if( RTOp_is_nan_inf( *phi_kp1 ) ) {
			// Cut back the step to min_frac * alpha_k
			*alpha_k = min_frac()*(*alpha_k);
			best_alpha = 0.0;
			best_phi = phi_k;
		}
		else {		

			// Armijo condition
			if( *phi_kp1 < frac_phi ) {
				// We have found our point.
				success = true;
				break;	// get out of the loop
			}

			// Select a new alpha to try:
			//   alpha_k = ( min_frac*alpha_k <= quadratic interpolation <= max_frac*alpha_k )

			// Quadratic interpolation of alpha_k that minimizes phi.
			// We know the values of phi at the initail point and alpha_k and
			// the derivate of phi w.r.t alpha at the initial point and
			// that's enough information for a quandratic interpolation.
			
			value_type alpha_quad =		( -0.5 * Dphi_k * (*alpha_k) * (*alpha_k) )
										/ ( (*phi_kp1) - phi_k - (*alpha_k) * Dphi_k );

			*alpha_k = min( max(min_frac()*(*alpha_k),alpha_quad), max_frac()*(*alpha_k) );

		}
		
		// Evaluate the point

		*phi_kp1 = phi(*alpha_k);

		// Save the best point found
		if(*phi_kp1 < best_phi) {
			best_phi = *phi_kp1;
			best_alpha = *alpha_k;
		}

	}

	// Be nice and reset the precision
	if(out) {
		out->precision(prec_saved);
	}

	if( success ) {
		return true;
	}

	// Line search failure.  Return the best point found and let the 
	// client decide what to do.
	*alpha_k = best_alpha;
	*phi_kp1 = phi(best_alpha);	// Make this the last call to phi(x)
	return false; 
}

void ConstrainedOptimizationPack::DirectLineSearchArmQuad_Strategy::print_algorithm(
	std::ostream& out, const std::string& L) const
{
	out
		<< L << "*** start line search using the Armijo cord test and quadratic interpolation of alpha\n"
		<< L << "default: max_ls_iter = " << max_iter() << std::endl
		<< L << "         eta = " << eta() << std::endl
		<< L << "         min_frac = " << min_frac() << std::endl
		<< L << "         max_frac = " << max_frac() << std::endl
		<< L << "Dphi_k = phi.deriv()\n"
		<< L << "if Dphi_k >= 0\n"
		<< L << "    throw not_descent_direction\n"
		<< L << "    end line search\n"
		<< L << "end\n"
		<< L << "best_alpha = alpha_k\n"
		<< L << "best_phi = phi_kp1\n"
		<< L << "for num_iter = 0... max_ls_iter\n"
		<< L << "    frac_phi = phi_k + eta * alpha_k * Dphi_k\n"
		<< L << "    print iteration\n"
		<< L << "    if phi_kp1 is not a valid number then\n"
		<< L << "        *** Cut back the step so the NLP's functions\n"
		<< L << "        *** will hopefully be defined.\n"
		<< L << "        alpha_k = min_frac * alpha_k\n"
		<< L << "        best_alpha = 0\n"
		<< L << "        best_phi = phi_k\n"
		<< L << "    else\n"
		<< L << "        if phi_kp1 < frac_phi then\n"
		<< L << "            *** We have found our point\n"
		<< L << "            end line search\n"
		<< L << "        end\n"
		<< L << "        *** Use a quadratic interpoation to minimize phi(alpha)\n"
		<< L << "        alpha_quad = (-0.5 * Dphi_k * alpha_k^2) / ( phi_kp1 - phi_k - alpha_k*Dphi_k )\n"
		<< L << "        alpha_k = min( max( min_frac*alpha_k, alpha_quad ), max_frac*alpha_k )\n"
		<< L << "    end\n"
		<< L << "    phi_kp1 = phi(alpha_k)\n"
		<< L << "    if phi_kp1 < best_phi\n"
		<< L << "        best_phi = phi_kp1\n"
		<< L << "        best_alpha = alpha_k\n"
		<< L << "    end\n"
		<< L << "end\n"
		<< L << "*** If you get there the line search failed.\n"
		<< L << "alpha_k = best_alpha\n"
		<< L << "phi_kp1 = phi(alpha_k)\n"
		<< L << "linesearch_failure = true\n";
}

void ConstrainedOptimizationPack::DirectLineSearchArmQuad_Strategy::validate_parameters() const
{
	if( eta() < 0.0 || 1.0 < eta() ) {
		std::ostringstream omsg;
		omsg
			<< "DirectLineSearchArmQuad_Strategy::validate_parameters() : "
			<< "Error, eta = " << eta() << " is not in the range [0, 1]";
		throw std::invalid_argument( omsg.str() );
	}
	if( min_frac() < 0.0 || 1.0 < min_frac() ) {
		std::ostringstream omsg;
		omsg
			<< "DirectLineSearchArmQuad_Strategy::validate_parameters() : "
			<< "Error, min_frac = " << min_frac() << " is not in the range [0, 1]";
		throw std::invalid_argument( omsg.str() );
	}
	if( max_frac() < 0.0 || 1.0 < max_frac() ) {
		std::ostringstream omsg;
		omsg
			<< "DirectLineSearchArmQuad_Strategy::validate_parameters() : "
			<< "Error, max_frac = " << max_frac() << " is not in the range [0, 1]";
		throw std::invalid_argument( omsg.str() );
	}
	if( max_frac() < min_frac() ) {
		std::ostringstream omsg;
		omsg
			<< "DirectLineSearchArmQuad_Strategy::validate_parameters() : "
			<< "Error, max_frac = " << max_frac()
			<< " < min_frac = " << min_frac();;
		throw std::invalid_argument( omsg.str() );
	}
}
