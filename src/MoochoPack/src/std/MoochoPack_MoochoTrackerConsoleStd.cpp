// ////////////////////////////////////////////////////////////////////////////
// rSQPTrackConsoleStd.cpp
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

#include <iomanip>

#include "ReducedSpaceSQPPack/src/std/rSQPTrackConsoleStd.hpp"
#include "ReducedSpaceSQPPack/src/rSQPState.hpp"
#include "ReducedSpaceSQPPack/src/rsqp_algo_conversion.hpp"
#include "NLPInterfacePack/src/NLPFirstOrderInfo.hpp"
#include "AbstractLinAlgPack/src/Vector.hpp"
#include "dynamic_cast_verbose.hpp"

namespace ReducedSpaceSQPPack {

using std::endl;
using std::setw;
using std::left;
using std::right;
using std::setprecision;

// Static members
int		rSQPTrackConsoleStd::w_i2_		= 2;
char	rSQPTrackConsoleStd::ul_i2_[]	= "--";
int		rSQPTrackConsoleStd::w_i4_		= 4;
char	rSQPTrackConsoleStd::ul_i4_[]	= "----";
int		rSQPTrackConsoleStd::p2_		= 1;
int		rSQPTrackConsoleStd::w_p2_		= 8;
char	rSQPTrackConsoleStd::ul_p2_[]	= "--------";
int		rSQPTrackConsoleStd::p3_		= 2;
int		rSQPTrackConsoleStd::w_p3_		= 9;
char	rSQPTrackConsoleStd::ul_p3_[]	= "---------";

rSQPTrackConsoleStd::rSQPTrackConsoleStd(
	const ostream_ptr_t&   o
	,const ostream_ptr_t&  journal_out
	)
	:rSQPTrack(journal_out)
	,o_(o)
	,printed_lines_(NUM_PRINT_LINES)
{}

void rSQPTrackConsoleStd::set_output_stream(const ostream_ptr_t& o)
{
	o_ = o;
}

const rSQPTrackConsoleStd::ostream_ptr_t&
rSQPTrackConsoleStd::get_output_stream() const
{
	return o_;
}

void rSQPTrackConsoleStd::initialize()
{
	timer_.reset();
	timer_.start();
}

void rSQPTrackConsoleStd::output_iteration(const Algorithm& p_algo) const
{
	const rSQPAlgo  &algo = rsqp_algo(p_algo);
	const rSQPState &s    = algo.rsqp_state();
	const NLP       &nlp  = algo.nlp(); 
	
	const size_type
		m = nlp.m();

	if(s.k() == 0) {
		print_top_header(s,algo);
		printed_lines_ = NUM_PRINT_LINES;
	}
	
	// Output the table's header
	if(printed_lines_ == NUM_PRINT_LINES) {
		printed_lines_ = 0;
		print_header(s,algo);
	}

	// ///////////////////////////////
	// Output a row for the iteration
	
	// Get a Quasi-Newton statistics.
	const QuasiNewtonStats	*quasi_newt_stats =
		( quasi_newton_stats_.exists_in(s) && quasi_newton_stats_(s).updated_k(0)
		  ? &quasi_newton_stats_(s).get_k(0)
		  : NULL );

	std::ostream& o = this->o();

	// k
	o << " " << right << setw(w_i4_) << s.k();
	// f
	if( s.f().updated_k(0) )
		o << " " << setprecision(p3_) << right << setw(w_p3_) << s.f().get_k(0);
	else
		o << " " << right << setw(w_p3_) << "-";
	// ||c||s
	if( m && s.feas_kkt_err().updated_k(0) )
		o << " " << setprecision(p3_) << right << setw(w_p3_) << s.feas_kkt_err().get_k(0);
	else
		o << " " << right << setw(w_p3_) << "-";
	// ||rGL||s
	if( s.opt_kkt_err().updated_k(0) )
		o << " " << setprecision(p3_) << right << setw(w_p3_) << s.opt_kkt_err().get_k(0);
	else
		o << " " << right << setw(w_p3_) << "-";
	// QN
	if( quasi_newt_stats ) {
		o << " " << right << setw(w_i2_);
		switch( quasi_newt_stats->updated() ) {
			case QuasiNewtonStats::UNKNOWN:
				o << "-";
				break;
			case QuasiNewtonStats:: REINITIALIZED:
				o << "IN";
				break;
			case QuasiNewtonStats::DAMPENED_UPDATED:
				o << "DU";
				break;
			case QuasiNewtonStats::UPDATED:
				o << "UP";
				break;
			case QuasiNewtonStats::SKIPED:
				o << "SK";
				break;
			case QuasiNewtonStats::INDEF_SKIPED:
				o << "IS";
				break;
			default:
				assert(0);
		}
	}
	else {
		o << " " << right << setw(w_i2_) << "-";
	}
	// #act
	if( s.nu().updated_k(0) )
		o << " " << right << setw(w_i4_) << s.nu().get_k(0).nz();
	else
		o	<< " " << right << setw(w_i4_) << "-";
	// ||Ypy||2
	if( m && s.Ypy().updated_k(0) )
		o << " "<< setprecision(p2_)  << right << setw(w_p2_) << s.Ypy().get_k(0).norm_2();
	else
		o << " " << right << setw(w_p2_) << "-";
	// ||Zpz||2
	if( s.Zpz().updated_k(0) )
		o << " " << setprecision(p2_) << right << setw(w_p2_) << s.Zpz().get_k(0).norm_2();
	else
		o << " " << right << setw(w_p2_) << "-";
	// ||d||inf
	if( s.d().updated_k(0) )
		o << " " << setprecision(p2_) << right << setw(w_p2_) << s.d().get_k(0).norm_inf();
	else
		o << " " << right << setw(w_p2_) << "-";
	// alpha
	if( s.alpha().updated_k(0) )
		o << " " << setprecision(p2_) << right << setw(w_p2_) << s.alpha().get_k(0);
	else
		o << " " << right << setw(w_p2_) << "-";

	o << std::endl;

	++printed_lines_;
}

void rSQPTrackConsoleStd::output_final( const Algorithm& p_algo
	, EAlgoReturn algo_return ) const
{
	using DynamicCastHelperPack::dyn_cast;

	const rSQPAlgo           &algo    = rsqp_algo(p_algo);
	const rSQPState          &s       = algo.rsqp_state();
	const NLPObjGradient     &nlp     = dyn_cast<const NLPObjGradient>(algo.nlp()); 
	const NLPFirstOrderInfo  *nlp_foi = dynamic_cast<const NLPFirstOrderInfo*>(&nlp); 
	
	const size_type
		m = nlp.m();

	std::ostream& o = this->o();

	// Output the table's header for the first iteration
	if(s.k() == 0) {
		print_top_header(s,algo);
		print_header(s,algo);
	}
	else {
		o
			<< " " << right << ul_i4_ 		// "k"
			<< " " << right << ul_p3_		// "f"
			<< " " << right << ul_p3_		// "||c||s"
			<< " " << right << ul_p3_		// "||rGL||s"
			<< " " << right << ul_i2_		// "QN"
			<< " " << right << ul_i4_		// "#act"
			<< endl;
	}

	// //////////////////////////////////////////
	// Output a row for the final iteration
	
	// Get a Quasi-Newton statistics.
	const QuasiNewtonStats	*quasi_newt_stats =
		( quasi_newton_stats_.exists_in(s) && quasi_newton_stats_(s).updated_k(0)
		  ? &quasi_newton_stats_(s).get_k(0)
		  : NULL );

	// k
	o << " " << right << setw(w_i4_) << s.k();
	// f
	if( s.f().updated_k(0) )
		o << " " << setprecision(p3_) << right << setw(w_p3_) << s.f().get_k(0);
	else
		o << " " << right << setw(w_p3_) << "-";
	// ||c||s
	if( m && s.feas_kkt_err().updated_k(0) )
		o << " " << setprecision(p3_) << right << setw(w_p3_) << s.feas_kkt_err().get_k(0);
	else
		o << " " << right << setw(w_p3_) << "-";
	// ||rGL||s
	if( s.opt_kkt_err().updated_k(0) )
		o << " " << setprecision(p3_) << right << setw(w_p3_) << s.opt_kkt_err().get_k(0);
	else
		o << " " << right << setw(w_p3_) << "-";
	// QN
	if( quasi_newt_stats ) {
		o << " " << right << setw(w_i2_);
		switch( quasi_newt_stats->updated() ) {
			case QuasiNewtonStats::UNKNOWN:
				o << "-";
				break;
			case QuasiNewtonStats:: REINITIALIZED:
				o << "IN";
				break;
			case QuasiNewtonStats::DAMPENED_UPDATED:
				o << "DU";
				break;
			case QuasiNewtonStats::UPDATED:
				o << "UP";
				break;
			case QuasiNewtonStats::SKIPED:
				o << "SK";
				break;
			case QuasiNewtonStats::INDEF_SKIPED:
				o << "IS";
				break;
			default:
				assert(0);
		}
	}
	else {
		o << " " << right << setw(w_i2_) << "-";
	}
	// #act
	if( s.nu().updated_k(0) )
		o << " " << right << setw(w_i4_) << s.nu().get_k(0).nz();
	else
		o	<< " " << right << setw(w_i4_) << "-";
	// ||Ypy||2
	if( m && s.Ypy().updated_k(0) )
		o << " "<< setprecision(p2_)  << right << setw(w_p2_) << s.Ypy().get_k(0).norm_2();
	else
		o << " " << right << setw(w_p2_) << "-";
	// ||Zpz||2
	if( s.Zpz().updated_k(0) )
		o << " " << setprecision(p2_) << right << setw(w_p2_) << s.Zpz().get_k(0).norm_2();
	else
		o << " " << right << setw(w_p2_) << "-";
	// ||d||inf
	if( s.d().updated_k(0) )
		o << " " << setprecision(p2_) << right << setw(w_p2_) << s.d().get_k(0).norm_inf();
	else
		o << " " << right << setw(w_p2_) << "-";

	o << endl;

	// Print total time
	o << "\nTotal time = " << timer_.read() << " sec\n";

	switch( algo_return ) {
		case IterationPack::TERMINATE_TRUE:
			o << "\nJackpot! You have found the solution!!!!!!\n";
			break;
		case IterationPack::TERMINATE_FALSE:
			o << "\nOops!  Not the solution.  Some error has occured!\n";
			break;
		case IterationPack::MAX_ITER_EXCEEDED:
			o << "\nOops!  Not the solution.  Maximum number of SQP iteration exceeded!\n";
			break;
		case IterationPack::MAX_RUN_TIME_EXCEEDED:
			o << "\nOops!  Not the solution.  Maximum runtime exceeded!\n";
			break;
		default:
			assert(0);
	}

	o	<< "\nNumber of function evaluations:\n"
		<<     "-------------------------------\n"
		<< "f(x)  : " << nlp.num_f_evals() << endl
		<< "c(x)  : " << ( m ? nlp.num_c_evals() : 0 ) << endl
		<< "Gf(x) : " << nlp.num_Gf_evals() << endl
		<< "Gc(x) : ";
	if(m){
		if( nlp_foi )
			o << nlp_foi->num_Gc_evals();
		else
			o << "?";
	}
	else {
		o << 0;
	}
	o << endl;
}

void rSQPTrackConsoleStd::print_top_header(const rSQPState &s
	, const rSQPAlgo &algo) const
{
	std::ostream& o = this->o();

	rSQPState::space_c_ptr_t
		space_c = s.get_space_c();

	o	<< "\n\n********************************\n"
		<< "*** Start of rSQP Iterations ***\n"
		<< "n = " << s.space_x().dim()
		<< ", m = " << ( space_c.get() ? space_c->dim() : 0 )
		<< ", nz = ";
	try {
		if(space_c.get()) {
			if( s.Gc().updated_k(0) )
				o	<< s.Gc().get_k(0).nz() << endl;
			else
				o	<< "?\n";
		}
		else {
			o	<< 0 << endl;
		}
	}
	catch( const AlgorithmState::DoesNotExist& ) {
			o	<< "?\n";
	}
	if( algo.nlp().scale_f() != 1.0 ) {
		o	<< "f(x) is scaled by : " << algo.nlp().scale_f() << endl;
	}
}

void rSQPTrackConsoleStd::print_header(const rSQPState &s
	, const rSQPAlgo &algo) const
{
	std::ostream& o = this->o();

	o
		<< endl
		<< " " << left << setw(w_i4_) << "k"
		<< " " << left << setw(w_p3_) << "f"
		<< " " << left << setw(w_p3_) << "||c||s"
		<< " " << left << setw(w_p3_) << "||rGL||s"
		<< " " << left << setw(w_i2_) << "QN"
		<< " " << left << setw(w_i4_) << "#act"
		<< " " << left << setw(w_p2_) << "||Ypy||2"
		<< " " << left << setw(w_p2_) << "||Zpz||2"
		<< " " << left << setw(w_p2_) << "||d||inf"
		<< " " << left << setw(w_p2_) << "alpha"
		<< endl
		<< " " << right << ul_i4_ 		// "k"
		<< " " << right << ul_p3_		// "f"
		<< " " << right << ul_p3_		// "||c||s"
		<< " " << right << ul_p3_		// "||rGL||s"
		<< " " << right << ul_i2_		// "QN"
		<< " " << right << ul_i4_		// "#act"
		<< " " << right << ul_p2_		// "||Ypy||2"
		<< " " << right << ul_p2_		// "||Zpz||2"
		<< " " << right << ul_p2_		// "||d||inf"
		<< " " << right << ul_p2_		// "alpha"
		<< endl;
}

} // end namespace ReducedSpaceSQPPack
