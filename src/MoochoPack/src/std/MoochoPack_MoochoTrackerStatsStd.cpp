// ////////////////////////////////////////////////////////////////////////////
// rSQPTrackStatsStd.cpp
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

#include "ReducedSpaceSQPPack/src/std/rSQPTrackStatsStd.h"
#include "ReducedSpaceSQPPack/src/rSQPState.h"
#include "ReducedSpaceSQPPack/src/rsqp_algo_conversion.h"
#include "NLPInterfacePack/src/NLPFirstOrderInfo.h"
#include "AbstractLinAlgPack/src/VectorWithOp.h"
#include "dynamic_cast_verbose.h"

namespace {
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
} // end namespace

namespace ReducedSpaceSQPPack {

using std::endl;
using std::setw;
using std::left;
using std::right;
using std::setprecision;

rSQPTrackStatsStd::rSQPTrackStatsStd(
	const ostream_ptr_t& o, const ostream_ptr_t& journal_out
	)
	: rSQPTrack(journal_out)
{
	set_output_stream(o);
}

void rSQPTrackStatsStd::set_output_stream(const ostream_ptr_t& o)
{
	o_ = o;
}

const rSQPTrackStatsStd::ostream_ptr_t&
rSQPTrackStatsStd::get_output_stream() const
{
	return o_;
}

void rSQPTrackStatsStd::initialize()
{
	num_QN_updates_ = 0;
	timer_.reset();
	timer_.start();
}

void rSQPTrackStatsStd::output_iteration(const Algorithm& p_algo) const
{
	const rSQPAlgo  &algo = rsqp_algo(p_algo);
	const rSQPState &s    = algo.rsqp_state();

	// All we have to do here is to just to count the number of quasi-newton updates
	const QuasiNewtonStats	*quasi_newt_stats =
		( quasi_newton_stats_.exists_in(s) && quasi_newton_stats_(s).updated_k(0)
		  ? &quasi_newton_stats_(s).get_k(0)
		  : NULL );
	if( quasi_newt_stats ) {
		QuasiNewtonStats::EUpdate updated = quasi_newt_stats->updated();
		if( updated == QuasiNewtonStats::DAMPENED_UPDATED || updated == QuasiNewtonStats::UPDATED )
			num_QN_updates_++;
	}
}

void rSQPTrackStatsStd::output_final( const Algorithm& p_algo
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

	// Stop the timer
	timer_.stop();

	// Formating info
	const int
		p      = 18,
		stat_w = 15,
		val_w  = p + 10;

	// Get a Quasi-Newton statistics.
	const QuasiNewtonStats	*quasi_newt_stats =
		( quasi_newton_stats_.exists_in(s) && quasi_newton_stats_(s).updated_k(0)
		  ? &quasi_newton_stats_(s).get_k(0)
		  : NULL );
	if( quasi_newt_stats ) {
		QuasiNewtonStats::EUpdate updated = quasi_newt_stats->updated();
		if( updated == QuasiNewtonStats::DAMPENED_UPDATED || updated == QuasiNewtonStats::UPDATED )
			num_QN_updates_++;
	}

	// status
	o << left << setw(stat_w) << "status" << "= "
		<< right << setw(val_w);
	switch( algo_return ) {
		case GeneralIterationPack::TERMINATE_TRUE:
			o << "solved";
			break;
		case GeneralIterationPack::TERMINATE_FALSE:
			o << "except";
			break;
		case GeneralIterationPack::MAX_ITER_EXCEEDED:
			o << "max_iter";
			break;
		case GeneralIterationPack::MAX_RUN_TIME_EXCEEDED:
			o << "max_run_time";
			break;
		default:
			assert(0);
	}
	o << "; # solved, except, max_iter, max_run_time\n";
	// niter
	o << left << setw(stat_w) << "niter" << "= "
		<< right << setw(val_w) << s.k()
		<< "; # Number of rSQP iterations (plus 1?)\n";
	// nfunc
	o << left << setw(stat_w) << "nfunc" << "= "
		<< right << setw(val_w) << my_max(nlp.num_f_evals(),(m? nlp.num_c_evals():0) )
		<< "; # max( number f(x) evals, number c(x) evals )\n";
	// ngrad
	o << left << setw(stat_w) << "ngrad" << "= "
		<< right << setw(val_w) << my_max(nlp.num_Gf_evals(),(m?(nlp_foi?nlp_foi->num_Gc_evals():s.k()+1):0))
		<< "; # max( number Gf(x) evals, number Gc(x) evals )\n";
	// CPU
	o << left << setw(stat_w) << "CPU" << "= "
		<< right << setw(val_w) << timer_.read()
		<< "; # Number of CPU seconds total\n";
	// obj_func
	o << left << setw(stat_w) << "obj_func" << "= "
		<< right << setw(val_w);
	if(s.f().updated_k(0))
		o << s.f().get_k(0);
	else
		o << "-";
	o << "; # Objective function value f(x) at final point\n";
	// feas_kkt_err
	o << left << setw(stat_w) << "feas_kkt_err" << "= "
		<< right << setw(val_w);
	if(s.feas_kkt_err().updated_k(0))
		o << s.feas_kkt_err().get_k(0);
	else if(s.c().updated_k(0))
		o << s.c().get_k(0).norm_inf();
	else
		o << "-";
	o << "; # Feasibility error at final point (scaled ||c(x)||inf, feas_err_k)\n";
	// opt_kkt_err
	o << left << setw(stat_w) << "opt_kkt_err" << "= "
		<< right << setw(val_w);
	if(s.opt_kkt_err().updated_k(0))
		o << s.opt_kkt_err().get_k(0);
	else if(s.rGL().updated_k(0))
		o << s.rGL().get_k(0).norm_inf();
	else if(s.rGL().updated_k(-1))
		o << s.rGL().get_k(-1).norm_inf();
	else
		o << "-";
	o << "; # Optimality error at final point (scaled ||rGL||inf, opt_err_k)\n";
	// nact
	o << left << setw(stat_w) << "nact" << "= "
		<< right << setw(val_w);
	if(s.nu().updated_k(0))
		o << s.nu().get_k(0).nz();
	else if(s.nu().updated_k(-1))
		o << s.nu().get_k(-1).nz();
	else
		o << "-";
	o << "; # Number of total active constraints at the final point\n";
	// nbasis_change
	const IterQuantityAccess<index_type> &num_basis = s.num_basis();
	const int lu_k = num_basis.last_updated();
	o << left << setw(stat_w) << "nbasis_change" << "= "
		<< right << setw(val_w) << ( lu_k != IterQuantity::NONE_UPDATED
									 ? num_basis.get_k(lu_k)
									 : 0 ) 
		<< "; # Number of basis changes\n";
	// nquasi_newton
	o << left << setw(stat_w) << "nquasi_newton" << "= "
		<< right << setw(val_w) << num_QN_updates_
		<< "; # Number of quasi-newton updates\n";

}

} // end namespace ReducedSpaceSQPPack
