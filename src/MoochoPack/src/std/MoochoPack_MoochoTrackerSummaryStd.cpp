// ////////////////////////////////////////////////////////////////////////////
// rSQPTrackSummaryStd.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>
#include <iomanip>

#include "../../include/std/rSQPTrackSummaryStd.h"
#include "../../include/rSQPState.h"
#include "../../include/rsqp_algo_conversion.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "NLPInterfacePack/include/NLPFirstOrderInfo.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/GenMatrixOut.h"
#include "LinAlgPack/include/VectorOut.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/PermOut.h"
#include "Misc/include/dynamic_cast_verbose.h"

using std::endl;
using std::setw;

void ReducedSpaceSQPPack::rSQPTrackSummaryStd::output_iteration(const Algorithm& algo) const
{
	using LinAlgPack::norm_2;
	using LinAlgPack::norm_inf;

	const rSQPState &s = rsqp_algo(algo).rsqp_state();
	
	int w = 15;
	int prec = 6;
	o().precision(prec);

	// Output the table's header for the first iteration
	if(s.k() == 0) {
		print_header(s);
	}

	// ///////////////////////////////
	// Output a row for the iteration
	
	// Get active set and QP solver statistics.
	const ActSetStats		*act_stats =
		( act_set_stats_(s).updated_k(0) ? &act_set_stats_(s).get_k(0) : 0 );
	const QPSolverStats		*qp_stats =
		( qp_solver_stats_(s).updated_k(0) ? &qp_solver_stats_(s).get_k(0) : 0 );
	const QuasiNewtonStats	*quasi_newt_stats =
		( quasi_newton_stats_(s).updated_k(0) ? &quasi_newton_stats_(s).get_k(0) : 0 );

	// Get the norms of Ypy and Zpz
	value_type norm_2_Ypy = -1.0, norm_2_Zpz = -1.0;
	bool Ypy_exists, Zpz_exists;
	if( Ypy_exists = s.Ypy().updated_k(0) )
		norm_2_Ypy = s.Ypy().get_k(0).norm_2();
	if( Zpz_exists = s.Zpz().updated_k(0) )
		norm_2_Zpz = s.Zpz().get_k(0).norm_2();

	o()	<< std::right
		<< setw(5) << s.k();

	if( s.f().updated_k(0) )
		o() << setw(w) << s.f().get_k(0);
	else
		o() << setw(w) << "-";

	if( s.Gf().updated_k(0) )
		o() << setw(w) << s.Gf().get_k(0).norm_inf();
	else
		o() << setw(w) << "-";

	if( s.c().updated_k(0) )
		o() << setw(w)
			<< s.c().get_k(0).norm_inf();
	else
		o() << setw(w) << "-";

	{
		const rSQPState::IQA_Vector
			&rGL_GL = ( opt_error_ == OPT_ERROR_REDUCED_GRADIENT_LAGR
							? s.rGL() : s.GL()  );
		if( rGL_GL.updated_k(0) )
			o() << setw(w) << rGL_GL.get_k(0).norm_inf();
		else
			o() << setw(w) << "-";
	}

	if( quasi_newt_stats ) {
		o() << setw(w);
		switch( quasi_newt_stats->updated() ) {
			case QuasiNewtonStats::UNKNOWN:
				o() << "-";
				break;
			case QuasiNewtonStats:: REINITIALIZED:
				o() << "initialized";
				break;
			case QuasiNewtonStats::DAMPENED_UPDATED:
				o() << "damp.updated";
				break;
			case QuasiNewtonStats::UPDATED:
				o() << "updated";
				break;
			case QuasiNewtonStats::SKIPED:
				o() << "skiped";
				break;
			case QuasiNewtonStats::INDEF_SKIPED:
				o() << "indef skiped";
				break;
			default:
				assert(0);
		}
	}
	else {
		o() << setw(w) << "-";
	}

	if( act_stats ) {
		o() << setw(7) << act_stats->num_active();
		// don't know num_add and num_drops on first iteration.
		if( act_stats->num_adds() == ActSetStats::NOT_KNOWN ) { 
			o()	<< setw(7) << "-";
		}
		else {		
			o()	<< setw(7) << act_stats->num_adds();
		}
		if( act_stats->num_drops() == ActSetStats::NOT_KNOWN ) {
			o()	<< setw(7) << "-";
		}
	else {		
		o()	<< setw(7) << act_stats->num_drops();
	}
	}
	else {
		o()	<< setw(7) << "-"
			<< setw(7) << "-"
			<< setw(7) << "-";
	}

	if( qp_stats ) {
		o()	<< setw(7) << qp_stats->num_qp_iter()
			<< setw(3) << ( qp_stats->warm_start() ? 'w' : 'c')
			<< setw(2) << ( qp_stats->infeasible_qp() ? 'i' : 'f');
		num_total_qp_iter_ += qp_stats->num_qp_iter();
	}
	else {
	o()	<< setw(7) << "-"
		<< setw(3) << "-"
		<< setw(2) << "-";
	}

	if(Ypy_exists)
		o()	<< setw(w) << norm_2_Ypy;
	else
		o() << setw(w) << "-";

	if(Zpz_exists)
		o()	<< setw(w) << norm_2_Zpz;
	else
		o() << setw(w) << "-";

	if( s.d().updated_k(0) )
		o() << setw(w)
			<< s.d().get_k(0).norm_inf();
	else
		o() << setw(w) << "-";

	if( s.alpha().updated_k(0) )
		o() << setw(w) << s.alpha().get_k(0);
	else
		o() << setw(w) << "-";

	o() << std::endl;
}

void ReducedSpaceSQPPack::rSQPTrackSummaryStd::output_final(const Algorithm& algo
	, EAlgoReturn algo_return) const
{
	using DynamicCastHelperPack::const_dyn_cast;
	using LinAlgPack::norm_inf;

	const rSQPAlgo			&_algo = rsqp_algo(algo);
	const rSQPState			&s =_algo.rsqp_state();
	const NLPFirstOrderInfo
#ifdef _WINDOWS
		&nlp = dynamic_cast<const NLPFirstOrderInfo&>(_algo.nlp()); 
#else
		&nlp = const_dyn_cast<NLPFirstOrderInfo>(_algo.nlp()); 
#endif
	int w = 15;
	int prec = 6;
	o().precision(prec);

	// Get active set and QP solver statistics.
	const IterQuantityAccess<ActSetStats>		&act_stats_iq = act_set_stats_(s);
	const IterQuantityAccess<QPSolverStats>		&qp_stats_iq = qp_solver_stats_(s);
	const IterQuantityAccess<QuasiNewtonStats>	&quasi_newt_stats_iq = quasi_newton_stats_(s);

	// Output the table's header for the first iteration
	if(s.k() == 0) {
		print_header(s);
	}
	else {
		o()	<< " ----"
			<< "   ------------"
			<< "   ------------"
			<< "   ------------"
			<< "   ------------"
			<< "   ------------"
			<< " ------"
			<< " ------"
			<< " ------"
			<< " ------"
			<< " ----\n";
	}

	o()	<< std::right
		<< setw(5) << s.k();

	if( s.f().updated_k(0) )
		o() << setw(w) << s.f().get_k(0);
	else
		o() << setw(w) << "-";

	if( s.Gf().updated_k(0) )
		o() << setw(w) << s.Gf().get_k(0).norm_inf();
	else
		o() << setw(w) << "-";

	if( s.c().updated_k(0) )
		o() << setw(w)
			<< s.c().get_k(0).norm_inf();
	else
		o() << setw(w) << "-";

	{
		const rSQPState::IQA_Vector
			&rGL_GL = ( opt_error_ == OPT_ERROR_REDUCED_GRADIENT_LAGR
							? s.rGL() : s.GL()  );
		if( rGL_GL.updated_k(0) )
			o() << setw(w) << rGL_GL.get_k(0).norm_inf();
		else
			o() << setw(w) << "-";
	}

	if( quasi_newt_stats_iq.updated_k(0) ) {
		const QuasiNewtonStats &stats = quasi_newt_stats_iq.get_k(0);
		o() << setw(w);
		switch( stats.updated() ) {
			case QuasiNewtonStats::UNKNOWN:
				o() << "-";
				break;
			case QuasiNewtonStats:: REINITIALIZED:
				o() << "initialized";
				break;
			case QuasiNewtonStats::DAMPENED_UPDATED:
				o() << "damp.updated";
				break;
			case QuasiNewtonStats::UPDATED:
				o() << "updated";
				break;
			case QuasiNewtonStats::SKIPED:
				o() << "skiped";
				break;
			case QuasiNewtonStats::INDEF_SKIPED:
				o() << "indef skiped";
				break;
			default:
				assert(0);
		}
	}
	else {
		o()	<< setw(w) << "-";;
	}

	if( act_stats_iq.updated_k(0) ) {
		const ActSetStats &act_stats = act_stats_iq.get_k(0);
		o()	<< setw(7) << act_stats.num_active()
			<< setw(7) << act_stats.num_adds()
			<< setw(7) << act_stats.num_drops();
	}
	else {
		o()	<< setw(7) << "-"
			<< setw(7) << "-"
			<< setw(7) << "-";
	}

	if( qp_stats_iq.updated_k(0) ) {
		const QPSolverStats &qp_stats = qp_stats_iq.get_k(0);
		o() << setw(7) << qp_stats.num_qp_iter()
			<< setw(3) << ( qp_stats.warm_start() ? 'w' : 'c')
			<< setw(2) << ( qp_stats.infeasible_qp() ? 'i' : 'f');
		num_total_qp_iter_ += qp_stats.num_qp_iter();
	}
	else {
		o()	<< setw(5) << "-";
	}

	if(s.Ypy().updated_k(0))
		o()	<< setw(w) 
			<< s.Ypy().get_k(0).norm_2();
	else
		o() << setw(w) << "-";

	if(s.Zpz().updated_k(0))
		o()	<< setw(w)
			<< s.Zpz().get_k(0).norm_2();
	else
		o() << setw(w) << "-";

	if( s.d().updated_k(0) )
		o() << setw(w)
			<< s.d().get_k(0).norm_inf();
	else
		o() << setw(w) << "-";

	o()	<< "\n\nnumber of function evaluations:\n"
		<<     "-------------------------------\n"
		<< "f(x)  : " << nlp.num_f_evals() << endl
		<< "c(x)  : " << nlp.num_c_evals() << endl
		<< "Gf(x) : " << nlp.num_Gf_evals() << endl
		<< "Gc(x) : " << nlp.num_Gc_evals() << endl;
}

void ReducedSpaceSQPPack::rSQPTrackSummaryStd::print_header(const rSQPState &s) const
{
	// Reset the count of total QP iterations
	num_total_qp_iter_ = 0;

	int w = 15;
	int prec = 6;

	o()	<< "\n\n********************************\n"
		<< "*** Start of rSQP Iterations ***\n"
		<< "n = " << s.x().get_k(0).cv().size()
		<< ", m = " << s.c().get_k(0).cv().size()
		<< ", nz = ";
	if( s.Gc().updated_k(0) )
		o()	<< s.Gc().get_k(0).nz() << endl;
	else
		o()	<< "?\n";
	o()
		<< "\n k   "
		<< "   f           "
		<< "   ||Gf||inf   "
		<< "   ||c||inf    ";
	switch(opt_error_) {
		case OPT_ERROR_REDUCED_GRADIENT_LAGR:
			o()	<< "   ||rGL||inf  ";
			break;
		case OPT_ERROR_GRADIENT_LAGR:
			o()	<< "   ||GL||inf  ";
			break;
		default:
			assert(0);
	}
	o()	<< "   quasi-Newton"
		<< " #act  "
		<< " #adds "
		<< " #drops"
		<< " #qpitr"
		<< " wcfi  "
		<< "   ||Ypy||2    "
		<< "   ||Zpz||2    "
		<< "   ||d||inf    "
		<< "   alpha\n"
		<< " ----"
		<< "   ------------"
		<< "   ------------"
		<< "   ------------"
		<< "   ------------"
		<< "   ------------"
		<< " ------"
		<< " ------"
		<< " ------"
		<< " ------"
		<< " ----"
		<< "   ------------"
		<< "   ------------"
		<< "   ------------"
		<< "   ------------\n";
}