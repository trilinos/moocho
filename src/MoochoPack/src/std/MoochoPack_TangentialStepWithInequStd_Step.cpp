// ////////////////////////////////////////////////////////////////////////////
// NullSpaceStepWithInequStd_Step.cpp
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
#include <sstream>

#include "ReducedSpaceSQPPack/include/std/NullSpaceStepWithInequStd_Step.h"
#include "ReducedSpaceSQPPack/include/rsqp_algo_conversion.h"
#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackExceptions.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/MatrixIdentConcat.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/VectorAuxiliaryOps.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "dynamic_cast_verbose.h"

namespace ReducedSpaceSQPPack {

NullSpaceStepWithInequStd_Step::NullSpaceStepWithInequStd_Step(
	const qp_solver_ptr_t       &qp_solver
	,const qp_tester_ptr_t      &qp_tester
	,value_type                 warm_start_frac
	,EQPTesting                 qp_testing
	,bool                       primal_feasible_point_error
	,bool                       dual_feasible_point_error
	)
	:qp_solver_(qp_solver)
	,qp_tester_(qp_tester)
	,warm_start_frac_(warm_start_frac)
	,qp_testing_(qp_testing)
	,primal_feasible_point_error_(primal_feasible_point_error)
	,dual_feasible_point_error_(dual_feasible_point_error)
	,dl_iq_(dl_name)
	,du_iq_(du_name)
{}

bool NullSpaceStepWithInequStd_Step::do_step(
	Algorithm& _algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss
	)
{
	namespace mmp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	using ::fabs;
	using LinAlgOpPack::Vt_S;
	using LinAlgOpPack::V_VpV;
	using LinAlgOpPack::V_VmV;
	using LinAlgOpPack::Vp_StV;
	using LinAlgOpPack::Vp_V;
	using LinAlgOpPack::V_StV;
	using LinAlgOpPack::V_MtV;
//	using ConstrainedOptimizationPack::min_abs;
	using AbstractLinAlgPack::max_near_feas_step;
	typedef VectorWithOpMutable::vec_mut_ptr_t   vec_mut_ptr_t;

	rSQPAlgo             &algo         = rsqp_algo(_algo);
	rSQPState            &s            = algo.rsqp_state();
	EJournalOutputLevel  olevel        = algo.algo_cntr().journal_output_level();
	std::ostream         &out          = algo.track().journal_out();
	const bool           check_results = algo.algo_cntr().check_results();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// problem dimensions
	const size_type
		n  = algo.nlp().n(),
		m  = algo.nlp().m(),
		mI = algo.nlp().mI(),
		r  = s.equ_decomp().size();

	// Get the iteration quantity container objects
	IterQuantityAccess<value_type>
		&alpha_iq = s.alpha(),
		&zeta_iq  = s.zeta(),
		&eta_iq   = s.eta();
	IterQuantityAccess<VectorWithOpMutable>
		&x_iq       = s.x(),
		&dl_iq      = dl_iq_(s),
		&du_iq      = du_iq_(s),
		&nu_iq      = s.nu(),
		*c_iq       = m > 0  ? &s.c()       : NULL,
		*lambda_iq  = m > 0  ? &s.lambda()  : NULL,
		*h_iq       = mI > 0 ? &s.h()       : NULL,
		*lambdaI_iq = mI > 0 ? &s.lambdaI() : NULL,
		&rGf_iq     = s.rGf(),
		&w_iq       = s.w(),
		&qp_grad_iq = s.qp_grad(),
		&py_iq      = s.py(),
		&pz_iq      = s.pz(),
		&Ypy_iq     = s.Ypy(),
		&Zpz_iq     = s.Zpz();
	IterQuantityAccess<MatrixWithOp>
		&Z_iq   = s.Z(),
		*Y_iq   = m        ? &s.Y()  : NULL,
		*Uz_iq  = (m > r)  ? &s.Uz() : NULL,
		*Uy_iq  = (m > r)  ? &s.Uy() : NULL,
		*Vz_iq  = (mI > 0) ? &s.Vz() : NULL,
		*Vy_iq  = (mI > 0) ? &s.Vy() : NULL;
	IterQuantityAccess<MatrixSymWithOp>
		&rHL_iq = s.rHL();
	IterQuantityAccess<ActSetStats>
		&act_set_stats_iq = act_set_stats_(s);
	
	// Accessed/modified/updated (just some)
	VectorWithOpMutable  *Ypy_k = (m ? &Ypy_iq.get_k(0) : NULL);
	const MatrixWithOp   &Z_k   = Z_iq.get_k(0);
	VectorWithOpMutable  &pz_k  = pz_iq.set_k(0);
	VectorWithOpMutable  &Zpz_k = Zpz_iq.set_k(0);

	// Comupte qp_grad which is an approximation to rGf + Z'*HL*Y*py

	// qp_grad = rGf
	VectorWithOpMutable
		&qp_grad_k = ( qp_grad_iq.set_k(0) = rGf_iq.get_k(0) );

	// qp_grad += zeta * w
	if( w_iq.updated_k(0) ) {
		if(zeta_iq.updated_k(0))
			Vp_StV( &qp_grad_k, zeta_iq.get_k(0), w_iq.get_k(0) );
		else
			Vp_V( &qp_grad_k, w_iq.get_k(0) );
	}

	// Set the bounds for:
	//   dl <= Z*pz + Y*py <= du  ->  dl - Ypy <= Z*pz <= du - Ypz

	vec_mut_ptr_t
		bl = s.space_x().create_member(),
		bu = s.space_x().create_member();

	if(m) {
		// bl = dl_k - Ypy_k
		V_VmV( bl.get(), dl_iq.get_k(0), *Ypy_k );
		// bu = du_k - Ypy_k
		V_VmV( bu.get(), du_iq.get_k(0), *Ypy_k );
	}
	else {
		*bl = dl_iq.get_k(0);
		*bu = du_iq.get_k(0);
	}

	// Print out the QP bounds for the constraints
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\nqp_grad_k = \n" << qp_grad_k;
		out << "\nbl = \n" << *bl;
		out << "\nbu = \n" << *bu;
	}

	//
	// Determine if we should perform a warm start or not.
	//
	bool do_warm_start = false;
	if( act_set_stats_iq.updated_k(-1) ) {
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out	<< "\nDetermining if the QP should use a warm start ...\n";
		}
		// We need to see if we should preform a warm start for the next iteration
		ActSetStats &stats = act_set_stats_iq.get_k(-1);
		const size_type
			num_active = stats.num_active(),
			num_adds   = stats.num_adds(),
			num_drops  = stats.num_drops();
		const value_type
			frac_same
			= ( num_adds == ActSetStats::NOT_KNOWN || num_active == 0
				? 0.0
				: std::_MAX(((double)(num_active)-num_adds-num_drops) / num_active, 0.0 ) );
		do_warm_start = ( num_active > 0 && frac_same >= warm_start_frac() );
		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << "\nnum_active = " << num_active;
			if( num_active ) {
				out	<< "\nmax(num_active-num_adds-num_drops,0)/(num_active) = "
					<< "max("<<num_active<<"-"<<num_adds<<"-"<<num_drops<<",0)/("<<num_active<<") = "
					<< frac_same;
				if( do_warm_start )
					out << " >= ";
				else
					out << " < ";
				out << "warm_start_frac = " << warm_start_frac();
			}
			if( do_warm_start )
				out << "\nUse a warm start this time!\n";
			else
				out << "\nDon't use a warm start this time!\n";
		}
	}

	// Use active set from last iteration as an estimate for current active set
	// if we are to use a warm start.
	// 
	// ToDo: If the selection of dependent and independent variables changes
	// then you will have to adjust this or not perform a warm start at all!
	if( do_warm_start ) {
		nu_iq.set_k(0,-1);
	}
	else {
		nu_iq.set_k(0) = 0.0; // No guess of the active set
	}
	VectorWithOpMutable  &nu_k  = nu_iq.get_k(0);

	//
	// Setup the reduced QP subproblem
	//
	// The call to the QP is setup for the more flexible call to the QPSolverRelaxed
	// interface to deal with the three independent variabilities: using simple
	// bounds for pz or not, general inequalities included or not, and extra equality
	// constraints included or not.
	// If this method of calling the QP solver were not used then 4 separate
	// calls to solve_qp(...) would have to be included to handle the four possible
	// QP formulations.
	//

	// The numeric arguments for the QP solver (in the nomenclatrue of QPSolverRelaxed)

	const value_type  qp_bnd_inf = NLP::infinite_bound();

	const VectorWithOp      &qp_g       = qp_grad_k;
	const MatrixSymWithOp   &qp_G       = rHL_iq.get_k(0);
	const value_type		qp_etaL     = 0.0;
	vec_mut_ptr_t           qp_dL       = mmp::null;
	vec_mut_ptr_t           qp_dU       = mmp::null;
	mmp::ref_count_ptr<const MatrixWithOp>
                            qp_E        = mmp::null;
	BLAS_Cpp::Transp        qp_trans_E  = BLAS_Cpp::no_trans;
	vec_mut_ptr_t           qp_b        = mmp::null;
	vec_mut_ptr_t           qp_eL       = mmp::null;
	vec_mut_ptr_t           qp_eU       = mmp::null;
	mmp::ref_count_ptr<const MatrixWithOp>
	                        qp_F        = mmp::null;
	BLAS_Cpp::Transp        qp_trans_F  = BLAS_Cpp::no_trans;
	vec_mut_ptr_t           qp_f        = mmp::null;
	value_type				qp_eta      = 0.0;
	VectorWithOpMutable     &qp_d       = pz_k;  // pz_k will be updated directly!
	vec_mut_ptr_t           qp_nu       = mmp::null;
	vec_mut_ptr_t           qp_mu       = mmp::null;
	vec_mut_ptr_t           qp_Ed       = mmp::null;
	vec_mut_ptr_t           qp_lambda   = mmp::null;

	//
	// Determine if we can use simple bounds on pz.
	// 
	// If we have a variable-reduction null-space matrix
	// (with any choice for Y) then:
	// 
	// d = Z*pz + (1-eta) * Y*py
	// 
	// [ d(var_dep)   ]  = [ D ] * pz  + (1-eta) * [ Ypy(var_dep)   ]
	// [ d(var_indep) ]    [ I ]                   [ Ypy(var_indep) ]
	// 
	// For a cooridinate decomposition (Y = [ I ; 0 ]) then Ypy(var_indep) == 0.0 and
	// in this case the bounds on d(var_indep) become simple bounds on pz even
	// with the relaxation.
	// 
	// Otherwise, we can not use simple variable bounds and implement the
	// relaxation properly.
	// 

	const MatrixIdentConcat
		*Zvr = dynamic_cast<const MatrixIdentConcat*>( &Z_k );
	const Range1D
		var_dep   = Zvr ? Zvr->D_rng() : Range1D::Invalid,
		var_indep = Zvr ? Zvr->I_rng() : Range1D();

	const value_type Ypy_indep_norm_inf = ( m ? Ypy_k->sub_view(var_indep)->norm_inf() : 0.0);

	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
		out
			<< "\nDetermine if we can use simple bounds on pz ...\n"
			<< "    m = " << m << std::endl
			<< "    dynamic_cast<const MatrixIdentConcat*>(&Z_k) = " << Zvr << std::endl
			<< "    ||Ypy_k(var_indep)||inf = " << Ypy_indep_norm_inf << std::endl;

	const bool
		use_simple_pz_bounds = ( m == 0 || ( Zvr != NULL && Ypy_indep_norm_inf == 0.0 ) );

	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
		out
			<< (use_simple_pz_bounds ? "\nUsing simple bounds on pz ...\n" : "\nUsing bounds on full Z*pz ...\n");

	if( mI && use_simple_pz_bounds ) {
		assert(0); // ToDo: Implement general inequalities!
	}
	else if( mI && !use_simple_pz_bounds ) {	
		assert(0); // ToDo: Implement general inequalities!
	}
	else if( !mI && use_simple_pz_bounds ) {
		// Set simple bound constraints on pz
		qp_dL = bl->sub_view(var_indep);
		qp_dU = bu->sub_view(var_indep);
		qp_nu = nu_k.sub_view(var_indep); // nu_k(var_indep) will be updated directly!
		if(m) {
			// Set general inequality constraints for D*pz
			qp_E   = mmp::rcp(&Zvr->D(),false);
			qp_b   = Ypy_k->sub_view(var_dep);
			qp_eL  = bl->sub_view(var_dep);
			qp_eU  = bu->sub_view(var_dep);
			qp_mu  = nu_k.sub_view(var_dep);  // nu_k(var_dep) will be updated directly!
			qp_Ed  = Zpz_k.sub_view(var_dep); // Zpz_k(var_dep) will be updated directly!
		}
	}
	else if( !mI && !use_simple_pz_bounds ) {
		// There are no simple bounds! (leave qp_dL, qp_dU and qp_nu as null)
		// Set general inequality constraints for Z*pz
		qp_E   = mmp::rcp(&Z_k,false);
		qp_b   = mmp::rcp(Ypy_k,false);
		qp_eL  = bl;
		qp_eU  = bu;
		qp_mu  = mmp::rcp(&nu_k,false);
		qp_Ed  = mmp::rcp(&Zpz_k,false); // Zpz_k will be updated directly!
	}
	else {
		assert(0);
	}

	// Set the general equality constriants (if they exist)
	Range1D equ_undecomp = s.equ_undecomp();
	if( m > r && m > 0 ) {
		// qp_f = Uy_k * py_k + c_k(equ_undecomp)
		qp_f = s.space_c().sub_space(equ_undecomp)->create_member();
		V_MtV( qp_f.get(), Uy_iq->get_k(0), BLAS_Cpp::no_trans, py_iq.get_k(0) );
		Vp_V( qp_f.get(), *c_iq->get_k(0).sub_view(equ_undecomp) );
		// Must resize for the undecomposed constriants if it has not already been
		qp_F       = mmp::rcp(&Uy_iq->get_k(0),false);
		qp_lambda  = lambda_iq->set_k(0).sub_view(equ_undecomp); // lambda_k(equ_undecomp), will be updated directly!
	}

	// Setup the rest of the arguments

	QPSolverRelaxed::EOutputLevel
		qp_olevel;
	switch( olevel ) {
		case PRINT_NOTHING:
			qp_olevel = QPSolverRelaxed::PRINT_NONE;
			break;
		case PRINT_BASIC_ALGORITHM_INFO:
			qp_olevel = QPSolverRelaxed::PRINT_NONE;
			break;
		case PRINT_ALGORITHM_STEPS:
			qp_olevel = QPSolverRelaxed::PRINT_BASIC_INFO;
			break;
		case PRINT_ACTIVE_SET:
			qp_olevel = QPSolverRelaxed::PRINT_ITER_SUMMARY;
			break;
		case PRINT_VECTORS:
			qp_olevel = QPSolverRelaxed::PRINT_ITER_VECTORS;
			break;
		case PRINT_ITERATION_QUANTITIES:
			qp_olevel = QPSolverRelaxed::PRINT_EVERY_THING;
			break;
		default:
			assert(0);
	}

	//
	// Solve the QP
	// 
	qp_solver().infinite_bound(qp_bnd_inf);
	const QPSolverStats::ESolutionType
		solution_type =
		qp_solver().solve_qp(
			int(olevel) == int(PRINT_NOTHING) ? NULL : &out
			,qp_olevel
			,( algo.algo_cntr().check_results()
			   ? QPSolverRelaxed::RUN_TESTS :  QPSolverRelaxed::NO_TESTS )
			,qp_g, qp_G, qp_etaL, qp_dL.get(), qp_dU.get()
			,qp_E.get(), qp_trans_E, qp_E.get() ? qp_b.get() : NULL
			,qp_E.get() ? qp_eL.get() : NULL, qp_E.get() ? qp_eU.get() : NULL
			,qp_F.get(), qp_trans_F, qp_F.get() ? qp_f.get() : NULL
			,NULL // obj_d
			,&qp_eta, &qp_d
			,qp_nu.get()
			,qp_mu.get(), qp_E.get() ? qp_Ed.get() : NULL
			,qp_F.get() ? qp_lambda.get() : NULL
			,NULL // qp_Fd
			);

	//
	// Check the optimality conditions for the QP
	//
	std::ostringstream omsg;
	bool throw_qp_failure = false;
	if(		qp_testing() == QP_TEST
		|| ( qp_testing() == QP_TEST_DEFAULT && algo.algo_cntr().check_results() )  )
	{
		if( int(olevel) >= int(PRINT_ALGORITHM_STEPS) ) {
			out	<< "\nChecking the optimality conditions of the reduced QP subproblem ...\n";
		}
		if(!qp_tester().check_optimality_conditions(
			solution_type,qp_solver().infinite_bound()
			,int(olevel) == int(PRINT_NOTHING) ? NULL : &out
			,int(olevel) >= int(PRINT_VECTORS) ? true : false
			,int(olevel) >= int(PRINT_ITERATION_QUANTITIES) ? true : false
			,qp_g, qp_G, qp_etaL, qp_dL.get(), qp_dU.get()
			,qp_E.get(), qp_trans_E, qp_E.get() ? qp_b.get() : NULL
			,qp_E.get() ? qp_eL.get() : NULL, qp_E.get() ? qp_eU.get() : NULL
			,qp_F.get(), qp_trans_F, qp_F.get() ? qp_f.get() : NULL
			,NULL // obj_d
			,&qp_eta, &qp_d
			,qp_nu.get()
			,qp_mu.get(), qp_E.get() ? qp_Ed.get() : NULL
			,qp_F.get() ? qp_lambda.get() : NULL
			,NULL // qp_Fd
			))
		{
			omsg << "\n*** Alert! at least one of the QP optimality conditions did not check out.\n";
			if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out << omsg.str();
			}
			throw_qp_failure = true;
		}
	}

	//
	// Set the solution
	//
	if( mI && use_simple_pz_bounds ) {
		assert(0); // ToDo: Implement general inequalities!
	}
	else if( mI && !use_simple_pz_bounds ) {	
		assert(0); // ToDo: Implement general inequalities!
	}
	else if( !mI && !use_simple_pz_bounds ) {
		// Everything is already updated!
	}
	else if( !mI && use_simple_pz_bounds ) {
		// Just have to set Zpz_k(var_indep) = pz_k
		*Zpz_k.sub_view(var_indep) = pz_k;
	}
	else {
		assert(0);
	}

	// Set the solution statistics
	qp_solver_stats_(s).set_k(0) = qp_solver().get_qp_stats();

	// Cut back Ypy_k = (1-eta) * Ypy_k
	const value_type eps = std::numeric_limits<value_type>::epsilon();
	if( fabs(qp_eta - 0.0) > eps ) {
		if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out
				<< "\n*** Alert! the QP was infeasible (eta = "<<qp_eta<<").  Cutting back Ypy_k = (1.0 - eta)*Ypy_k  ...\n";
		}
		Vt_S( Ypy_k, 1.0 - qp_eta );
	}

	// eta_k
	eta_iq.set_k(0) = qp_eta;

	//
	// Modify the solution if we have to!
	// 
	switch(solution_type) {
		case QPSolverStats::OPTIMAL_SOLUTION:
			break;	// we are good!
		case QPSolverStats::PRIMAL_FEASIBLE_POINT:
		{
			omsg
				<< "\n*** Alert! the returned QP solution is PRIMAL_FEASIBLE_POINT but not optimal!\n";
			if( primal_feasible_point_error() )
				omsg
					<< "\n*** primal_feasible_point_error == true, this is an error!\n";
			if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out << omsg.str();
			}
			throw_qp_failure = primal_feasible_point_error();
			break;
		}	
		case QPSolverStats::DUAL_FEASIBLE_POINT:
		{
			omsg
				<< "\n*** Alert! the returned QP solution is DUAL_FEASIBLE_POINT"
				<< "\n*** but not optimal so we cut back the step ...\n";
			if( dual_feasible_point_error() )
				omsg
					<< "\n*** dual_feasible_point_error == true, this is an error!\n";
			if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out << omsg.str();
			}
			// Cut back the step to fit in the bounds
			// 
			// dl <= u*(Ypy_k+Zpz_k) <= du
			//
			vec_mut_ptr_t
				zero  = s.space_x().create_member(0.0),
				d_tmp = s.space_x().create_member();
			V_VpV( d_tmp.get(), *Ypy_k, Zpz_k );
			const std::pair<value_type,value_type>
				u_steps = max_near_feas_step( *zero, *d_tmp, dl_iq.get_k(0), du_iq.get_k(0), 0.0 );
			const value_type
				u = std::_MIN( u_steps.first, 1.0 ); // largest positive step size
			alpha_iq.set_k(0) = u;
			if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out	<< "\nFinding u s.t. dl <= u*(Ypy_k+Zpz_k) <= du\n"
					<< "max step length u = " << u << std::endl
					<< "alpha_k = u = " << alpha_iq.get_k(0) << std::endl;
			}
			throw_qp_failure = dual_feasible_point_error();
			break;
		}	
		case QPSolverStats::SUBOPTIMAL_POINT:
		{
			omsg
				<< "\n*** Alert!, the returned QP solution is SUBOPTIMAL_POINT!\n";
			if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
				out << omsg.str();
			}
			throw_qp_failure = true;
			break;
		}
		default:
			assert(0);	// should not happen!
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\n||pz||inf    = " << s.pz().get_k(0).norm_inf()
			<< "\nnu.nz()      = " << s.nu().get_k(0).nz()
			<< "\nmax(|nu(i)|) = " << s.nu().get_k(0).norm_inf()
//			<< "\nmin(|nu(i)|) = " << min_abs( s.nu().get_k(0)() )
			;
		if( m > r ) out << "\n||lambda_k(undecomp)||inf = " << s.lambda().get_k(0).norm_inf();
		out	<< "\n||Zpz||2     = " << s.Zpz().get_k(0).norm_2()
			;
		if(qp_eta > 0.0) out << "\n||Ypy||2 = " << s.Ypy().get_k(0).norm_2();
		out << std::endl;
	}

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
		out << "\npz_k = \n" << s.pz().get_k(0);
		out << "\nnu_k = \n" << s.nu().get_k(0);
		if( m > r )
			out << "\nlambda_k(equ_undecomp) = \n" << *s.lambda().get_k(0).sub_view(equ_undecomp);
		out << "\nZpz_k = \n" << s.Zpz().get_k(0);
		if(qp_eta > 0.0) out << "\nYpy = \n" << s.Ypy().get_k(0);
	}

	if( qp_eta == 1.0 ) {
		omsg
			<< "NullSpaceStepWithInequStd_Step::do_step(...) : Error, a QP relaxation parameter\n"
			<< "of eta = " << qp_eta << " was calculated and therefore it must be assumed\n"
			<< "that the NLP's constraints are infeasible\n"
			<< "Throwing an InfeasibleConstraints exception!\n";
		if(  static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
			out << omsg.str();
		}
		throw InfeasibleConstraints(omsg.str());
	}

	if( throw_qp_failure )
		throw QPFailure( omsg.str(), qp_solver().get_qp_stats() );

	return true;
}

void NullSpaceStepWithInequStd_Step::print_step(
	const Algorithm& algo, poss_type step_poss, GeneralIterationPack::EDoStepType type
	,poss_type assoc_step_poss, std::ostream& out, const std::string& L
	) const
{
return;
	out
		<< L << "*** Calculate the null-space step by solving a constrained QP\n"
		<< L << "qp_grad_k = rGf_k\n"
		<< L << "if w_k is updated then set qp_grad_k = qp_grad_k + zeta_k * w_k\n"
		<< L << "bl = dl_k - Ypy_k\n"
		<< L << "bu = du_k - Ypy_k\n"
		<< L << "if (mI > 0) then\n"
		<< L << "           / h_k - hl  : if h_k < hl\n"
		<< L << "  h_eta =  | 0         : if hl <= h_k <= hu\n"
		<< L << "           \\ h_k - hu  : if h_k > hu\n"
		<< L << "end\n"
		<< L << "etaL = 0.0\n"
		<< L << "*** Determine if we can use simple bounds on pz or not\n"
		<< L << "if m==0 or (Z_k is a variable reduction null space matrix and ||Ypy_k(var_indep)||inf == 0) then\n"
		<< L << "  use_simple_pz_bounds = true\n"
		<< L << "else\n"
		<< L << "  use_simple_pz_bounds = false\n"
		<< L << "end\n"
		<< L << "*** Setup QP arguments\n"
		<< L << "qp_g = qp_grad_k\n"
		<< L << "qp_G = rHL_k\n"
		<< L << "if (mI > 0 ) and (use_simple_pz_bounds == true) then\n"
		<< L << "  qp_dL = bl(var_indep)\n"
		<< L << "  qp_dU = bu(var_indep)\n"
		<< L << "  qp_E  = [ Z_k.D ]\n"
		<< L << "          [ Vz_k  ]\n"
		<< L << "  qp_b  = [ Ypy_k(var_dep)    ]\n"
		<< L << "          [ Vy_k*py_k + h_eta ]\n" 
		<< L << "  qp_eL = [ bl(var_dep)          ]\n"
		<< L << "          [ hl - h_k - Vy_k*py_k ]\n"
		<< L << "  qp_eU = [ bu(var_dep)          ]\n"
		<< L << "          [ hu - h_k - Vy_k*py_k ]\n"
		<< L << "elseif (mI > 0) and (use_simple_pz_bounds == false) then\n"
		<< L << "  qp_dL = -inf\n"
		<< L << "  qp_dU = +inf\n"
		<< L << "  qp_E  = [ Z_k  ]\n"
		<< L << "          [ Vz_k ]\n"
		<< L << "  qp_b  = [ Ypy_k             ]\n"
		<< L << "          [ Vy_k*py_k + h_eta ]\n" 
		<< L << "  qp_eL = [ bl                   ]\n"
		<< L << "          [ hl - h_k - Vy_k*py_k ]\n"
		<< L << "  qp_eU = [ bu                   ]\n"
		<< L << "          [ hu - h_k - Vy_k*py_k ]\n"
		<< L << "elseif (mI == 0) and (use_simple_pz_bounds == true) then\n"
		<< L << "  qp_dL = bl(var_indep),  qp_dU = bu(var_indep)\n"
		<< L << "  if (m > 0) then\n"
		<< L << "    qp_E  = Z_k.D,          qp_b  = Ypy_k(var_dep)\n"
		<< L << "    qp_eL = bl(var_dep),    qp_eU = bu(var_dep)\n"
		<< L << "  end\n"
		<< L << "elseif (mI == 0) and (use_simple_pz_bounds == false) then\n"
		<< L << "  qp_dL = -inf,           qp_dU = +inf\n"
		<< L << "  qp_E  = Z_k,            qp_b  = Ypy_k\n"
		<< L << "  qp_eL = bl,             qp_eU = bu\n"
		<< L << "end\n"
		<< L << "if (m > r) then\n"
		<< L << "  qp_F  = V_k,            qp_f  = Uy_k*py_k + c_k(equ_undecomp)\n"
		<< L << "else\n"
		<< L << "  qp_F  = empty,          qp_f  = empty\n"
		<< L << "end\n"
		<< L << "Given active-set statistics (act_set_stats_km1)\n"
		<< L << "  frac_same = max(num_active-num_adds-num_drops,0)/(num_active)\n"
		<< L << "Use a warm start when frac_same >= warm_start_frac\n"
		<< L << "Solve the following QP to compute qp_d, qp_eta, qp_Ed = qp_E * qp_d\n"
		<< L << ",qp_nu, qp_mu and qp_lambda (" << typeid(qp_solver()).name() << "):\n"
		<< L << "  min  qp_g' * qp_d + 1/2 * qp_d' * qp_G * qp_d + M(eta)\n"
		<< L << "  qp_d <: R^(n-r)\n"
		<< L << "  s.t.\n"
		<< L << "       etaL  <=  eta\n"
		<< L << "       qp_dL <= qp_d                         <= qp_dU   [qp_nu]\n"
		<< L << "       qp_eL <= qp_E * qp_d + (1-eta)*qp_b   <= qp_eU   [qp_mu]\n"
		<< L << "                qp_F * d_qp + (1-eta) * qp_f  = 0       [qp_lambda]\n"
		<< L << "if (qp_testing==QP_TEST) or (fd_deriv_testing==QP_TEST_DEFAULT\n"
		<< L << "and check_results==true) then\n"
		<< L << "  Check the optimality conditions of the above QP\n"
		<< L << "  if the optimality conditions do not check out then\n"
		<< L << "    set throw_qp_failure = true\n"
		<< L << "  end\n"
		<< L << "end\n"
		<< L << "*** Set the solution to the QP subproblem\n"
		<< L << "pz_k  = qp_d\n"
		<< L << "eta_k = qp_eta\n"
		<< L << "if (mI > 0) and (use_simple_pz_bounds == true) then\n"
		<< L << "  nu_k(var_dep)  = qp_mu(1:r),  nu_k(var_indep) = qp_nu\n"
		<< L << "  lambdaI_k      = qp_mu(r+1:r+mI)\n"
		<< L << "  Zpz_k(var_dep) = qp_Ed(1:r), Zpz_k(var_indep) = pz_k\n"
		<< L << "elseif (mI > 0) and (use_simple_pz_bounds == false) then\n"
		<< L << "  nu_k           = qp_mu(1:n)\n"
		<< L << "  lambdaI_k      = qp_mu(n+1:n+mI)\n"
		<< L << "  Zpz_k          = qp_Ed(1:n)\n"
		<< L << "if (mI == 0) and (use_simple_pz_bounds == true) then\n"
		<< L << "  nu_k(var_dep)  = qp_mu,  nu_k(var_indep)  = qp_nu\n"
		<< L << "  Zpz_k(var_dep) = qp_Ed,  Zpz_k(var_indep) = pz_k\n"
		<< L << "elseif (mI == 0) and (use_simple_pz_bounds == false) then\n"
		<< L << "  nu_k  = qp_mu\n"
		<< L << "  Zpz_k = qp_Ed\n"
		<< L << "end\n"
		<< L << "if m > r then\n"
		<< L << "  lambda_k(equ_undecomp) = qp_lambda\n"
		<< L << "end\n"
		<< L << "if (eta_k > 0) then set Ypy_k = (1-eta_k) * Ypy_k\n"
		<< L << "if QP solution is suboptimal then\n"
		<< L << "  throw_qp_failure = true\n"
		<< L << "elseif QP solution is primal feasible (not optimal) then\n"
		<< L << "  throw_qp_failure = primal_feasible_point_error\n"
		<< L << "elseif QP solution is dual feasible (not optimal) then\n"
		<< L << "  find max u s.t.\n"
		<< L << "    dl_k <= u*(Ypy_k+Zpz_k) <= du_k\n"
		<< L << "  alpha_k = u\n"
		<< L << "  throw_qp_failure = dual_feasible_point_error\n"
		<< L << "end\n"
		<< L << "if (eta_k == 1.0) then\n"
		<< L << "  The constraints are infeasible!\n"
		<< L << "  throw InfeasibleConstraints(...)\n"
		<< L << "end\n"
		<< L << "if (throw_qp_failure == true) then\n"
		<< L << "  throw QPFailure(...)\n"
		<< L << "end\n"
		;
}

}	// end namespace ReducedSpaceSQPPack
