// /////////////////////////////////////////////////////////////////////////
// NLPAlgoConfigIP.cpp
//
// Copyright (C) 2001
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

#include <sstream>
#include <typeinfo>
#include <iostream>

#include "debug.hpp"

#include "NLPAlgoConfigIP.hpp"
#include "NLPInterfacePack/src/abstract/tools/NLPBarrier.hpp"
#include "MoochoPack/src/NLPAlgo.hpp"
#include "MoochoPack/src/IpState.hpp"
#include "MoochoPack/src/NLPAlgoContainer.hpp"
#include "AbstractLinAlgPack/src/serial/implementations/MatrixSymPosDefCholFactor.hpp"                     // rHL 
//#include "ConstrainedOptPack/src/matrices/MatrixSymPosDefInvCholFactor.hpp"		// .
#include "ConstrainedOptPack/src/matrices/MatrixSymPosDefLBFGS.hpp"				// .
//#include "ConstrainedOptPack/src/matrices/MatrixHessianSuperBasicInitDiagonal.hpp"// | rHL (super basics)
#include "AbstractLinAlgPack/src/abstract/tools/MatrixSymDiagStd.hpp"                          // |

#include "NLPInterfacePack/src/abstract/interfaces/NLPDirect.hpp"
#include "NLPInterfacePack/src/abstract/interfaces/NLPVarReductPerm.hpp"
#include "NLPInterfacePack/src/abstract/tools/CalcFiniteDiffProd.hpp"

// line search
#include "ConstrainedOptPack/src/globalization/DirectLineSearchArmQuad_Strategy.hpp"
#include "ConstrainedOptPack/src/globalization/DirectLineSearchArmQuad_StrategySetOptions.hpp"
#include "ConstrainedOptPack/src/globalization/MeritFuncNLPL1.hpp"
#include "ConstrainedOptPack/src/globalization/MeritFuncNLPModL1.hpp"

// Basis permutations and direct sparse solvers
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
#include "ConstrainedOptPack/src/decompositions/DecompositionSystemVarReductPerm.hpp"
#endif

#include "MoochoPack/src/std/MoochoAlgorithmStepNames.hpp"

#include "MoochoPack/src/std/UpdateBarrierParameter_Step.hpp"

#include "MoochoPack/src/std/PreEvalNewPointBarrier_Step.hpp"
#include "MoochoPack/src/std/PostEvalNewPointBarrier_Step.hpp"
#include "MoochoPack/src/std/ReducedGradientStd_Step.hpp"
//#include "MoochoPack/src/std/InitFinDiffReducedHessian_Step.hpp"
//#include "MoochoPack/src/std/InitFinDiffReducedHessian_StepSetOptions.hpp"
#include "MoochoPack/src/std/ReducedHessianSecantUpdateStd_Step.hpp"
#include "MoochoPack/src/std/ReducedHessianSecantUpdateBFGSFull_Strategy.hpp"
//#include "MoochoPack/src/std/ReducedHessianSecantUpdateBFGSProjected_Strategy.hpp"
//#include "MoochoPack/src/std/ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions.hpp"
//#include "MoochoPack/src/std/ReducedHessianSecantUpdateLPBFGS_Strategy.hpp"
//#include "MoochoPack/src/std/ReducedHessianSecantUpdateLPBFGS_StrategySetOptions.hpp"
#include "MoochoPack/src/std/BFGSUpdate_Strategy.hpp"
#include "MoochoPack/src/std/BFGSUpdate_StrategySetOptions.hpp"
#include "MoochoPack/src/std/QuasiNormalStepStd_Step.hpp"
#include "MoochoPack/src/std/CheckDescentRangeSpaceStep_Step.hpp"
#include "MoochoPack/src/std/CheckDecompositionFromPy_Step.hpp"
#include "MoochoPack/src/std/CheckDecompositionFromRPy_Step.hpp"
#include "MoochoPack/src/std/TangentialStepIP_Step.hpp"
//#include "MoochoPack/src/std/TangentialStepWithoutBounds_Step.hpp"
//#include "MoochoPack/src/std/TangentialStepWithInequStd_Step.hpp"
//#include "MoochoPack/src/std/TangentialStepWithInequStd_StepSetOptions.hpp"
//#include "MoochoPack/src/std/SetDBoundsStd_AddedStep.hpp"
#include "MoochoPack/src/std/QPFailureReinitReducedHessian_Step.hpp"
#include "MoochoPack/src/std/CalcDFromYPYZPZ_Step.hpp"
#include "MoochoPack/src/std/CalcD_vStep_Step.hpp"

#include "MoochoPack/src/std/PreProcessBarrierLineSearch_Step.hpp"
#include "MoochoPack/src/std/PostProcessBarrierLineSearch_Step.hpp"
#include "MoochoPack/src/std/LineSearchFailureNewDecompositionSelection_Step.hpp"
#include "MoochoPack/src/std/LineSearchFilter_Step.hpp"
#include "MoochoPack/src/std/LineSearchFilter_StepSetOptions.hpp"
#include "MoochoPack/src/std/LineSearchFullStep_Step.hpp"
#include "MoochoPack/src/std/LineSearchDirect_Step.hpp"
//#include "MoochoPack/src/std/LineSearch2ndOrderCorrect_Step.hpp"
//#include "MoochoPack/src/std/LineSearch2ndOrderCorrect_StepSetOptions.hpp"
//#include "MoochoPack/src/std/FeasibilityStepReducedStd_Strategy.hpp"
//#include "MoochoPack/src/std/FeasibilityStepReducedStd_StrategySetOptions.hpp"
//#include "MoochoPack/src/std/QuasiRangeSpaceStepStd_Strategy.hpp"
//#include "MoochoPack/src/std/QuasiRangeSpaceStepTailoredApproach_Strategy.hpp"
//#include "MoochoPack/src/std/LineSearchWatchDog_Step.hpp"
//#include "MoochoPack/src/std/LineSearchWatchDog_StepSetOptions.hpp"
//#include "MoochoPack/src/std/LineSearchFullStepAfterKIter_Step.hpp"
//#include "MoochoPack/src/std/CalcLambdaIndepStd_AddedStep.hpp"
#include "MoochoPack/src/std/CalcReducedGradLagrangianStd_AddedStep.hpp"
#include "MoochoPack/src/std/CheckConvergenceStd_AddedStep.hpp"
#include "MoochoPack/src/std/CheckConvergenceIP_Strategy.hpp"
#include "MoochoPack/src/std/CheckSkipBFGSUpdateStd_StepSetOptions.hpp"
#include "MoochoPack/src/std/MeritFunc_PenaltyParamUpdate_AddedStepSetOptions.hpp"
#include "MoochoPack/src/std/MeritFunc_PenaltyParamUpdateMultFree_AddedStep.hpp"
//#include "MoochoPack/src/std/MeritFunc_PenaltyParamUpdateWithMult_AddedStep.hpp"
//#include "MoochoPack/src/std/MeritFunc_PenaltyParamsUpdateWithMult_AddedStep.hpp"
//#include "MoochoPack/src/std/MeritFunc_ModifiedL1LargerSteps_AddedStep.hpp"
//#include "MoochoPack/src/std/MeritFunc_ModifiedL1LargerSteps_AddedStepSetOptions.hpp"
//#include "MoochoPack/src/std/ActSetStats_AddedStep.hpp"
//#include "MoochoPack/src/std/NumFixedDepIndep_AddedStep.hpp"
#include "MoochoPack/src/std/UpdateReducedSigma_Step.hpp"

#include "MoochoPack/src/std/quasi_newton_stats.hpp"

// Misc utilities
#include "AbstractFactoryStd.hpp"
#include "dynamic_cast_verbose.hpp"
#include "ReleaseResource_ref_count_ptr.hpp"
#include "ThrowException.hpp"

// Stuff to read in options
#include "StringToIntMap.hpp"
#include "StringToBool.hpp"

// Stuff for exact reduced hessian
//#include "MoochoPack/src/std/ReducedHessianExactStd_Step.hpp"
//#include "MoochoPack/src/std/CrossTermExactStd_Step.hpp"
//#include "MoochoPack/src/std/DampenCrossTermStd_Step.hpp"

namespace {
	const double INF_BASIS_COND_CHANGE_FRAC      = 1e+20;
}

namespace MoochoPack {

//
// Here is where we define the default values for the algorithm.  These
// should agree with what are in the Moocho.opt.NLPAlgoConfigIP file.
//
NLPAlgoConfigIP::SOptionValues::SOptionValues()
	:max_basis_cond_change_frac_(-1.0)
	,exact_reduced_hessian_(false)
	,quasi_newton_(QN_AUTO)
	,num_lbfgs_updates_stored_(-1)
	,lbfgs_auto_scaling_(true)
	,hessian_initialization_(INIT_HESS_AUTO)
	,qp_solver_type_(QP_AUTO)
	,reinit_hessian_on_qp_fail_(true)
	,line_search_method_(LINE_SEARCH_AUTO)
	,merit_function_type_(MERIT_FUNC_AUTO)
	,l1_penalty_param_update_(L1_PENALTY_PARAM_AUTO)
	,full_steps_after_k_(-1)
{}

NLPAlgoConfigIP::NLPAlgoConfigIP()
{}

NLPAlgoConfigIP::~NLPAlgoConfigIP()
{}

// overridden from NLPAlgoConfig

void NLPAlgoConfigIP::set_options( const options_ptr_t& options )
{
	options_ = options;
	decomp_sys_step_builder_.set_options(options);
}

const NLPAlgoConfig::options_ptr_t&
NLPAlgoConfigIP::get_options() const
{
	return options_;
}

void NLPAlgoConfigIP::config_algo_cntr(
	NLPAlgoContainer   *algo_cntr
	,std::ostream       *trase_out
	)
{
	namespace afp = MemMngPack;
	namespace mmp = MemMngPack;
	using mmp::ref_count_ptr;
	using DynamicCastHelperPack::dyn_cast;
	using DynamicCastHelperPack::dyn_cast;

	if(trase_out) {
		*trase_out
			<< std::endl
			<< "*****************************************************************\n"
			<< "*** NLPAlgoConfigIP configuration                               ***\n"
			<< "***                                                           ***\n"
			<< "*** Here, summary information about how the algorithm is      ***\n"
			<< "*** configured is printed so that the user can see how the    ***\n"
			<< "*** properties of the NLP and the set options influence       ***\n"
			<< "*** how an algorithm is configured.                           ***\n"
			<< "*****************************************************************\n";
	}

	// ////////////////////////////////////////////////////////////
	// A. ???

	// /////////////////////////////////////////////////////////////////////////
	// B. Create an algo object, give to algo_cntr, then give algo_cntr to algo

	if(trase_out)
		*trase_out << "\n*** Creating the NLPAlgo algo object ...\n";

	typedef mmp::ref_count_ptr<NLPAlgo>	algo_ptr_t;
	algo_ptr_t algo = mmp::rcp(new NLPAlgo);
	assert(algo.get());
	algo_cntr->set_algo(algo);
	algo->set_algo_cntr(algo_cntr);

	// /////////////////////////////////////////////
	// C. Configure algo

	// /////////////////////////////////////////////////////
	// C.0 Set the nlp and track objects

	if(trase_out)
		*trase_out << "\n*** Setting the NLP and track objects to the algo object ...\n";

	algo->set_nlp( algo_cntr->get_nlp().get() );
	algo->set_track( algo_cntr->get_track() );

	// ////////////////////////////////////////////////
	// Determine what the options are:

	// Readin the options
	if(options_.get()) {
		readin_options( *options_, &uov_, trase_out );
	}
	else {
		if(trase_out) {
			*trase_out
				<< "\n*** Warning, no OptionsFromStream object was set so a default set"
					" of options will be used!\n";
		}
	}	

	NLP &nlp = algo->nlp();
	nlp.initialize(algo->algo_cntr().check_results());
	// Get the dimensions of the NLP
	const size_type
		n   = nlp.n(),
		m   = nlp.m(),
		mI  = nlp.mI(),
		r   = m, // ToDo: Compute this for real!
		dof = n - r,
		nb  = nlp.num_bounded_x();

	// Process the NLP
	NLPFirstOrder    *nlp_foi = NULL;
	NLPSecondOrder   *nlp_soi = NULL;
	NLPDirect  *nlp_fod = NULL;
	bool                 tailored_approach = false;
	decomp_sys_step_builder_.process_nlp_and_options(
		trase_out, nlp
		,&nlp_foi, &nlp_soi, &nlp_fod, &tailored_approach
		);

	const int max_dof_quasi_newton_dense
		= decomp_sys_step_builder_.current_option_values().max_dof_quasi_newton_dense_;

	// Make sure that we can handle this type of NLP currently
	THROW_EXCEPTION(
		m == 0, std::logic_error
		,"NLPAlgoConfigIP::config_algo_cntr(...) : Error, "
		"can not currently solve an unconstrained NLP!" );
	THROW_EXCEPTION(
		n == m, std::logic_error
		,"NLPAlgoConfigIP::config_algo_cntr(...) : Error, "
		"can not currently solve a square system of equations!" );
	THROW_EXCEPTION(
		mI, std::logic_error
		,"NLPAlgoConfigIP::config_algo_cntr(...) : Error, "
		"can not currently solve an NLPcd .. with general inequalities!" );

	// //////////////////////////////////////////////////////
	// C.1.  Sort out the options

	if(trase_out)
		*trase_out
			<< "\n*** Sorting out some of the options given input options ...\n";

	if( tailored_approach ) {
		// Change the options for the tailored approach. 
		if(trase_out) {
			*trase_out
				<< "\nThis is a tailored approach NLP (NLPDirect) which forces the following options:\n"
				<< "merit_function_type         = L1;\n"
				<< "l1_penalty_parameter_update = MULT_FREE;\n"
				<< "null_space_matrix           = EXPLICIT;\n"
				;
		}
		cov_.merit_function_type_
			= MERIT_FUNC_L1;
		cov_.l1_penalty_param_update_
			= L1_PENALTY_PARAM_MULT_FREE;
		decomp_sys_step_builder_.current_option_values().null_space_matrix_type_
			= DecompositionSystemStateStepBuilderStd::NULL_SPACE_MATRIX_EXPLICIT;
	}

	if( !tailored_approach && uov_.merit_function_type_ != MERIT_FUNC_L1  ) {
		if(trase_out) {
			*trase_out
				<< "\nThe only merit function currently supported is:\n"
				<< "merit_function_type         = L1;\n"
				;
		}
		cov_.merit_function_type_		= MERIT_FUNC_L1;
	}

	// Decide what type of quasi newton update to use
	switch( uov_.quasi_newton_ ) {
		case QN_AUTO: {
			if(trase_out)
				*trase_out
					<< "\nquasi_newton == AUTO:"
					<< "\nnlp.num_bounded_x() == " << nlp.num_bounded_x() << ":\n";
			//if( n - r > cov_.max_dof_quasi_newton_dense_ ) {
			//	if(trase_out)
			//		*trase_out
			//			<< "n-r = " << n-r << " > max_dof_quasi_newton_dense = "
			//			<< cov_.max_dof_quasi_newton_dense_ <<  ":\n"
			//			<< "setting quasi_newton == LBFGS\n";
			//	cov_.quasi_newton_ = QN_LBFGS;
			//}
			//else {
				if(trase_out)
					*trase_out
						<< "n-r = " << n-r << " <= max_dof_quasi_newton_dense = "
						<< max_dof_quasi_newton_dense << ":\n"
						<< "setting quasi_newton == BFGS\n";
				cov_.quasi_newton_ = QN_BFGS;
				//}
			break;
		}
		case QN_BFGS:
		case QN_PBFGS:
		case QN_LBFGS:
		case QN_LPBFGS:
			cov_.quasi_newton_ = uov_.quasi_newton_;
			break;
	    default:
			assert(0); // Invalid option!
	}

	// ToDo: Sort out the rest of the options!

	// Set the default options that where not already set yet
	set_default_options(uov_,&cov_,trase_out);

	// ToDo: Implement the 2nd order correction linesearch
	if( cov_.line_search_method_ == LINE_SEARCH_2ND_ORDER_CORRECT ) {
		if(trase_out)
			*trase_out <<
				"\nline_search_method == 2ND_ORDER_CORRECT:\n"
				"Sorry, the second order corrrection linesearch is not updated yet!\n"
				"setting line_search_method = DIRECT ...\n";
		cov_.line_search_method_ = LINE_SEARCH_DIRECT;
	}
	if( cov_.line_search_method_ == LINE_SEARCH_WATCHDOG ) {
		if(trase_out)
			*trase_out <<
				"\nline_search_method ==WATCHDOG:\n"
				"Sorry, the watchdog linesearch is not updated yet!\n"
				"setting line_search_method = DIRECT ...\n";
		cov_.line_search_method_ = LINE_SEARCH_DIRECT;
	}
	
	// /////////////////////////////////////////////////////
	// C.1. Create the decomposition system object

	typedef ref_count_ptr<DecompositionSystem> decomp_sys_ptr_t;
	decomp_sys_ptr_t decomp_sys;
	decomp_sys_step_builder_.create_decomp_sys(
		trase_out, nlp, nlp_foi, nlp_soi, nlp_fod, tailored_approach
		,&decomp_sys
		);

#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
	ref_count_ptr<DecompositionSystemVarReductPerm>
		decomp_sys_perm = mmp::rcp_dynamic_cast<DecompositionSystemVarReductPerm>(decomp_sys);
#endif

	// /////////////////////////////////////////////////////
	// C.2. Create and set the state object

	if(trase_out)
		*trase_out
			<< "\n*** Creating the state object and setting up iteration quantity objects ...\n";

	{
		//
		// Create the state object with the vector spaces
		//

		typedef ref_count_ptr<IpState>   state_ptr_t;
		state_ptr_t
			state = mmp::rcp(
				new IpState(
					decomp_sys
					,nlp.space_x()
					,nlp.space_c()
					,nlp.space_h()
					,( tailored_approach
					   ? ( nlp_fod->var_dep().size() 
						   ? nlp.space_x()->sub_space(nlp_fod->var_dep())->clone()
						   : mmp::null )
					   : decomp_sys->space_range() // could be NULL for BasisSystemPerm
						)
					,( tailored_approach
					   ?( nlp_fod->var_indep().size()
						  ? nlp.space_x()->sub_space(nlp_fod->var_indep())->clone()
						  : mmp::null )
					   : decomp_sys->space_null() // could be NULL for BasisSystemPerm
						)
					)
				);

		//
		// Set the iteration quantities for the NLP matrix objects
		//

		decomp_sys_step_builder_.add_iter_quantities(
			trase_out, nlp, nlp_foi, nlp_soi, nlp_fod, tailored_approach, decomp_sys
			,state
			);

		///*****************************************************
		// Set the iteration quantities for the barrier terms
		///*****************************************************
		state->set_iter_quant(
		  Vu_name
		  ,mmp::rcp(
			new IterQuantityAccessContiguous<MatrixSymDiagStd>(
			  1,
			  Vu_name,
			  mmp::rcp( new afp::AbstractFactoryStd<MatrixSymDiagStd,MatrixSymDiagStd,
						MatrixSymDiagStd::PostMod>( nlp.space_x() ) 
				)
			  )
			)
		  );

		state->set_iter_quant(
		  Vl_name
		  ,mmp::rcp(
			new IterQuantityAccessContiguous<MatrixSymDiagStd>(
			  1,
			  Vl_name,
			  mmp::rcp( new afp::AbstractFactoryStd<MatrixSymDiagStd,MatrixSymDiagStd,
						MatrixSymDiagStd::PostMod>( nlp.space_x() ) 
				)
			  )
			)
		  );

		state->set_iter_quant(
		  invXu_name
		  ,mmp::rcp(
			new IterQuantityAccessContiguous<MatrixSymDiagStd>(
			  1,
			  invXu_name,
			  mmp::rcp( new afp::AbstractFactoryStd<MatrixSymDiagStd,MatrixSymDiagStd,
						MatrixSymDiagStd::PostMod>( nlp.space_x() ) 
				)
			  )
			)
		  );

		state->set_iter_quant(
		  invXl_name
		  ,mmp::rcp(
			new IterQuantityAccessContiguous<MatrixSymDiagStd>(
			  1,
			  invXl_name,
			  mmp::rcp( new afp::AbstractFactoryStd<MatrixSymDiagStd,MatrixSymDiagStd,
						MatrixSymDiagStd::PostMod>( nlp.space_x() ) 
				)
			  )
			)
		  );

		state->set_iter_quant(
		  rHB_name
		  ,mmp::rcp(
			new IterQuantityAccessContiguous<MatrixSymOp>(
			  1,
			  rHB_name,
			  mmp::rcp(
				new afp::AbstractFactoryStd<MatrixSymOp,MatrixSymPosDefCholFactor,MatrixSymPosDefCholFactor::PostMod>(
				  MatrixSymPosDefCholFactor::PostMod(
					true      // maintain_original
					,false    // maintain_factor
					,true     // allow_factor      (always!)
					)
				  )
				)
			  )
			)
		  );
				
		state->set_iter_quant(
		  B_name
		  ,mmp::rcp(
			new IterQuantityAccessContiguous<MatrixSymOp>(
			  1,
			  B_name,
			  mmp::rcp(
				new afp::AbstractFactoryStd<MatrixSymOp,MatrixSymPosDefCholFactor,MatrixSymPosDefCholFactor::PostMod>(
				  MatrixSymPosDefCholFactor::PostMod(
					true      // maintain_original
					,false    // maintain_factor
					,true     // allow_factor      (always!)
					)
				  )
				)
			  )
			)
		  );

		state->set_iter_quant(
		  Sigma_name
		  ,mmp::rcp(
			new IterQuantityAccessContiguous<MatrixSymDiagStd>(
			  1,
			  Sigma_name,
			  mmp::rcp( new afp::AbstractFactoryStd<MatrixSymDiagStd,MatrixSymDiagStd,
						MatrixSymDiagStd::PostMod>( nlp.space_x() ) 
				)
			  )
			)
		  );

		// These iteration quantities are defined in IpState, 		
		// force their creation and resize them
		dyn_cast< IterQuantityAccessContiguous<value_type> >(state->barrier_obj()).resize(2);
		dyn_cast< IterQuantityAccessContiguous<VectorMutable> >(state->grad_barrier_obj()).resize(2);

		// Add reduced Hessian of the Lagrangian

		if( !cov_.exact_reduced_hessian_ ) {
			ref_count_ptr<afp::AbstractFactory<MatrixSymOp> >
				abstract_factory_rHL = mmp::rcp(
					new afp::AbstractFactoryStd<MatrixSymOp,MatrixSymPosDefCholFactor,MatrixSymPosDefCholFactor::PostMod>(
						MatrixSymPosDefCholFactor::PostMod(
							true    // maintain_original
							,false  // maintain_factor
							,true   // allow_factor      (always!)
							)
						)
					);
			state->set_iter_quant(
				rHL_name
				,mmp::rcp(
					new IterQuantityAccessContiguous<MatrixSymOp>(
						1
						,rHL_name
						,abstract_factory_rHL
						)
					)
				);
		}
		else {
			assert(0); // ToDo: Add rHL for an exact reduced Hessian!
		}
		
		//
		// Set the NLP merit function 
		//

		if( cov_.line_search_method_ != LINE_SEARCH_NONE 
			 && cov_.line_search_method_ != LINE_SEARCH_FILTER) {
			ref_count_ptr<afp::AbstractFactory<MeritFuncNLP> >
				merit_func_factory = mmp::null;
			switch( cov_.merit_function_type_ ) {
				case MERIT_FUNC_L1:
					merit_func_factory = mmp::rcp(
						new afp::AbstractFactoryStd<MeritFuncNLP,MeritFuncNLPL1>());
					break;
				case MERIT_FUNC_MOD_L1:
				case MERIT_FUNC_MOD_L1_INCR:
					merit_func_factory = mmp::rcp(
						new afp::AbstractFactoryStd<MeritFuncNLP,MeritFuncNLPModL1>());
					break;
				default:
					assert(0);	// local programming error
			}
			state->set_iter_quant(
				merit_func_nlp_name
				,mmp::rcp(
					new IterQuantityAccessContiguous<MeritFuncNLP>(
						1
						,merit_func_nlp_name
						,merit_func_factory
						)
					)
				);
		}

		if (cov_.line_search_method_ == LINE_SEARCH_FILTER)
		  {
		    // Add the filter iteration quantity
		    state->set_iter_quant(
				FILTER_IQ_STRING
				,mmp::rcp(
				     new IterQuantityAccessContiguous<Filter_T>(1,FILTER_IQ_STRING)
				)
		    );
		  }

		//
		// Resize the number of storage locations (these can be changed later).
		//
		// Also, touch all of the value_type, index_type and vector iteration quantities
		// that we know about so that when state.dump_iter_quant() is called, all of the
		// iteration quantities will be included.
		//
		
		typedef IterQuantityAccessContiguous<value_type>            IQ_scalar_cngs;
		typedef IterQuantityAccessContiguous<VectorMutable>   IQ_vector_cngs;

		dyn_cast<IQ_vector_cngs>(state->x()).resize(2);
		dyn_cast<IQ_scalar_cngs>(state->f()).resize(2);
		if(m) dyn_cast<IQ_vector_cngs>(state->c()).resize(2);
		dyn_cast<IQ_vector_cngs>(state->Gf()).resize(2);
		if(m && nlp_foi) state->Gc();
		if(mI) state->Gh();

		if( m
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
		   && decomp_sys_perm.get() == NULL
#endif
			) state->py();
		if(m) dyn_cast<IQ_vector_cngs>(state->Ypy()).resize(2);
		if( m
#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
			&& decomp_sys_perm.get() == NULL
#endif
			) state->pz();
		if(m) dyn_cast<IQ_vector_cngs>(state->Zpz()).resize(2);
		dyn_cast<IQ_vector_cngs>(state->d()).resize(2);

		if( n > m ) {
			dyn_cast<IQ_vector_cngs>(state->rGf()).resize(2);
			state->w();
			state->zeta();
			state->qp_grad();
		}
		state->eta();
		
		dyn_cast<IQ_scalar_cngs>(state->barrier_parameter()).resize(2);

		dyn_cast<IQ_scalar_cngs>(state->alpha()).resize(2);
		dyn_cast<IQ_scalar_cngs>(state->mu()).resize(2);
		dyn_cast<IQ_scalar_cngs>(state->phi()).resize(2);

		dyn_cast<IQ_scalar_cngs>(state->opt_kkt_err()).resize(2);
		dyn_cast<IQ_scalar_cngs>(state->feas_kkt_err()).resize(2);
		if( n > m ) {
			dyn_cast<IQ_vector_cngs>(state->rGL()).resize(2);
		}
		if(m) dyn_cast<IQ_vector_cngs>(state->lambda()).resize(2);
		if(mI) dyn_cast<IQ_vector_cngs>(state->lambdaI()).resize(2);
		dyn_cast<IQ_vector_cngs>(state->nu()).resize(2);

		// Set the state object
		algo->set_state( state );
	}

	// /////////////////////////////////////////////////////
	// C.3  Create and set the step objects

	if(trase_out)
		*trase_out << "\n*** Creating and setting the step objects ...\n";

	{

		//
		// Create some standard step objects that will be used by many different
		// specific algorithms
		//
		
		typedef ref_count_ptr<AlgorithmStep>   algo_step_ptr_t;

		// Create the EvalNewPoint step and associated objects
		algo_step_ptr_t                                    eval_new_point_step           = mmp::null;
		ref_count_ptr<CalcFiniteDiffProd>                  calc_fd_prod                  = mmp::null;
		ref_count_ptr<VariableBoundsTester>                bounds_tester                 = mmp::null;
		ref_count_ptr<NewDecompositionSelection_Strategy>  new_decomp_selection_strategy = mmp::null;
		decomp_sys_step_builder_.create_eval_new_point(
			trase_out, nlp, nlp_foi, nlp_soi, nlp_fod, tailored_approach, decomp_sys
			,&eval_new_point_step, &calc_fd_prod, &bounds_tester, &new_decomp_selection_strategy
			);

		// UpdateBarrierParameter_Step
		mmp::ref_count_ptr<UpdateBarrierParameter_Step>	updateBarrierParameter_step  = mmp::null;
		updateBarrierParameter_step = mmp::rcp(new UpdateBarrierParameter_Step());
		if(options_.get()) 
			{
			UpdateBarrierParameter_StepSetOptions options_setter(updateBarrierParameter_step.get());
			options_setter.set_options(*options_);
			}
		
		// PreEvalNewPointBarrier_Step
		mmp::ref_count_ptr<PreEvalNewPointBarrier_Step>  preEvalNewPointBarrier_step  = mmp::null;
		preEvalNewPointBarrier_step = mmp::rcp(new PreEvalNewPointBarrier_Step());
		if(options_.get()) 
			{
			PreEvalNewPointBarrier_StepSetOptions
				options_setter(preEvalNewPointBarrier_step.get());
			options_setter.set_options(*options_);
			}

		// PostEvalNewPointBarrier_Step
		algo_step_ptr_t postEvalNewPointBarrier_step = mmp::rcp(new PostEvalNewPointBarrier_Step());

		// ReducedGradient_Step
		algo_step_ptr_t    reduced_gradient_step = mmp::null;
		if( !tailored_approach ) {
			reduced_gradient_step = mmp::rcp(new ReducedGradientStd_Step());
		}

		// RangeSpace_Step
		algo_step_ptr_t    range_space_step_step = mmp::null;
		if( !tailored_approach ) {
			range_space_step_step = mmp::rcp(new QuasiNormalStepStd_Step());
		}

		// Check and change decomposition 
		algo_step_ptr_t    check_decomp_from_py_step  = mmp::null;
		algo_step_ptr_t    check_decomp_from_Rpy_step = mmp::null;
		if( new_decomp_selection_strategy.get() && cov_.max_basis_cond_change_frac_ < INF_BASIS_COND_CHANGE_FRAC ) {
			check_decomp_from_py_step = mmp::rcp(
				new CheckDecompositionFromPy_Step(
					new_decomp_selection_strategy
					,cov_.max_basis_cond_change_frac_
					) );
			check_decomp_from_Rpy_step = mmp::rcp(
				new CheckDecompositionFromRPy_Step(
					new_decomp_selection_strategy
					,cov_.max_basis_cond_change_frac_
					) );
		}

		// CheckDescentRangeSpaceStep
		algo_step_ptr_t    check_descent_range_space_step_step = mmp::null;
		if( algo->algo_cntr().check_results() ) {
			check_descent_range_space_step_step = mmp::rcp(new CheckDescentRangeSpaceStep_Step(calc_fd_prod));
		}

		// ReducedGradient_Step
		//algo_step_ptr_t    reduced_gradient_step = mmp::null;
		//if( !tailored_approach ) {
		//	reduced_gradient_step = mmp::rcp(new ReducedGradientStd_Step());
		//}

		// CheckSkipBFGSUpdate
		algo_step_ptr_t    check_skip_bfgs_update_step = mmp::null;
		if(!cov_.exact_reduced_hessian_) {
			ref_count_ptr<CheckSkipBFGSUpdateStd_Step>
				step = mmp::rcp(new CheckSkipBFGSUpdateStd_Step());
			if(options_.get()) {
				CheckSkipBFGSUpdateStd_StepSetOptions
					opt_setter( step.get() );
				opt_setter.set_options( *options_ );
			}
			check_skip_bfgs_update_step = step;
		}

		// ReducedHessian_Step
		algo_step_ptr_t    reduced_hessian_step = mmp::null;
		{
			// Get the strategy object that will perform the actual secant update.
			ref_count_ptr<ReducedHessianSecantUpdate_Strategy>
				secant_update_strategy = mmp::null;
			switch( cov_.quasi_newton_ )
			{
			    case QN_BFGS:
			    case QN_PBFGS:
			    case QN_LBFGS:
			    case QN_LPBFGS:
				{
					// create and setup the actual BFGS strategy object
					typedef ref_count_ptr<BFGSUpdate_Strategy> bfgs_strategy_ptr_t;
					bfgs_strategy_ptr_t
						bfgs_strategy = mmp::rcp(new BFGSUpdate_Strategy);
					if(options_.get()) { 
						BFGSUpdate_StrategySetOptions
							opt_setter( bfgs_strategy.get() );
						opt_setter.set_options( *options_ );
					}
					switch( cov_.quasi_newton_ ) {
					    case QN_BFGS:
					    case QN_LBFGS:
						{
							secant_update_strategy = mmp::rcp(new ReducedHessianSecantUpdateBFGSFull_Strategy(bfgs_strategy));
							break;
						}
					    case QN_PBFGS:
					    case QN_LPBFGS:
						{
							THROW_EXCEPTION(
								true, std::logic_error
								,"NLPAlgoConfigIP::config_algo_cntr(...) : Error, "
								"The quansi_newton options of PBFGS and LPBFGS have not been updated yet!" );
							break;
						}
					}
					break;
				}
			    default:
					assert(0);
			}
			
			// Finally build the step object
			reduced_hessian_step = mmp::rcp(
				new ReducedHessianSecantUpdateStd_Step( secant_update_strategy ) );
			// Add the QuasiNewtonStats iteration quantity
			algo->state().set_iter_quant(
				quasi_newton_stats_name
				,mmp::rcp(new IterQuantityAccessContiguous<QuasiNewtonStats>(
					1
					,quasi_newton_stats_name
#ifdef _MIPS_CXX
					,mmp::ref_count_ptr<MemMngPack::AbstractFactoryStd<QuasiNewtonStats,QuasiNewtonStats> >(
						new MemMngPack::AbstractFactoryStd<QuasiNewtonStats,QuasiNewtonStats>())
#endif
					)
				));
		}

		// UpdateReducedSigma_Step
		MemMngPack::ref_count_ptr<UpdateReducedSigma_Step> updateReducedSigma_step = mmp::null;
		updateReducedSigma_step = mmp::rcp(new UpdateReducedSigma_Step());
		if(options_.get()) 
			{
			UpdateReducedSigma_StepSetOptions
				options_setter(updateReducedSigma_step.get());
			options_setter.set_options(*options_);
			}
		
		// NullSpace_Step
		algo_step_ptr_t    null_space_step = mmp::rcp(new TangentialStepIP_Step());
		/*		algo_step_ptr_t    set_d_bounds_step    = mmp::null;
		algo_step_ptr_t    null_space_step_step = mmp::null;
		if( mI == 0 && nb == 0 ) {
			null_space_step_step = mmp::rcp(new TangentialStepWithoutBounds_Step());
		}
		else {
			// Step object that sets bounds for QP subproblem
			set_d_bounds_step = mmp::rcp(new SetDBoundsStd_AddedStep());
			// QP Solver object
			mmp::ref_count_ptr<QPSolverRelaxed>  qp_solver = mmp::null;
				// ToDo: Create the QP solver!
			// QP solver tester
			mmp::ref_count_ptr<QPSolverRelaxedTester> 
				qp_solver_tester = mmp::rcp(new QPSolverRelaxedTester());
			if(options_.get()) {
				QPSolverRelaxedTesterSetOptions
					opt_setter( qp_solver_tester.get() );
				opt_setter.set_options( *options_ );
			}
			// The null-space step
			mmp::ref_count_ptr<TangentialStepWithInequStd_Step>
				null_space_step_with_inequ_step = mmp::rcp(
					new TangentialStepWithInequStd_Step(
						qp_solver, qp_solver_tester ) );
			if(options_.get()) {
				TangentialStepWithInequStd_StepSetOptions
					opt_setter( null_space_step_with_inequ_step.get() );
				opt_setter.set_options( *options_ );
			}
			null_space_step_step = null_space_step_with_inequ_step;
			// Step for reinitialization reduced Hessian on QP failure
			null_space_step_step = mmp::rcp(
				new QPFailureReinitReducedHessian_Step(null_space_step_step)
				);
				}*/

		// CalcDFromYPYZPZ_Step
		algo_step_ptr_t calc_d_from_Ypy_Zpy_step = mmp::null;
			{
			calc_d_from_Ypy_Zpy_step = mmp::rcp(new CalcDFromYPYZPZ_Step());
			}

		// CalcD_vStep_Step
		algo_step_ptr_t calc_d_v_step_step = mmp::rcp(new CalcD_vStep_Step());

		// build the barrier nlp decorator to be used by the line search
		MemMngPack::ref_count_ptr<NLPInterfacePack::NLPBarrier> barrier_nlp = mmp::rcp(new NLPInterfacePack::NLPBarrier());
		barrier_nlp->InitializeFromNLP( algo_cntr->get_nlp() );

		// PreProcessBarrierLineSearch_Step
		algo_step_ptr_t preprocess_barrier_linesearch_step = mmp::rcp(new PreProcessBarrierLineSearch_Step(barrier_nlp));

		// PostProcessBarrierLineSearch_Step
		algo_step_ptr_t postprocess_barrier_linesearch_step = mmp::rcp(new PostProcessBarrierLineSearch_Step(barrier_nlp));

		// CalcReducedGradLagrangianStd_AddedStep
		algo_step_ptr_t    calc_reduced_grad_lagr_step = mmp::null;
		{
			calc_reduced_grad_lagr_step = mmp::rcp(
				new CalcReducedGradLagrangianStd_AddedStep() );
		}

		// CheckConvergence_Step
		algo_step_ptr_t    check_convergence_step = mmp::null;
			{
			// Create the strategy object
			ref_count_ptr<CheckConvergenceIP_Strategy>
				check_convergence_strategy = mmp::rcp(new CheckConvergenceIP_Strategy());

			if(options_.get()) 
				{
				CheckConvergence_StrategySetOptions
					opt_setter( check_convergence_strategy.get() );
				opt_setter.set_options( *options_ );
				}
			
			ref_count_ptr<CheckConvergenceStd_AddedStep>
				_check_convergence_step = mmp::rcp(new CheckConvergenceStd_AddedStep(check_convergence_strategy));
			
			check_convergence_step = _check_convergence_step;
			}

		// MeritFuncPenaltyParamUpdate_Step
		algo_step_ptr_t    merit_func_penalty_param_update_step = mmp::null;
		if( cov_.line_search_method_ == LINE_SEARCH_FILTER ) {
		  // We don't need to update a penalty parameter for the filter method :-)
		}
		else if( cov_.line_search_method_ != LINE_SEARCH_NONE ) {
			ref_count_ptr<MeritFunc_PenaltyParamUpdate_AddedStep>
				param_update_step = mmp::null;
			switch( cov_.merit_function_type_ ) {
				case MERIT_FUNC_L1: {
					switch(cov_.l1_penalty_param_update_) {
						case L1_PENALTY_PARAM_WITH_MULT:
//							param_update_step
//								= mmp::rcp(new  MeritFunc_PenaltyParamUpdateWithMult_AddedStep());
							THROW_EXCEPTION(
								true, std::logic_error
								,"NLPAlgoConfigIP::config_algo_cntr(...) : Error, "
								"The l1_penalty_parameter_update option of MULT_FREE has not been updated yet!" );
							break;
						case L1_PENALTY_PARAM_MULT_FREE:
							param_update_step
								= mmp::rcp(new  MeritFunc_PenaltyParamUpdateMultFree_AddedStep());
							break;
						default:
							assert(0);
					}
					break;
				}
				case MERIT_FUNC_MOD_L1:
				case MERIT_FUNC_MOD_L1_INCR:
//					param_update_step = new  MeritFunc_PenaltyParamsUpdateWithMult_AddedStep(
//											mmp::rcp_implicit_cast<MeritFuncNLP>(merit_func) );
					THROW_EXCEPTION(
						true, std::logic_error
						,"NLPAlgoConfigIP::config_algo_cntr(...) : Error, "
						"The merit_function_type options of MODIFIED_L1 and MODIFIED_L1_INCR have not been updated yet!" );
					break;
				default:
					assert(0);	// local programming error
			}
			if(options_.get()) {
				MeritFunc_PenaltyParamUpdate_AddedStepSetOptions
					ppu_options_setter( param_update_step.get() );
				ppu_options_setter.set_options( *options_ );
			}
			merit_func_penalty_param_update_step = param_update_step;
		}

		// LineSearch_Step
		algo_step_ptr_t    line_search_full_step_step = mmp::null;
		{
			line_search_full_step_step = mmp::rcp(new LineSearchFullStep_Step(bounds_tester));
		}

		// LineSearch_Step
		algo_step_ptr_t    line_search_step = mmp::null;
		if( cov_.line_search_method_ != LINE_SEARCH_NONE ) {
			ref_count_ptr<DirectLineSearchArmQuad_Strategy>
				direct_line_search = mmp::rcp(new  DirectLineSearchArmQuad_Strategy());
			if(options_.get()) {
				ConstrainedOptPack::DirectLineSearchArmQuad_StrategySetOptions
					ls_options_setter( direct_line_search.get(), "DirectLineSearchArmQuadSQPStep" );
				ls_options_setter.set_options( *options_ );
			}
			switch( cov_.line_search_method_ ) {
				case LINE_SEARCH_DIRECT: {
					line_search_step = mmp::rcp(new LineSearchDirect_Step(direct_line_search));
					break;
				}
				case LINE_SEARCH_2ND_ORDER_CORRECT: {
					THROW_EXCEPTION(
						true, std::logic_error
						,"NLPAlgoConfigIP::config_algo_cntr(...) : Error, "
						"The line_search_method option of 2ND_ORDER_CORRECT has not been updated yet!" );
					break;
				}
				case LINE_SEARCH_WATCHDOG: {
					THROW_EXCEPTION(
						true, std::logic_error
						,"NLPAlgoConfigIP::config_algo_cntr(...) : Error, "
						"The line_search_method option of WATCHDOG has not been updated yet!" );
					break;
				}
				case LINE_SEARCH_FILTER: 
					{
					mmp::ref_count_ptr<LineSearchFilter_Step> 
						line_search_filter_step = mmp::rcp(new LineSearchFilter_Step(barrier_nlp, barrier_obj_name, grad_barrier_obj_name));

					if(options_.get()) 
						{
						LineSearchFilter_StepSetOptions options_setter(line_search_filter_step.get());
						options_setter.set_options(*options_);
						}

					line_search_step = line_search_filter_step;					
					break; 
					}
			}
		}

		// LineSearchFailure
		if( new_decomp_selection_strategy.get() ) {
			line_search_step = mmp::rcp(
				new LineSearchFailureNewDecompositionSelection_Step(
					line_search_step
					,new_decomp_selection_strategy
					)
				);
		}
	
		//
		// Create the algorithm depending on the type of NLP we are trying to solve.
		//

		if( m == 0 && mI == 0 ) {
			if( nb == 0 ) {
				//
				// Unconstrained NLP (m == 0, mI == 0, num_bounded_x == 0)
				//
				if(trase_out)
					*trase_out 
						<< "\nConfiguring an algorithm for an unconstrained "
						<< "NLP (m == 0, mI == 0, num_bounded_x == 0) ...\n";
				THROW_EXCEPTION(
					m == 0 && mI == 0 && nb == 0, std::logic_error
					,"NLPAlgoConfigIP::config_alg_cntr(...) : Error, "
					"Unconstrained NLPs are not supported yet!" );
			}
			else {
				//
				// Simple bound constrained NLP (m == 0, mI == 0, num_bounded_x > 0)
				//
				if(trase_out)
					*trase_out 
						<< "\nConfiguring an algorithm for a simple bound constrained "
						<< "NLP (m == 0, mI == 0, num_bounded_x > 0) ...\n";
				THROW_EXCEPTION(
					m == 0 && mI == 0 && nb == 0, std::logic_error
					,"NLPAlgoConfigIP::config_alg_cntr(...) : Error, "
					"Bound constrained NLPs are not supported yet!" );
			}
		}
		else if( n == m ) {
			//
			// System of Nonlinear equations (n == m)
			//
			if(trase_out)
				*trase_out 
					<< "\nConfiguring an algorithm for a system of nonlinear equations "
					<< "NLP (n == m) ...\n";
			THROW_EXCEPTION(
				n == m, std::logic_error
				,"NLPAlgoConfigIP::config_alg_cntr(...) : Error, "
				"Nonlinear equation (NLE) problems are not supported yet!" );
			assert(0); // ToDo: add the step objects for this algorithm
		}
		else if ( m > 0 || mI > 0 || nb > 0 ) {
			//
			// General nonlinear NLP ( m > 0 || mI > 0 )
			//
			if(  mI == 0 && nb == 0 ) {
				//
				// Nonlinear equality constrained NLP ( m > 0 && mI == 0 && num_bounded_x == 0 )
				//
				if(trase_out)
					*trase_out 
						<< "\nConfiguring an algorithm for a nonlinear equality constrained "
						<< "NLP ( m > 0 && mI == 0 && num_bounded_x == 0) ...\n";
			}
			else {
				//
				// Nonlinear inequality constrained NLP ( mI > 0 || num_bounded_x > 0 )
				//
				if(trase_out)
					*trase_out 
						<< "\nConfiguring an algorithm for a nonlinear generally constrained "
						<< "NLP ( mI > 0 || num_bounded_x > 0 ) ...\n";
			}



			///***********************************************************
			// Add all the steps to the algorithm
			///***********************************************************
			
			int step_num       = 0;
			int assoc_step_num = 0;
	
	 		// UpdateBarrierParameter
			//algo->insert_step( ++step_num, "UpdateBarrierParameter", updateBarrierParameter_step );

			// EvalNewPoint
			algo->insert_step( ++step_num, EvalNewPoint_name, eval_new_point_step );

			//* EvalNewPoint pre steps
			// PreEvalNewPointBarrier
			algo->insert_assoc_step( step_num, IterationPack::PRE_STEP, 1, "PreEvalNewPointBarrier", preEvalNewPointBarrier_step);

			//* EvalNewPoint post steps
			if( check_descent_range_space_step_step.get() && tailored_approach && algo->algo_cntr().check_results() )
				{
				algo->insert_assoc_step(
				  step_num
				  ,IterationPack::POST_STEP
				  ,++assoc_step_num
				  ,"CheckDescentRangeSpaceStep"
				  ,check_descent_range_space_step_step
				  );
				}

			// PostEvalNewPointBarrier
			algo->insert_assoc_step( step_num, IterationPack::POST_STEP, ++assoc_step_num, "PostEvalNewPointBarrier", postEvalNewPointBarrier_step);
			assoc_step_num = 0;
			
			// ReducedGradient
			if( !tailored_approach ) {
				algo->insert_step( ++step_num, ReducedGradient_name, reduced_gradient_step );
				} 

			//			if( mI == 0 && nb == 0 ) 
			//{

 				// CalcReducedGradLagrangian
				algo->insert_step( ++step_num, CalcReducedGradLagrangian_name, calc_reduced_grad_lagr_step );

				// CalcLagrangeMultDecomposed
				// Compute these here so that in case we converge we can report them
				if( !tailored_approach ) {
	 				// ToDo: Insert this step
				}

				// CheckConvergence
				algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );

				//} 

	 		// UpdateBarrierParameter
			algo->insert_step( ++step_num, "UpdateBarrierParameter", updateBarrierParameter_step );

			// RangeSpaceStep
			if( !tailored_approach ) {
				algo->insert_step( ++step_num, RangeSpaceStep_name, range_space_step_step );
				assoc_step_num = 0;
				if( check_decomp_from_py_step.get() )
					algo->insert_assoc_step(
						step_num
						,IterationPack::POST_STEP
						,++assoc_step_num
						,"CheckDecompositionFromPy"
						,check_decomp_from_py_step
						);
				if( check_decomp_from_Rpy_step.get() )
					algo->insert_assoc_step(
						step_num
						,IterationPack::POST_STEP
						,++assoc_step_num
						,"CheckDecompositionFromRPy"
						,check_decomp_from_Rpy_step
						);
				if( check_descent_range_space_step_step.get() )
					algo->insert_assoc_step(
						step_num
						,IterationPack::POST_STEP
						,++assoc_step_num
						,"CheckDescentRangeSpaceStep"
						,check_descent_range_space_step_step
						);
			}

			/*			// ReducedGradient
 			if( !tailored_approach ) {
	 			algo->insert_step( ++step_num, ReducedGradient_name, reduced_gradient_step );
		 	}

			if( mI == 0 && nb == 0 ) {

				// CalcReducedGradLagrangian
				algo->insert_step( ++step_num, CalcReducedGradLagrangian_name, calc_reduced_grad_lagr_step );

				// CalcLagrangeMultDecomposed
				// Compute these here so that in case we converge we can report them
				if( !tailored_approach ) {
					// ToDo: Insert this step
				}

				// CheckConvergence
				algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );
				}*/

			// ReducedHessian
			algo->insert_step( ++step_num, ReducedHessian_name, reduced_hessian_step );

			// (.-1) CheckSkipBFGSUpdate
			/*algo->insert_assoc_step(
				step_num
				,IterationPack::PRE_STEP
				,1
				,CheckSkipBFGSUpdate_name
				,check_skip_bfgs_update_step
				);*/

			// UpdateReducedSigma_Step
			algo->insert_step( ++step_num, "UpdateReducedSigma", updateReducedSigma_step);

			// NullSpaceStep
			algo->insert_step( ++step_num, "NullSpaceStepIP", null_space_step);
			/*algo->insert_step( ++step_num, NullSpaceStep_name, null_space_step_step );
			if( mI > 0 || nb > 0 ) {
				// SetDBoundsStd
				algo->insert_assoc_step(
					step_num
					,IterationPack::PRE_STEP
					,1
					,"SetDBoundsStd"
					,set_d_bounds_step
				  );
				  }*/

			// CalcDFromYPYZPZ
			algo->insert_step( ++step_num, CalcDFromYPYZPZ_name, calc_d_from_Ypy_Zpy_step );
			
			// CalcD_vStep_Step
			algo->insert_step( ++step_num, "CalcD_vStep_Step", calc_d_v_step_step );

			// PreProcessBarrierLineSearch_Step
			algo->insert_step( ++step_num, "PreProcessBarrierLineSearch_Step", preprocess_barrier_linesearch_step );

			/*			if( mI > 0 || nb > 0 ) {

				// CalcReducedGradLagrangian
				algo->insert_step( ++step_num, CalcReducedGradLagrangian_name, calc_reduced_grad_lagr_step );

				// CalcLagrangeMultDecomposed
				// Compute these here so that in case we converge we can report them
				if( !tailored_approach ) {
					// ToDo: Insert this step
				}

				// CheckConvergence
				algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );
				}*/

			// LineSearch
			if( cov_.line_search_method_ == LINE_SEARCH_NONE ) {
				algo->insert_step( ++step_num, LineSearch_name, line_search_full_step_step );
			}
			else {
				// Main line search step
				algo->insert_step( ++step_num, LineSearch_name, line_search_step );
				// Insert presteps
				Algorithm::poss_type
					pre_step_i = 0;
				// (.-?) LineSearchFullStep
				//algo->insert_assoc_step(
				//	step_num
				//	,IterationPack::PRE_STEP
				//	,++pre_step_i
				//	,"LineSearchFullStep"
				//	,line_search_full_step_step
				//	);
				// (.-?) MeritFunc_PenaltyPramUpdate
				if(merit_func_penalty_param_update_step.get()) {
				  algo->insert_assoc_step(
					  step_num
					  ,IterationPack::PRE_STEP
					  ,++pre_step_i
					  ,"MeritFunc_PenaltyParamUpdate"
					  ,merit_func_penalty_param_update_step
					  );
				}
			}

			// PostProcessBarrierLineSearch_Step
			algo->insert_step( ++step_num, "PostProcessBarrierLineSearch_Step", postprocess_barrier_linesearch_step );

			// CheckConvergence
			//algo->insert_step( ++step_num, CheckConvergence_name, check_convergence_step );

		}
		else {
			assert(0); // Error, this should not ever be called!
		}
	}
	
}

void NLPAlgoConfigIP::init_algo(NLPAlgoInterface* _algo)
{
	using DynamicCastHelperPack::dyn_cast;
	namespace mmp = MemMngPack;

	THROW_EXCEPTION(
		_algo == NULL, std::invalid_argument
		,"NLPAlgoConfigIP::init_algo(_algo) : Error, "
		"_algo can not be NULL" );

	NLPAlgo             &algo    = dyn_cast<NLPAlgo>(*_algo);
	NLPAlgoState	         &state   = algo.rsqp_state();
	NLP			         &nlp     = algo.nlp();
	NLPVarReductPerm     *nlp_vrp = dynamic_cast<NLPVarReductPerm*>(&nlp);
	NLPDirect  *nlp_fod = dynamic_cast<NLPDirect*>(&nlp);

	algo.max_iter( algo.algo_cntr().max_iter() );
	algo.max_run_time( algo.algo_cntr().max_run_time() );

	// Reset the iteration count to zero
	state.k(0);

	// Get organized output of vectors and matrices even if setw is not used by Step objects.
	algo.track().journal_out()
		<< std::setprecision(algo.algo_cntr().journal_print_digits())
		<< std::scientific;

	// set the first step
	algo.do_step_first(1);

	// The rest of the algorithm should initialize itself
}

// private

void NLPAlgoConfigIP::readin_options(
	  const OptionsFromStreamPack::OptionsFromStream     &options
	, SOptionValues                                      *ov
	, std::ostream                                       *trase_out
	)
{
	namespace	ofsp = OptionsFromStreamPack;
	using		ofsp::OptionsFromStream;
	typedef		OptionsFromStream::options_group_t		options_group_t;
	using		ofsp::StringToIntMap;
	using		ofsp::StringToBool;

	assert(ov);	// only a local class error

	// Get the options group for "NLPAlgoConfigIP"
	const std::string opt_grp_name = "NLPAlgoConfigIP";
	const OptionsFromStream::options_group_t optgrp = options.options_group( opt_grp_name );
	if( OptionsFromStream::options_group_exists( optgrp ) ) {

		// Define map for options group "IpConfig".
		const int num_opts = 11;
		enum EIpConfig {
			MAX_BASIS_COND_CHANGE_FRAC
			,EXACT_REDUCED_HESSIAN
			,QUASI_NEWTON
			,NUM_LBFGS_UPDATES_STORED
			,LBFGS_AUTO_SCALING
			,HESSIAN_INITIALIZATION
			,QP_SOLVER
			,REINIT_HESSIAN_ON_QP_FAIL
			,LINE_SEARCH_METHOD
			,MERIT_FUNCTION_TYPE
			,L1_PENALTY_PARAM_UPDATE
		};
		const char* SIpConfig[num_opts]	= {
			"max_basis_cond_change_frac"
			,"exact_reduced_hessian"
			,"quasi_newton"
			,"num_lbfgs_updates_stored"
			,"lbfgs_auto_scaling"
			,"hessian_initialization"
			,"qp_solver"
			,"reinit_hessian_on_qp_fail"
			,"line_search_method"
			,"merit_function_type"
			,"l1_penalty_parameter_update"
		};
		StringToIntMap	map(	opt_grp_name, num_opts, SIpConfig );

		options_group_t::const_iterator itr = optgrp.begin();
		for( ; itr != optgrp.end(); ++itr ) {
			switch( (EIpConfig)map( ofsp::option_name(itr) ) ) {
				case MAX_BASIS_COND_CHANGE_FRAC:
					ov->max_basis_cond_change_frac_ = ::atof( ofsp::option_value(itr).c_str() );
					break;
				case EXACT_REDUCED_HESSIAN:
					ov->exact_reduced_hessian_ = StringToBool( "exact_reduced_hessian", ofsp::option_value(itr).c_str() );
					break;
				case QUASI_NEWTON:
				{
					const std::string &opt_val = ofsp::option_value(itr);
					if( opt_val == "AUTO" )
						ov->quasi_newton_ = QN_AUTO;
					else if( opt_val == "BFGS" )
						ov->quasi_newton_ = QN_BFGS;
					else if( opt_val == "PBFGS" )
						ov->quasi_newton_ = QN_PBFGS;
					else if( opt_val == "LBFGS" )
						ov->quasi_newton_ = QN_LBFGS;
					else if( opt_val == "LPBFGS" )
						ov->quasi_newton_ = QN_LPBFGS;
					else
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"NLPAlgoConfigIP::readin_options(...) : "
							"Error, incorrect value for \"quasi_newton\" "
							", Only options of BFGS, PBFGS"
							", LBFGS, LPBFGS and AUTO are avalible."
							);
					break;
				}
				case NUM_LBFGS_UPDATES_STORED:
					ov->num_lbfgs_updates_stored_ = ::atoi( ofsp::option_value(itr).c_str() );
					break;
				case LBFGS_AUTO_SCALING:
					ov->lbfgs_auto_scaling_
						= StringToBool( "lbfgs_auto_scaling", ofsp::option_value(itr).c_str() );
					break;
				case HESSIAN_INITIALIZATION:
				{
					const std::string &opt_val = ofsp::option_value(itr);
					if( opt_val == "IDENTITY" )
						ov->hessian_initialization_ = INIT_HESS_IDENTITY;
					else if( opt_val == "FINITE_DIFF_SCALE_IDENTITY" )
						ov->hessian_initialization_ = INIT_HESS_FIN_DIFF_SCALE_IDENTITY;
					else if( opt_val == "FINITE_DIFF_DIAGONAL" )
						ov->hessian_initialization_ = INIT_HESS_FIN_DIFF_SCALE_DIAGONAL;
					else if( opt_val == "FINITE_DIFF_DIAGONAL_ABS" )
						ov->hessian_initialization_ = INIT_HESS_FIN_DIFF_SCALE_DIAGONAL_ABS;
					else if( opt_val == "AUTO" )
						ov->hessian_initialization_ = INIT_HESS_AUTO;
					else
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"NLPAlgoConfigIP::readin_options(...) : "
							"Error, incorrect value for \"hessian_initialization\" "
							", Only options of IDENTITY, FINITE_DIFF_SCALE_IDENTITY,"
							" FINITE_DIFF_DIAGONAL, FINITE_DIFF_DIAGONAL_ABS and AUTO"
							" are available"  );
					break;
				}
				case QP_SOLVER:
				{
					const std::string &qp_solver = ofsp::option_value(itr);
					if( qp_solver == "AUTO" ) {
						ov->qp_solver_type_ = QP_AUTO;
					} else if( qp_solver == "QPSOL" ) {
						ov->qp_solver_type_ = QP_QPSOL;
					} else if( qp_solver == "QPOPT" ) {
#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT
						ov->qp_solver_type_ = QP_QPOPT;
#else
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"NLPAlgoConfigIP::readin_options(...) : QPOPT is not supported,"
							" must define CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT!" );
#endif
					} else if( qp_solver == "QPKWIK" ) {
						ov->qp_solver_type_ = QP_QPKWIK;
					} else if( qp_solver == "QPSCHUR" ) {
						ov->qp_solver_type_ = QP_QPSCHUR;
					} else {
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"NLPAlgoConfigIP::readin_options(...) : "
							"Error, incorrect value for \"qp_solver\" "
							"Only qp solvers QPOPT, QPSOL, QPKWIK, QPSCHUR and AUTO are avalible."	);
					}
					break;
				}
				case REINIT_HESSIAN_ON_QP_FAIL:
					ov->reinit_hessian_on_qp_fail_ = StringToBool( "reinit_hessian_on_qp_fail", ofsp::option_value(itr).c_str() );
					break;
				case LINE_SEARCH_METHOD:
				{
					const std::string &option = ofsp::option_value(itr);
					if( option == "NONE" ) {
						ov->line_search_method_ = LINE_SEARCH_NONE;
					} else if( option == "DIRECT" ) {
						ov->line_search_method_ = LINE_SEARCH_DIRECT;
					} else if( option == "2ND_ORDER_CORRECT" ) {
						ov->line_search_method_ = LINE_SEARCH_2ND_ORDER_CORRECT;
					} else if( option == "WATCHDOG" ) {
						ov->line_search_method_ = LINE_SEARCH_WATCHDOG;
					} else if( option == "AUTO" ) {
						ov->line_search_method_ = LINE_SEARCH_AUTO;
					} else if( option == "FILTER" ) {
						ov->line_search_method_ = LINE_SEARCH_FILTER;
					} else {
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"NLPAlgoConfigIP::readin_options(...) : "
							"Error, incorrect value for \"line_search_method\".\n"
							"Only the options NONE, DIRECT, 2ND_ORDER_CORRECT, FILTER, WATCHDOG "
							"and AUTO are avalible." );
					}
					break;
				}
				case MERIT_FUNCTION_TYPE:
				{
					const std::string &option = ofsp::option_value(itr);
					if( option == "L1" )
						ov->merit_function_type_ = MERIT_FUNC_L1;
					else if( option == "MODIFIED_L1" )
						ov->merit_function_type_ = MERIT_FUNC_MOD_L1;
					else if( option == "MODIFIED_L1_INCR" )
						ov->merit_function_type_ = MERIT_FUNC_MOD_L1_INCR;
					else if( option == "AUTO" )
						ov->merit_function_type_ = MERIT_FUNC_AUTO;
					else
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"NLPAlgoConfigIP::readin_options(...) : "
							"Error, incorrect value for \"merit_function_type\".\n"
							"Only the options L1, MODIFIED_L1, MODIFIED_L1_INCR "
							"and AUTO are avalible." );
					break;
				}
				case L1_PENALTY_PARAM_UPDATE:
				{
					const std::string &option = ofsp::option_value(itr);
					if( option == "WITH_MULT" )
						ov->l1_penalty_param_update_
							= L1_PENALTY_PARAM_WITH_MULT;
					else if( option == "MULT_FREE" )
						ov->l1_penalty_param_update_
							= L1_PENALTY_PARAM_MULT_FREE;
					else if( option == "AUTO" )
						ov->l1_penalty_param_update_
							= L1_PENALTY_PARAM_AUTO;
					else
						THROW_EXCEPTION(
							true, std::invalid_argument
							,"NLPAlgoConfigIP::readin_options(...) : "
							"Error, incorrect value for \"l1_penalty_param_update\".\n"
							"Only the options WITH_MULT, MULT_FREE and AUTO"
							"are avalible."  );
					break;
				}
				default:
					assert(0);	// this would be a local programming error only.
			}
		}
	}
	else {
		if(trase_out)
			*trase_out
				<< "\n\n*** Warning!  The options group \"NLPAlgoConfigIP\" was not found.\n"
				<< "Using a default set of options instead ... \n";
	}
}

//
// This is where some of the default options are set and the user is alerted to what their
// value is.
//
void NLPAlgoConfigIP::set_default_options( 
	const SOptionValues     &uov
	,SOptionValues          *cov
	,std::ostream           *trase_out
	)
{
	if(trase_out)
		*trase_out
			<< "\n*** Setting option defaults for options not set by the user or determined some other way ...\n";

	if( cov->max_basis_cond_change_frac_ < 0.0 &&  uov.max_basis_cond_change_frac_ < 0.0 ) {
		if(trase_out)
			*trase_out
				<< "\nmax_basis_cond_change_frac < 0 : setting max_basis_cond_change_frac = 1e+4 \n";
		cov->max_basis_cond_change_frac_ = 1e+4;
	}
	else {
		cov->max_basis_cond_change_frac_ = uov.max_basis_cond_change_frac_;
	}
	cov->exact_reduced_hessian_ = uov.exact_reduced_hessian_;
	if( cov->quasi_newton_ == QN_AUTO && uov.quasi_newton_ == QN_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nquasi_newton == AUTO: setting quasi_newton = BFGS\n";
		cov->quasi_newton_ = QN_BFGS;
	}
	else if(cov->quasi_newton_ == QN_AUTO) {
		cov->quasi_newton_ = uov.quasi_newton_;
	}
	if( cov->num_lbfgs_updates_stored_ < 0 && uov.num_lbfgs_updates_stored_ < 0 ) {
		if(trase_out)
			*trase_out
				<< "\nnum_lbfgs_updates_stored < 0 : setting num_lbfgs_updates_stored = 10\n";
		cov->num_lbfgs_updates_stored_ = 10;
	}
	else if(cov->num_lbfgs_updates_stored_ < 0) {
		cov->num_lbfgs_updates_stored_ = uov.num_lbfgs_updates_stored_;
	}
	cov->lbfgs_auto_scaling_ = uov.lbfgs_auto_scaling_;
	if( cov->hessian_initialization_ == INIT_HESS_AUTO && uov.hessian_initialization_ == INIT_HESS_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nhessian_initialization == AUTO: setting hessian_initialization = FINITE_DIFF_DIAGONAL_ABS\n";
		cov->hessian_initialization_ = INIT_HESS_FIN_DIFF_SCALE_DIAGONAL_ABS;
	}
	else if(cov->hessian_initialization_ == INIT_HESS_AUTO) {
		cov->hessian_initialization_ = uov.hessian_initialization_;
	}
	if( cov->qp_solver_type_ == QP_AUTO && uov.qp_solver_type_ == QP_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nqp_solver_type == AUTO: setting qp_solver_type = QPSCHUR\n";
		cov->qp_solver_type_ = QP_QPSCHUR;
	}
	else if(cov->qp_solver_type_ == QP_AUTO) {
		cov->qp_solver_type_ = uov.qp_solver_type_;
	}
	cov->reinit_hessian_on_qp_fail_ = uov.reinit_hessian_on_qp_fail_;
	if( cov->line_search_method_ == LINE_SEARCH_AUTO && uov.line_search_method_ == LINE_SEARCH_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nline_search_method == AUTO: setting line_search_method = FILTER\n";
		cov->line_search_method_ = LINE_SEARCH_FILTER;
	}
	else if(cov->line_search_method_ == LINE_SEARCH_AUTO) {
		cov->line_search_method_ = uov.line_search_method_;
	}
	if( cov->merit_function_type_ == MERIT_FUNC_AUTO && uov.merit_function_type_ == MERIT_FUNC_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nmerit_function_type == AUTO: setting merit_function_type = MODIFIED_L1_INCR\n";
		cov->merit_function_type_ = MERIT_FUNC_MOD_L1_INCR;
	}
	else if(cov->merit_function_type_ == MERIT_FUNC_AUTO) {
		cov->merit_function_type_ = uov.merit_function_type_;
	}
	if( cov->l1_penalty_param_update_ == L1_PENALTY_PARAM_AUTO && uov.l1_penalty_param_update_ == L1_PENALTY_PARAM_AUTO ) {
		if(trase_out)
			*trase_out
				<< "\nl1_penalty_param_update == AUTO: setting l1_penalty_param_update = MULT_FREE\n";
		cov->l1_penalty_param_update_ = L1_PENALTY_PARAM_MULT_FREE;
	}
	else if(cov->l1_penalty_param_update_ == L1_PENALTY_PARAM_AUTO) {
		cov->l1_penalty_param_update_ = uov.l1_penalty_param_update_;
	}
	if( cov->full_steps_after_k_ < 0 && uov.full_steps_after_k_ < 0 ) {
		if(trase_out)
			*trase_out
				<< "\nfull_steps_after_k < 0 : the line search will never be turned off after so many iterations\n";
	}
	else {
		cov->full_steps_after_k_ = uov.full_steps_after_k_;
	}
	if(trase_out)
		*trase_out
			<< "\n*** End setting default options\n";
}

}	// end namespace MoochoPack
