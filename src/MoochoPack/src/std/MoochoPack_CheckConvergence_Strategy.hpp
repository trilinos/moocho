// ////////////////////////////////////////////////////////////////////////////
// CheckConvergence_Strategy.hpp
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

#ifndef CHECK_CONVERGENCE_STRATEGY_H
#define CHECK_CONVERGENCE_STRATEGY_H

#include "MoochoPack/src/MoochoPackTypes.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "SetOptionsFromStreamNode.hpp"
#include "SetOptionsToTargetBase.hpp"

namespace MoochoPack {

///
/** Strategy interface for performing convergence checks.
 *
 * This object can not change the flow of control or do anything fancy.  It just
 *  checks convergence by calculating norm errors and comparing with tolerance
 *  It can update iteration quantities if desired.
 *
 * See the printed documentation generated by \c this->print_step().
 */
class CheckConvergence_Strategy 
	{
	public:
		
		enum EOptErrorCheck 
			{ 
			OPT_ERROR_REDUCED_GRADIENT_LAGR, 
			OPT_ERROR_GRADIENT_LAGR 
			};

		///
		/** <<std member comp>> members for whether to check the reduced
		 * or full gradient of the Lagrangian. 
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( EOptErrorCheck, opt_error_check )

		enum EScaleKKTErrorBy 
			{ 
			SCALE_BY_ONE, 
			SCALE_BY_NORM_2_X, 
			SCALE_BY_NORM_INF_X 
			};

		///
		/** <<std member comp>> members for how the optimality condition should
		 *   be scaled
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( EScaleKKTErrorBy, scale_opt_error_by )
		

		///
		/** <<std member comp>> members for how the feasibility condition should
		 *   be scaled
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( EScaleKKTErrorBy, scale_feas_error_by )


		///
		/** <<std member comp>> members for how the complimentarity condition should
		 *   be scaled
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( EScaleKKTErrorBy, scale_comp_error_by )


		///
		/** <<std member comp>> members for whether the optimality conditions
		 * should be scaled by the 
		 */
		STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, scale_opt_error_by_Gf )

		///
		CheckConvergence_Strategy(
		  EOptErrorCheck opt_error_check = OPT_ERROR_REDUCED_GRADIENT_LAGR,
		  EScaleKKTErrorBy scale_opt_error_by = SCALE_BY_ONE,
		  EScaleKKTErrorBy scale_feas_error_by = SCALE_BY_ONE,
		  EScaleKKTErrorBy scale_comp_error_by = SCALE_BY_ONE,
		  bool scale_opt_error_by_Gf = true
		  );
		
		///
		virtual bool Converged( Algorithm& _algo)=0;

		///
		virtual void print_step( const Algorithm& _algo, std::ostream& out, const std::string& L ) const =0;


	}; // end interface CheckConvergence_Strategy


///
/** Set options for CheckConvergence_Strategy from an
  * OptionsFromStream object.
  *
  * The default options group name is CheckConvergenceStrategy.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group CheckConvergenceStrategy {
		scale_kkt_error_by   = SCALE_BY_ONE;
		scale_opt_error_by_Gf = true;
	}
  \end{verbatim}
  *
  * \begin{description}
  *	\item[scale_kkt_error_by] Determines if and how the optimality (opt_kkt_err)
  *		and feasiblity (feas_kkt_err)
  *		errors for the convergence check are scaled by for the unkowns x before
  *		comparing it to the set tolerances of opt_tol and feas_tol (see the
  *		class \Ref{CheckConvergenceStd_AddedStep} and its printed algorithm
  *		for more details).
  *		\begin{description}
  *		\item[SCALE_BY_ONE]			no scaling by x
  *		\item[SCALE_BY_NORM_2_X]    scale opt_kkt_err and feas_kkt_err by 1/||x||2
  *		\item[SCALE_BY_NORM_INF_X]  scale opt_kkt_err and feas_kkt_err by 1/||x||inf
  *		\end{description}
  *	\item[scale_opt_error_by_Gf] Determines if opt_kkt_err is scaled by
  *		||Gf_k||inf or not.
  *	\end{description}
  */
class CheckConvergence_StrategySetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode, 
	  public OptionsFromStreamPack::SetOptionsToTargetBase<
	          CheckConvergence_Strategy >
	{
	public:
		
		///
		CheckConvergence_StrategySetOptions(
		  CheckConvergence_Strategy* target = 0
		  , const char opt_grp_name[] = "CheckConvergenceStrategy" );
		
	protected:
		
		/// Overridden from SetOptionsFromStreamNode
		void setOption( int option_num, const std::string& option_value );

	};	// end class CheckConvergence_Strategy


} // end namespace MoochoPack

#endif // CHECK_CONVERGENCE_STRATEGY_H
