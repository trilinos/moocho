// ////////////////////////////////////////////////////////////////
// ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions.h

#ifndef REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_SET_OPTIONS_H
#define REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_SET_OPTIONS_H

#include "ReducedHessianSecantUpdateBFGSProjected_Strategy.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Set options for \Ref{ReducedHessianSecantUpdateBFGSProjected_Strategy} from a
  * OptionsFromStream object.
  *
  * The options group is (with the default name):
  *
  \begin{verbatim}
	options_group ReducedHessianSecantUpdateBFGSProjected {
		act_set_frac_proj_start   = 0.8;   *** (+dbl)
		project_error_tol         = 1e-5;  *** (+dbl) [0.0, 1.0]
        super_basic_mult_drop_tol = 1e-5;  *** (+dbl)
	}
  \end{verbatim}
  *
  * \begin{description}
  *	\item[act_set_frac_proj_start] ToDo : Finish.
  *	\item[project_error_tol] ToDo : Finish.
  * \item[super_basic_mult_drop_tol] ToDo: Finish.
  *	\end{description}
  */
class ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			ReducedHessianSecantUpdateBFGSProjected_Strategy >
{
public:

	///
	ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions(
		ReducedHessianSecantUpdateBFGSProjected_Strategy* target = 0
		, const char opt_grp_name[] = "ReducedHessianSecantUpdatePBFGS" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class ReducedHessianSecantUpdateBFGSProjected_StrategySetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// REDUCED_HESSIAN_SECANT_UPDATE_BFGS_PROJECTED_STRATEGY_SET_OPTIONS_H