// ////////////////////////////////////////////////////////////////
// CheckConvergenceStd_AddedStepSetOptions.h

#ifndef CHECK_CONVERGENCE_STD_ADDED_STEP_SET_OPTIONS_H
#define CHECK_CONVERGENCE_STD_ADDED_STEP_SET_OPTIONS_H

#include "CheckConvergenceStd_AddedStep.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Set options for CheckConvergenceStd_AddedStep from an
  * OptionsFromStream object.
  *
  * The default options group name is CheckConvergenceStd.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group CheckConvergenceStd {
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
class CheckConvergenceStd_AddedStepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			CheckConvergenceStd_AddedStep >
{
public:

	///
	CheckConvergenceStd_AddedStepSetOptions(
		  CheckConvergenceStd_AddedStep* target = 0
		, const char opt_grp_name[] = "CheckConvergenceStd" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class CheckConvergenceStd_AddedStepSetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// CHECK_CONVERGENCE_STD_ADDED_STEP_SET_OPTIONS_H