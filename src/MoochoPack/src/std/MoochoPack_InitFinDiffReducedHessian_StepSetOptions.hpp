// ////////////////////////////////////////////////////////////////
// InitFinDiffReducedHessian_StepSetOptions.h

#ifndef INIT_FIN_DIFF_REDUCED_HESSIAN_STEP_SET_OPTIONS_H
#define INIT_FIN_DIFF_REDUCED_HESSIAN_STEP_SET_OPTIONS_H

#include "InitFinDiffReducedHessian_Step.h"
#include "Misc/include/SetOptionsFromStreamNode.h"
#include "Misc/include/SetOptionsToTargetBase.h"

namespace ReducedSpaceSQPPack {

///
/** Set options for InitFinDiffReducedHessian_Step from an
  * OptionsFromStream object.
  *
  * The default options group name is InitFinDiffReducedHessian.
  *
  * The options group is:
  *
  \begin{verbatim}
	options_group InitFinDiffReducedHessian {
		initialization_method	= ?;
		max_cond				= ?;
		min_diag				= ?;
		step_scale				= ?;
	}
  \end{verbatim}
  *
  * \begin{description}
  *	\item[initialization_method] Determines how the diagonal is initialized.
  *		from the finite difference taken.
  *		\begin{description}
  *		\item[SCALE_IDENTITY]      diag(i) = max( ||rGf_fd||inf , smallest_ele )
  *		\item[SCALE_DIAGONAL]      diag(i) = max( rGf_fd(i)     , smallest_ele )
  *		\item[SCALE_DIAGONAL_ABS]  diag(i) = max( abs(rGf_fd(i)), smallest_ele )
  *		\end{description}
  *		where: smallest_ele = max( ||rGf_fd||inf / max_cond , min_diag )
  *	\item[max_cond] The maximum condition of the initialized matrix.
  *		See initialization_method.
  *		Example: max_cond = 1e+5.
  *	\item[min_diag] The smallest absolute diagonal element.
  *		Example: min_diag = 1e-8.
  *	\item[step_scale] scales the step for the finite difference by
  *		u = scale_step / ||Z*e||inf.
  *		The finite difference is then taken as:
  *		rGf_fd = ( Z_k * g(x_k + u * Z*e - rGf_k ) / u
  *		Example: step_scale = 1.0.
  */
class InitFinDiffReducedHessian_StepSetOptions
	: public OptionsFromStreamPack::SetOptionsFromStreamNode 
		, public OptionsFromStreamPack::SetOptionsToTargetBase<
			InitFinDiffReducedHessian_Step >
{
public:

	///
	InitFinDiffReducedHessian_StepSetOptions(
		  InitFinDiffReducedHessian_Step* target = 0
		, const char opt_grp_name[] = "InitFinDiffReducedHessian" );

protected:

	/// Overridden from SetOptionsFromStreamNode
	void set_option( int option_num, const std::string& option_value );

};	// end class InitFinDiffReducedHessian_StepSetOptions

}	// end namespace ReducedSpaceSQPPack

#endif	// INIT_FIN_DIFF_REDUCED_HESSIAN_STEP_SET_OPTIONS_H