// ////////////////////////////////////////////////////////////////////////////////
// FeasibilityStep_Strategy.h

#ifndef FEASIBILITY_STEP_STRATEGY_H
#define FEASIBILITY_STEP_STRATEGY_H

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackTypes.h"

namespace ReducedSpaceSQPPack {

///
/** Abstract interface for a strategy object that will compute a step that will
 * improve feasibility (at least descent) {abstract}.
 */
class FeasibilityStep_Strategy {
public:

	///
	virtual ~FeasibilityStep_Strategy() {}

	///
	/** Compute a step that improves feasibility (at least locally).
	 *
	 * This function will compute a step #w# that satisfies:
	 *
	 * #d_bounds_k.l <= xo + w - x_k <= d_bounds_k.u#
	 *
	 * and gives descent for ||c(xo + beta*w)|| for at least small beta > 0.
	 * This norm ||.|| could be any valid norm and the implementation is free to
	 * define what descent means any way it would like.  Any information being used
	 * in the algorithm can be used to compute this step.
	 *
	 * @param out           [out] Output stream journal data is written to.
	 * @param olevel        [in] Output level for printing to #out#
	 * @param algo          [in/out] The rSQPAlgo object.  This object can be queryed for
	 *                      information.
	 * @param s             [in/out] rSQPState object.  May be queried or modified if needed.
	 * @param xo            [in] Base point vector (size n) xo.
	 * @param c_xo          [in] c(xo).
	 * @param w             [out] Computed step vector (size n) w.  Must not be NULL.
	 *
	 * @return Returns true if a step that reduces feasibility subject to the bounds could
	 * be found and false otherwise.
	 */
 	virtual bool compute_feasibility_step(
		std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
		,const VectorSlice& xo, const VectorSlice& c_xo, VectorSlice* w
	  	) = 0;

	///
	/** This function will print a description of the computations and logic used.
	 */
	virtual void print_step( std::ostream& out, const std::string& leading_str ) const = 0;

}; // end class FeasibilityStep_Strategy

} // end namespace ReducedSpaceSQPPack

#endif // FEASIBILITY_STEP_STRATEGY_H
