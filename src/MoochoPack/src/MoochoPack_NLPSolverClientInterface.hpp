// ////////////////////////////////////////////////////////////////////////////
// rSQPSolverClientInterface.hpp
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

#ifndef RSQP_SOLVER_CLIENT_INTERFACE_H
#define RSQP_SOLVER_CLIENT_INTERFACE_H

#include <stdexcept>

#include "ReducedSpaceSQPPackTypes.hpp"
#include "NLPInterfacePack/src/NLP.hpp"
#include "StandardCompositionMacros.hpp"
#include "StandardMemberCompositionMacros.hpp"

namespace ReducedSpaceSQPPack {

///
/** This is the most basic interface that clients use to solve an NLP.
 *
 * ToDo: Finish documentaiton.
 */
class rSQPSolverClientInterface {
public:

	/** @name Public Types */
	//@{

	///
	enum EFindMinReturn { SOLUTION_FOUND, MAX_ITER_EXCEEDED, MAX_RUN_TIME_EXCEEDED
		, ALGORITHMIC_ERROR };

	/// Thrown if the setup is not valid
	class InvalidSetup : public std::logic_error
	{public: InvalidSetup(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}
	
	/** @name Constructors/initalizers */
	//@{

	/// <<std comp>> members for the nlp
	STANDARD_COMPOSITION_MEMBERS( NLP, nlp )

	/// <<std comp>> members for the track
	STANDARD_COMPOSITION_MEMBERS( AlgorithmTrack, track )

	///
	/** Construct with no references set to nlp or track objects.
	 */
	rSQPSolverClientInterface(
		int                    max_iter             = 10000
		,double                max_run_time         = 1e+10 // run forever
		,value_type            opt_tol              = 1e-6
		,value_type            feas_tol             = 1e-6
		,value_type            comp_tol             = 1e-6
		,value_type            step_tol             = 1e-2
		,EJournalOutputLevel   journal_output_level = PRINT_NOTHING
		,int                   journal_print_digits = 6
		,bool                  check_results        = false
		,bool                  calc_conditioning    = false
		);

	///
	virtual ~rSQPSolverClientInterface() {}

	//@}

	/** @name Solver Parameters */
	//@{

	/// Set the maximum number of iterations the rSQP algorithm can perform
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, max_iter )

	///
	/** Set the maximum run_time
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( double, max_run_time )

	///
	/** Set the termination tolerance for the relative (scaled) linear dependence of the
	 * gradients part of the first order necessary optimality conditions.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, opt_tol )

	///
	/** Set the termination tolerance for the (scaled) equality constraints ||c(x*)||inf
	 * which is part of the first order necessary optimality conditions.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, feas_tol )

	///
	/** Set the termination tolerance for the complementarity condition 
	 *  for the (scaled) bound constraints
	 *  which is part of the first order necessary optimality conditions.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, comp_tol )

	///
	/** Set the termination tolerance for the change in the estimate of the solution.
	 *
	 * The test is: <tt>|d(i)|/(1+|x(i)|) < step_tol</tt>. 
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, step_tol )

	///
	/** Determine the amount of output to a journal file.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EJournalOutputLevel, journal_output_level )

	///
	/** Set the precesion of the journal output.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, journal_print_digits )

	///
	/** Set whether computations will be double checked or not.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_results )

	///
	/** Set whether the condition numbers of important matrics is
	 * computed and printed or not.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, calc_conditioning )

	//@}

	/** @name Solve the NLP*/
	//@{

	///
	/** Find the minimun of the set NLP.
	 *
	 * This function returns <tt>SOLUTION_FOUND</tt> if the NLP has been solved to the desired
	 * tolerances.  In this case <tt>this->track().output_final(...,TERMINATE_TRUE)</tt>
	 * and <tt>this->nlp().report_final_solution(...,true)</tt> are called before this function
	 * returns..
	 *
	 * If the solution is not found, then <tt>this->nlp().report_final_solution(...,false)</tt> is
	 * called and one of the following occurs:
	 * <ul>
	 * <li>  If the maximum number of iterations has been exceeded then <tt>MAX_ITER_EXCEEDED</tt> will
	 * be returned.  In this case <tt>this->track().output_final(...,MAX_ITER_EXCEEDED)</tt> is called.
	 * <li> If the maximum runtime has been exceeded then <tt>MAX_RUN_TIME_EXCEEDED</tt> will be
	 * returned.  In this case <tt>this->track().output_final(...,MAX_RUN_TIME_EXCEEDED)</tt> is called.
	 * <li> An exception is thrown.  The client should be prepaired to catch any exceptions thrown
	 * from this function.
	 * All of the purposefully thrown exceptions are derived from std::exception so the
	 * client can check the what() function to see and description of the error.  If an
	 * exception is thrown then <tt>this->track().output_final(...,TERMINATE_FALSE)</tt> will
	 * be called before this exception is rethrown out of this functiion.  If the constraints
	 * are found to be infeasible, then the exception <tt>InfeasibleConstraints</tt> will be thrown.
	 * If a line search failure occurs then the exception <tt>LineSearchFailure</tt> will be thrown.
	 * If some test failed then the exception <tt>TestFailed</tt> will be thrown.  Many other exceptions
	 * may be thrown but these are the main ones that the SQP algorithms known about and will
	 * purposefully generate.
	 * </ul>
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->nlp() != 0</tt> (throw InvalidSetup)
	 * </ul>
	 *
	 * Postcondtions:<ul>
	 * <li> Minimum of %NLP is found to opt_tol, max_iter was reached
	 * or max_run_time reached (throw std::exection)
	 * </ul>
	 */
	virtual EFindMinReturn find_min() = 0;

	//@}

	/** @name Algorithm description */
	//@{

	///
	/** Prints a description of the algorithm.
	  */
	virtual void print_algorithm(std::ostream& out) const = 0;

	//@}

	/** @name Algorithm timing */
	//@{

	///
	/** Causes algorithm to be timed.
	 *
	 * Call with <tt>algo_timing == true</tt> before calling <tt>find_min()</tt>
	 * to have the algorithm timed.
	 */
	virtual void set_algo_timing( bool algo_timing ) = 0;

	///
	virtual bool algo_timing() const = 0;

	///
	/** Outputs table of times for each step and the cummulative times.
	 *
	 * Call after <tt>find_min()</tt> has executed to get a table
	 * of times.
	 */
	virtual void print_algorithm_times( std::ostream& out ) const = 0;

	//@}

private:

#ifdef DOXYGEN_COMPILE // Strictly for doxygen diagrams
	///
	NLPInterfacePack::NLP                  *nlp;
	///
	GeneralIterationPack::AlgorithmTrack   *track;
#endif

}; // end class rSQPSolverClientInterface

}  // end namespace ReducedSpaceSQPPack

#endif	// RSQP_SOLVER_CLIENT_INTERFACE_H
