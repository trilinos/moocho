// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxedLOQO.h

#ifndef QP_SOLVER_RELAXED_LOQO_H
#define QP_SOLVER_RELAXED_LOQO_H

#ifdef USE_QP_LOQO

#include <vector>

#include "QPSolverRelaxed.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "Misc/include/StandardMemberCompositionMacros.h"
#include "Misc/include/StandardCompositionMacros.h"

namespace ConstrainedOptimizationPack {

///
/** Solves Quadratic Programming (QP) problem using the primal-dual interior-point
  * solver LOQO (Vanderbei).
  *
  * In this implementation it is required that G support the \Ref{MatrixExtractInvCholFactor}
  * interface and is therefore quite restrictive on the likes of QPs it can solve.
  */
class QPSolverRelaxedLOQO : public QPSolverRelaxed
{
public:

	///
	/** Strategy interface for initializing the Hessian matrix Q and the
	 * Jacobian matrix A for LOQO.
	 *
	 * This interface allows clients to specialize how G, op(E), and op(F)
	 * are translated into the sparse compresssed column data structurs for
	 * Q and A for LOQO.
	 */
	class InitLOQOHessianJacobian {
	public:

		///
		virtual ~InitLOQOHessianJacobian() {}

		///
		/** Extract the Hessian Q and Jacobian A matrices.
		 *
		 * The matrices formed are:
		 *
		 \begin{verbatim}
		 Q = [ G,   0  ]     A = [ P*op(E), -P*b ]
		     [ 0, bigM ]         [  op(F),   -f  ]
		 \end{verbatim}
		 *
		 * The default implementation just assumes that G, E and F are all
		 * dense matrices.
		 *
		 * @param  G       [in] Original Hessian matrix (size nd x nd)
		 * @param  bigM    [in] Value of "big M" to use in the relaxed Hessian
		 * @param  E       [in] Pointer to inequality constriants matrix (may be NULL)
		 * @param  trans_E [in] Determines op(E) = E or op(E) = E'
		 * @param  b       [in] Relaxation vector (size m_in) for inequalities
		 * @param  loqo_b_stat
		 *                 [in] Array (size num_inequal) that determines what P is above as
		 *                 follows:
		 *                             / [ +op(E)(j,:), -b(j) ] : if j = loqo_b_stat(k) > 0
		 *                   A(k,:) =  |
		 *                             \ [ -op(E)(j,:), +b(j) ] : if j = loqo_b_stat(k) < 0
		 *                          , for k = 1...num_inequal
		 * @param  num_inequal
		 *                 [in] Number of finitly bounded inequalites in op(E).
		 *                 num_inequal <= m_in.
		 * @param  F       [in] Pointer to equality constriants matrix (may be NULL)
		 * @param  trans_E [in] Determines op(F) = F or op(F) = F'
		 * @param  f       [in] RHS vector (size m_eq) for equalities
		 * @param  loqo_lp [in/out] Pointer (void*) to LOQO data structure.  On input
		 *                 loqo_lp->n and oqo_lp->m should already be set.  On output
		 *                 loqo_lp->qnz, loqo_lp->Q, loqo_lp->iQ and loqo_lp->kQ should be
		 *                 set for Q and loqo_lp->nz, loqo_lp->A, loqo_lp->iA and loqo_lp->kA
		 *                 should be set for A.  The reason a void pointer is passed is
		 *                 because the declarations in loqo.h conflicit with may things
		 *                 because of bad software engineering.
		 */
		virtual void init_hess_jacob(
			const MatrixWithOp& G, const value_type bigM
			, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
			, const int loqo_b_stat[], const size_type num_inequal
			, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
			, void* loqo_lp
			) const;

	}; // end class InitLOQOHessianJacobian

	///
	/** Strategy object that sets up Q and A for LOQO
	  */
	STANDARD_COMPOSITION_MEMBERS( InitLOQOHessianJacobian, init_hess_jacob )

 	///
	/** <<std member comp>> Big M parameter used in the objective.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, bigM )

 	///
	/** <<std member comp>> Tolerance under which inequalities are
	 * considered nonbinding.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, nonbinding_lag_mult )

	///
	QPSolverRelaxedLOQO(
		const init_hess_jacob_ptr_t    init_hess_jacob      = new  InitLOQOHessianJacobian()
		,value_type                     bigM                 = 1e+10
		,value_type                     nonbinding_lag_mult  = 1e-12
		);

	///
	~QPSolverRelaxedLOQO();

	// /////////////////////////////////
	// Overridden from QPSolverRelaxed

	///
	QPSolverStats get_qp_stats() const;

	///
	void release_memory();

protected:

	// /////////////////////////////////
	// Overridden from QPSolverRelaxed

	///
	QPSolverStats::ESolutionType imp_solve_qp(
		  std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
			, const SpVectorSlice* eL, const SpVectorSlice* eU
		, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
		, value_type* obj_d
		, value_type* eta, VectorSlice* d
		, SpVector* nu
		, SpVector* mu, VectorSlice* Ed
		, VectorSlice* lambda, VectorSlice* Fd
		);

private:

	// //////////////////////////////////////////////////////////////
	// Private types

	// //////////////////////////////////////////////////////////////
	// Private Data Members.

	QPSolverStats			qp_stats_;

	// ////////////////////////////
	// Private member functions

};	// end class QPSolverRelaxedLOQO

}	// end namespace ConstrainedOptimizationPackTypes

#endif // USE_QP_LOQO

#endif	// QP_SOLVER_RELAXED_QPKWIK_H