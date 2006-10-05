// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_LOQO

#ifndef QP_SOLVER_RELAXED_LOQO_H
#define QP_SOLVER_RELAXED_LOQO_H

#include <vector>

#include "ConstrainedOptPack_QPSolverRelaxed.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Solves Quadratic Programming (QP) problem using the primal-dual interior-point
  * solver LOQO (Vanderbei).
  *
  * In this implementation it is required that G support the \Ref{MatrixExtractInvCholFactor}
  * interface and is therefore quite restrictive on the likes of QPs it can solve.
  */
class QPSolverRelaxedLOQO : public QPSolverRelaxed
{
public:

  /** \brief Strategy interface for initializing the Hessian matrix Q and the
   * Jacobian matrix A for LOQO.
   *
   * This interface allows clients to specialize how G, op(E), and op(F)
   * are translated into the sparse compresssed column data structurs for
   * Q and A for LOQO.
   */
  class InitLOQOHessianJacobian {
  public:

    /** \brief . */
    virtual ~InitLOQOHessianJacobian() {}

    /** \brief Extract the Hessian Q and Jacobian A matrices.
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
      const MatrixOp& G, const value_type bigM
      , const MatrixOp* E, BLAS_Cpp::Transp trans_E, const DVectorSlice* b
      , const int loqo_b_stat[], const size_type num_inequal
      , const MatrixOp* F, BLAS_Cpp::Transp trans_F, const DVectorSlice* f
      , void* loqo_lp
      ) const;

  }; // end class InitLOQOHessianJacobian

  /** \brief Strategy object that sets up Q and A for LOQO
    */
  STANDARD_COMPOSITION_MEMBERS( InitLOQOHessianJacobian, init_hess_jacob )

  /** \brief <<std member comp>> Big M parameter used in the objective.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, bigM )

  /** \brief <<std member comp>> Tolerance under which inequalities are
   * considered nonbinding.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, nonbinding_lag_mult )

  /** \brief . */
  QPSolverRelaxedLOQO(
    const init_hess_jacob_ptr_t    init_hess_jacob      = new  InitLOQOHessianJacobian()
    ,value_type                     bigM                 = 1e+10
    ,value_type                     nonbinding_lag_mult  = 1e-12
    );

  /** \brief . */
  ~QPSolverRelaxedLOQO();

  // /////////////////////////////////
  // Overridden from QPSolverRelaxed

  /** \brief . */
  QPSolverStats get_qp_stats() const;

  /** \brief . */
  void release_memory();

protected:

  // /////////////////////////////////
  // Overridden from QPSolverRelaxed

  /** \brief . */
  QPSolverStats::ESolutionType imp_solve_qp(
      std::ostream* out, EOutputLevel olevel, ERunTests test_what
    , const DVectorSlice& g, const MatrixOp& G
    , value_type etaL
    , const SpVectorSlice& dL, const SpVectorSlice& dU
    , const MatrixOp* E, BLAS_Cpp::Transp trans_E, const DVectorSlice* b
      , const SpVectorSlice* eL, const SpVectorSlice* eU
    , const MatrixOp* F, BLAS_Cpp::Transp trans_F, const DVectorSlice* f
    , value_type* obj_d
    , value_type* eta, DVectorSlice* d
    , SpVector* nu
    , SpVector* mu, DVectorSlice* Ed
    , DVectorSlice* lambda, DVectorSlice* Fd
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

#endif	// QP_SOLVER_RELAXED_QPKWIK_H

#endif // CONSTRAINED_OPTIMIZATION_PACK_USE_LOQO
