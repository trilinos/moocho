// //////////////////////////////////////////////////////////////
// ConstraintsRelaxedStd.h
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

#ifndef QP_SCHUR_CONSTRAINTS_RELAXED_STD_H
#define QP_SCHUR_CONSTRAINTS_RELAXED_STD_H

#include <list>

#include "QPSchur.h"
#include "AbstractLinAlgPack/include/MatrixWithOp.h"
#include "AbstractLinAlgPack/include/VectorSpaceCompositeStd.h"

namespace ConstrainedOptimizationPack {
namespace QPSchurPack {

///
/** Constraints subclass that is used to represent generic
 * varaible bounds, and general inequality and equality constraints.
 *
 * The generic constraints represented by this class are those
 * of the \c QPSolverRelaxed interface which are:
 \verbatim
  
  (1.2)               etaL <=  eta
  (1.3)               dL   <=  d                                 <= dU
  (1.4)               eL   <=  op(E)*d - b*eta                   <= eU
  (1.5)                        P_u'*op(F)*d + (1 - eta) * P_u'*f  = 0
 \endverbatim
 * These constraints are converted into the form:
 \verbatim

       [   dL    ]     [ I                    ]              [   dU    ]
       [   etaL  ] <=  [                   1  ] * [  d  ] <= [   inf   ]
  (2)  [   eL    ]     [ op(E)            -b  ]   [ eta ]    [   eU    ]
       [ -P_u'*f ]     [ P_u'*op(F)   -P_u'*f ]              [ -P_u'*f ]
       \_________/     \______________________/   \_____/    \_________/
         cL_bar                 A_bar'               x          cU_bar

       =>

  (3)   [     xl     ]    [  I        ]          [     xu     ]
        [  cl_breve  ] <= [  A_breve' ] * x  <=  [  cu_breve  ]

       =>

  (4)	cl_bar <= A_bar'*x <= cu_bar
 \endverbatim
 * The main responsibilities of this class are to expose a
 * \c MatrixWithOp object for \c A_bar shown in (2) and to pick
 * violated constraints.
 */
class ConstraintsRelaxedStd : public Constraints {
public:

	// /////////////////////////////////////////////
	// Public types

	///
	/** Matrix type for A_bar.
	 *
	 \verbatim

		A_bar = [  I   0  op(E')   op(F')*P_u  ]
		        [  0   1   -b'       -f'*P_u   ]

	 \endverbatim
	 *
	 */
	class MatrixConstraints : public MatrixWithOp {
	public:

		///
		/** Construct to unitinitialized.
		 *
		 * this->nd() == 0 after construction.
		 */
		MatrixConstraints();

		///
		/** Initialize.
		 *
		 * The sizes of the arguments are
		 * not checked.
		 * 
		 * It is expected that the objects being
		 * pointed to will not be resized or invalidated
		 * since copies of data are not made!
		 */
		void initialize(
			const VectorSpace::space_ptr_t   &space_d_eta
			,const size_type                 m_in
			,const size_type                 m_eq
			,const MatrixWithOp              *E
			,BLAS_Cpp::Transp                trans_E
			,const VectorWithOp              *b
			,const MatrixWithOp              *F
			,BLAS_Cpp::Transp                trans_F
			,const VectorWithOp              *f
			,size_type                       m_undecomp
			,const size_type                 j_f_undecomp[]
			);

		/** @name Access */
		//@{

		///
		size_type			nd() const
		{	return nd_;	}
		///
		size_type			m_in() const
		{	return m_in_;	}
		///
		size_type			m_eq() const
		{	return m_eq_;	}
		///
		const MatrixWithOp*	E() const
		{	return E_;	}
		///
		BLAS_Cpp::Transp	trans_E() const
		{	return trans_E_;	}
		///
		const VectorWithOp*	b() const
		{	return b_;	}
		///
		const MatrixWithOp*	F() const
		{	return F_;	}
		///
		BLAS_Cpp::Transp	trans_F() const
		{	return trans_F_;	}
		///
		const VectorWithOp*	f() const
		{	return f_;	}
		///
		const GenPermMatrixSlice& P_u() const
		{	return P_u_;	}

		//@}

		/** @name Overridden from MatrixWithOp */
		//@{

		///
		const VectorSpace& space_cols() const;
		///
		const VectorSpace& space_rows() const;
		///
		MatrixWithOp& operator=(const MatrixWithOp& m);
//		///
//		void Mp_StPtMtP(
//			GenMatrixSlice* gms_lhs, value_type alpha
//			,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
//			,BLAS_Cpp::Transp M_trans
//			,const GenPermMatrixSlice& P_rhs2, BLAS_Cpp::Transp P_rhs2_trans
//			) const ;
		///
		void Vp_StMtV(
			VectorWithOpMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
			,const VectorWithOp& vs_rhs2, value_type beta
			) const;
		///
		void Vp_StPtMtV(
			VectorWithOpMutable* vs_lhs, value_type alpha
			,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
			,BLAS_Cpp::Transp M_rhs2_trans
			,const SpVectorSlice& sv_rhs3, value_type beta
			) const;

		//@}
		
	private:
		typedef std::vector<size_type>		row_i_t;
		typedef std::vector<size_type>		col_j_t;
		size_type			nd_;	// # unknowns d
		size_type			m_in_;	// # op(E)*d inequality constraints
		size_type			m_eq_;	// # op(F)*d equality constraints
		const MatrixWithOp	*E_;	// If NULL then no general inequalities
		BLAS_Cpp::Transp	trans_E_;
		const VectorWithOp	*b_;
		const MatrixWithOp	*F_;	// If NULL then no general equalities
		BLAS_Cpp::Transp	trans_F_;
		const VectorWithOp	*f_;
		GenPermMatrixSlice  P_u_;
		row_i_t             P_u_row_i_;
		col_j_t             P_u_col_j_;
		VectorSpace::space_ptr_t   space_cols_;
		VectorSpaceCompositeStd    space_rows_;
	};	// end class MatrixConstraints

	///
	enum EInequalityPickPolicy {
		ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY
		,ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY
		,ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY
	};

	/** @name Public member functions */
	//@{

	///
	/** <<std comp>> members for feasibility tolerance for the bound constriants.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, bounds_tol )

	///
	/** <<std comp>> members for feasibility tolerance for the general inequality constraints.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, inequality_tol )

	///
	/** <<std comp>> members for feasibility tolerance for the general equality constriants.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, equality_tol )

	///
	/** <<std comp>> members for policy used to select a violated constraint.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EInequalityPickPolicy, inequality_pick_policy )

	/// Constructs to uninitialized
	ConstraintsRelaxedStd();

	///
	/** Initialize constriants.
	 *
	 * If there are no variable bounds then set:<br>
	 * <tt>void(dL) == void(dU) == NULL</tt>
	 * 
	 * If there are no general inequality constraints
	 * then set:<br>
	 * <tt>void(E) == void(b) == void(eL) == void(eU) == NULL</tt>
	 * 
	 * If there are no general equality constraints then
	 * set:<br>
	 * <tt>void(F) = void(f) == NULL</tt>
	 * 
	 * If <tt>check_F == false</tt>, then the equality constriants
	 * in <tt>op(F)</tt> will not be checked as violated constriants.
	 * This is to facilitate the addition of the equality
	 * constraints to the initial schur complement and therefore
	 * these constraints should never be violated (except for
	 * illconditioning).
	 * The tolerances below which a constriant will not be considered
	 * violated are given by <tt>bounds_tol</tt>, <tt>inequality_tol</tt> and <tt>equality_tol</tt>.
	 * 
	 * Here, <tt>Ed</tt> is updated (if <tt>Ed != NULL</tt>) within the function
	 * <tt>this->pick_violated}(...)</tt>.  This saves some computational work of
	 * having to compute <tt>op(E)*d</tt> again.  To skip computing this value, just set
	 * <tt>Ed == NULL</tt>.
	 *
	 * ToDo: Specify more concretely exactly what the criteria is for
	 * considering that a constraint is violated or in picking the most
	 * violated constraint.
	 *
	 * @param  m_undecomp
	 *                  [in] Number of undecomposed equality constraints.
	 * @param  j_f_undecomp
	 *                  [in] array (size m_undecomp) of indexes of constraints
	 *                  in op(F)*d + (1-eta)*f that are not decomposed and therefore
	 *                  should be considered when looking for violated constraints.
	 *                  This array is used to define the mapping matrix P_u.
	 *                  It is required that this be sorted and that:
	 *                  j_f_undecomp[k+1] >= j_f_undecomp[k], for k = 0...m_undecomp-2.
	 *                  If m_undecomp == f->size() then j_f_undecomp == NULL is allowed
	 *                  and the matrix P_u will be the identity matrix.
	 */
	void initialize(
		const VectorSpace::space_ptr_t   &space_d_eta
		,value_type                      etaL
		,const VectorWithOp              *dL
		,const VectorWithOp              *dU
		,const MatrixWithOp              *E
		,BLAS_Cpp::Transp                trans_E
		,const VectorWithOp              *b
		,const VectorWithOp              *eL
		,const VectorWithOp              *eU
		,const MatrixWithOp              *F
		,BLAS_Cpp::Transp                trans_F
		,const VectorWithOp              *f
		,size_type                       m_undecomp
		,const size_type                 j_f_undecomp[]
		,VectorWithOpMutable             *Ed
		,bool                            check_F           = true
		,value_type                      bounds_tol        = 1e-10
		,value_type                      inequality_tol    = 1e-10
		,value_type                      equality_tol      = 1e-10
		);

	///
	const MatrixConstraints& A_bar_relaxed() const;

	//@}

	/** @name Overridden from Constraints */
	//@{

	///
	size_type n() const;
	///
	size_type m_breve() const;
	///
	/** Represents the constraints matrix.
	 *
	 \verbatim

		A_bar = [  I   0  op(E')   op(F')*P_u  ]
		        [  0   1   -b'       -f'*P_u   ]
	 \endverbatim
	 */
	const MatrixWithOp& A_bar() const;
	///
	void pick_violated_policy( EPickPolicy pick_policy );
	///
	EPickPolicy pick_violated_policy() const;
	///
	/** Here the next violated constraint to add to the active set is selected.
	 *
	 * Violated constraints are selected to to add to the active set in the following
	 * order:
	 * <ul>
	 * <li> The equality constraints are added first, one at a time (if not already added
	 *		as part of the warm start).
	 * <li> Add inequality constraints according according to the following options:
	 *		<ul>   
	 *		<li> <tt>ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY</tt>
	 *			Check the variable bounds first and add the most violated.  If no
	 *			variable bounds are violated by more than <tt>this->bounds_tol()</tt> then check for
	 *			the most violated inequality constraint by computing <tt>r = op(E)*d+b*eta</tt> and
	 *			add the most violated bound (<tt>eL</tt>, <tt>eU</tt>) if one exists.
	 *		<li> <tt>ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY</tt>
	 *			Check the variable bounds first and add the most violated.  If no
	 *			variable bounds are violated by more than <tt>this->bounds_tol()</tt> then check for
	 *			the first violated inequality constraint by computing <tt>e(j)'*(op(E)*d+b*eta)</tt>
	 *			one or more constraints at a time.  This option may be
	 *			better if the cost of computing <tt>op(E)*d</tt> is significant.
	 *		<li> <tt>ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY</tt>
	 *			Select the most violated constraint from the variable bounds and the
	 *			general inequality constraints by computing  r = op(E)*d+b*eta then
	 *			add the most violated variable bound.  This option is always the most
	 *			expensive but may result in less QP iterations.
	 * </ul>
	 *
	 * As a side effect, the vector pointed to by <tt>Ed</tt> which was passed to
	 * <tt>this->initialize(...)</tt> will be guarrenteed to be updated for
	 * the current </tt>op(E)*d</tt>, where </tt>d = x(1,nd)</tt>, if any of the following is true:
	 * <ul>
	 * <li> <tt>j_viol == 0</tt>
	 * <li> <tt>this->pick_violated_policy() == ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY</tt>
	 * <li> <tt>this->pick_violated_policy() == ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY</tt>
	 *		<tt>&& j_viol > this->n()</tt>
	 *	</ul>
	 * If none of the above criteria applies then the client can not assume that
	 * <tt>Ed</tt> was updated and therefore the client must compute this value on its own.
	 */
	void pick_violated(
		const VectorSlice& x, size_type* j_viol, value_type* constr_val
		,value_type* viol_bnd_val, value_type* norm_2_constr, EBounds* bnd, bool* can_ignore
		) const;
	///
	void ignore( size_type j );
	///
	value_type get_bnd( size_type j, EBounds bnd ) const;

	//@}

private:

	// //////////////////////////////
	// Private types

	typedef std::list<size_type>  passed_over_equalities_t;

	// //////////////////////////////
	// Private data members

	MatrixConstraints   A_bar_;
	value_type          etaL_;
	const VectorWithOp  *dL_;	// If NULL then no simple bounds
	const VectorWithOp  *dU_;
	const VectorWithOp  *eL_;
	const VectorWithOp  *eU_;
	VectorWithOpMutable *Ed_;
	bool                check_F_;
	mutable size_type   last_added_j_;          // Remember the last bound added so that
	mutable value_type  last_added_bound_;      // we can save our selfs some work.
	mutable EBounds     last_added_bound_type_; // ...
	mutable size_type   next_undecomp_f_k_;
	   	// Determines the next constraint [P_u'*op(F)*d + (1 - eta) * P_u'*f](next_undecomp_f_k)
	    // to be checked to see if it is violated.  If next_undecomp_f_k > P_u.nz() then all
	    // of the constriants have been checked.
	mutable passed_over_equalities_t passed_over_equalities_;
	    // This is a list that keeps track of those equality constraints that were checked
	    // for being violated but were within tolerance and therefore passed over and not added.
	    // This list can be traversed again and again to check these constraints.  Specifically, the
	    // indexes of f(k) are sorted, not the indexes in P_u'.

	// //////////////////////////////
	// Private member functions

	///
	void cache_last_added(
		size_type last_added_j, value_type last_added_bound
		,EBounds last_added_bound_type
		) const;

}; // end class ConstraintsRelaxedStd

} // end namespace QPSchurPack 
} // end namespace ConstrainedOptimizationPack 

#endif // QP_SCHUR_CONSTRAINTS_RELAXED_STD_H
