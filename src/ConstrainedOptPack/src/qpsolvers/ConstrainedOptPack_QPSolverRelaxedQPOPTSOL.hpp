// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxedQPOPTSOL.hpp
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

#ifndef QP_SOLVER_RELAXED_QPOPTSOL_H
#define QP_SOLVER_RELAXED_QPOPTSOL_H

#include <vector>

#include "QPSolverRelaxed.hpp"
#include "DenseLinAlgPack/src/DVectorClass.hpp"
#include "DenseLinAlgPack/src/DMatrixClass.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

///
/** Node base clase for the primal QP solvers QPOPT and QPSOL.
  *
  * In this implementation it is required that G only support the \Ref{MatrixOp}
  * interface and is therefore quite flexible in the QPs it can solve.
  */
class QPSolverRelaxedQPOPTSOL : public QPSolverRelaxed
{
public:

	// /////////////////////////////////////
	/** @name Public Types */
	//@{

	///
	typedef FortranTypes::f_int			f_int;
	///
	typedef FortranTypes::f_dbl_prec	f_dbl_prec;
	///
	typedef FortranTypes::f_logical		f_logical;

	//@}

	///
	QPSolverRelaxedQPOPTSOL();

	///
	~QPSolverRelaxedQPOPTSOL();

	/// Return a pointer to the matrix G to be used in the calculation of H*x by QPOPT and QPSOL.
	virtual const MatrixOp* G() const;

	/// Return the value of the "big M" used in the relaxation (called by QPHESS functions).
	virtual value_type use_as_bigM() const;

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

	// //////////////////////////////////////////////////////////////
	// Protected types

	///
	typedef std::vector<f_int>			ISTATE_t;
	///
	typedef std::vector<f_int>			IWORK_t;
	///
	typedef std::vector<f_dbl_prec>		WORK_t;
public: // RAB: 2001/05/03: MS VC++ 6.0 must have this public ???
	///
	enum EInform {
		STRONG_LOCAL_MIN,
		WEAK_LOCAL_MIN,
		MAX_ITER_EXCEEDED,
		OTHER_ERROR
	};

protected:

	// //////////////////////////////////////////////////////////////
	// Protected Data Members.

	QPSolverStats			qp_stats_;

	/** @name Input/Output parameters common to both QPOPT and QPSOL
	  *
	  * These are access and updated by subclasses that call QPOPT
	  * and QPSOL.
	  */
	//@{

	///
	f_int		N_;
	///
	f_int		NCLIN_;
	///
	DMatrix	A_;
	///
	DVector		BL_;
	///
	DVector		BU_;
	///
	DVector		CVEC_;
	///
	ISTATE_t	ISTATE_;
	///
	DVector		X_;
	///
	DVector		AX_;
	///
	DVector		CLAMDA_;
	///
	f_int		ITER_;
	///
	f_dbl_prec	OBJ_;
	///
	f_int		LIWORK_;
	///
	IWORK_t		IWORK_;
	///
	f_int		LWORK_;
	///
	WORK_t		WORK_;

	//@}

	// /////////////////////////////////////////////////////////////
	// Template method primatives to be overridden.

	/// Length of integer workspace
	virtual f_int liwork(f_int N, f_int NCLIN) const = 0;

	/// Length of real workspace
	virtual f_int lrwork(f_int N, f_int NCLIN) const = 0;

	///
	/** Solve the QP defined in the protected input data members
	  * and set the solution in the protected output data members.
	  */
	virtual EInform call_qp_solver(bool warm_start) = 0;

private:

	// /////////////////////////
	// Private types

	typedef std::vector<f_int>	ibnds_t;

	// ///////////////////////////
	// Private data members

	size_type			n_inequ_bnds_;		// Used to record the number of bounds with at least
											// one bound existing in eL and eU.
	ibnds_t				i_inequ_bnds_;		// size(nbounds_). Remembers which bounds in
											// eL, eU had at least one bound present.  This is
											// needed to map from CLAMDA_ to mu.
	value_type			bigM_;				// Big M value used to construct relaxation.
	value_type			use_as_bigM_;		// Big M value used in QPHESS.
	const MatrixOp*	G_;					// used to compute HESS * x = [ G, 0; 0, bigM ] * x products.

	// ///////////////////////////
	// Private member functions

	// not defined and not to be called.
	QPSolverRelaxedQPOPTSOL(const QPSolverRelaxedQPOPTSOL&);
	QPSolverRelaxedQPOPTSOL& operator=(const QPSolverRelaxedQPOPTSOL&);

};	// end class QPSolverRelaxedQPOPTSOL

}	// end namespace ConstrainedOptimizationPackTypes

#endif	// QP_SOLVER_RELAXED_QPOPTSOL_H
