// //////////////////////////////////////////////////////////
// QPSolverRelaxedQPOPT.cpp
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

#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT

#include <assert.h>

#include <algorithm>
#include <ostream>
#include <iomanip>

#include "ConstrainedOptimizationPack/include/QPSolverRelaxedQPOPT.h"
#include "ConstrainedOptimizationPack/include/QPOPT_CppDecl.h"
#include "SparseLinAlgPack/include/MatrixWithOp.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "SparseLinAlgPack/include/EtaVector.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/GenMatrixOut.h"
#include "LinAlgPack/include/VectorOut.h"
#include "Misc/include/debug.h"


// //////////////////////////////////////////////////////////
// Local implementation functions.

namespace {

typedef FortranTypes::f_int			f_int;
typedef FortranTypes::f_dbl_prec	f_dbl_prec;
typedef FortranTypes::f_logical		f_logical;

// Compute:
//
// HESS * x = [ G  0 ] * [ X(1,N-1) ] = [ G * X(1,N-1) ]
//            [ 0  M ]   [   X(N)   ]   [ M * X(N)     ]
//
// The matrix vector product is implemented through the MatrixWithOp interface.
//
inline
void qphess_server_relax( const f_int& N, const f_int& LDH
	, const f_int& JTHCOL, const f_dbl_prec* H, const f_dbl_prec* X, f_dbl_prec* HX
	, f_int* IW, const f_int& LENIW, f_dbl_prec* W, const f_int& LENW )
{
	using LinAlgPack::VectorSlice;
	using SparseLinAlgPack::SpVector;
	using LinAlgOpPack::V_MtV;
	using ConstrainedOptimizationPack::QPSolverRelaxedQPOPT;

	// Here we have used some casting to pass on information about the qp solver
	// that called QPSOL.
	const QPSolverRelaxedQPOPT* qp_solver = reinterpret_cast<const QPSolverRelaxedQPOPT*>(H);

	VectorSlice hx(HX,N);

	if( JTHCOL == 0 ) {
		const VectorSlice x( const_cast<VectorSlice::value_type*>(X), N );
		// hx(1,N-1) = G * x(1,N-1)
		V_MtV( &hx(1,N-1), *qp_solver->G(), BLAS_Cpp::no_trans, x(1,N-1) );
		// hx(N) = bigM * x(N)
		hx(N) = qp_solver->use_as_bigM() * x(N);
	}
	else {
		// we are extracting the JTHCOL column of G so use sparse x
		if(JTHCOL == N) {
			// 0
			std::fill( HX, HX + (N-1), 0.0 );
			// bigM
			HX[N-1] = qp_solver->use_as_bigM();
		}
		else {
			// G(:,JTHCOL)
			SparseLinAlgPack::EtaVector e_j(JTHCOL,N-1);
			V_MtV( &hx(1,N-1), *qp_solver->G(), BLAS_Cpp::no_trans, e_j() );
			// 0
			hx(N) = 0.0;
		}
	}
}

}	// end namespace

// ///////////////////////////////////////////////////////////////////////////
// Fortran declarations.

extern "C" {

// These are declarations for the subroutines that preform the communication
// between C++ and Fortran.  There is no use in putting them in a
// namespace since the namespace name will not be used by the linker since
// we are using extern "C".

//
FORTRAN_FUNC_DECL_UL( void, QPHESS_SERVER_RELAX2, qphess_server_relax2 ) ( const f_int& N, const f_int& LDH
	, const f_int& JTHCOL, const f_dbl_prec* H, const f_dbl_prec* X, f_dbl_prec* HX
	, f_int* IW, const f_int& LENIW, f_dbl_prec* W, const f_int& LENW )
{
	qphess_server_relax( N, LDH, JTHCOL, H, X, HX, IW, LENIW, W, LENW );
}

}	// end extern "C"

namespace LinAlgOpPack {
	using SparseLinAlgPack::Mp_StM;
	using SparseLinAlgPack::Vp_StMtV;
}

// ///////////////////////////////////////
// QPSolverRelaxedQPOPT

namespace ConstrainedOptimizationPack {

QPSolverRelaxedQPOPT::QPSolverRelaxedQPOPT(
	value_type        max_qp_iter_frac
	)
	:max_qp_iter_frac_(max_qp_iter_frac)
	,ITMAX_(100)
	,BIGBND_(std::numeric_limits<value_type>::max())
	,FEATOL_( 1.0e-10 )
	,LDH_(1)
{}

QPSolverRelaxedQPOPT::~QPSolverRelaxedQPOPT()
{
	release_memory();
}

// Overridden from QPSolverRelaxed

void QPSolverRelaxedQPOPT::release_memory()
{
	// ToDo: Resize all arrays to zero!
	QPSolverRelaxedQPOPTSOL::release_memory();
}

// Overridden protected members

QPSolverRelaxedQPOPT::f_int QPSolverRelaxedQPOPT::liwork(f_int N, f_int NCLIN) const
{	return 2* N + 3; }

QPSolverRelaxedQPOPT::f_int QPSolverRelaxedQPOPT::lrwork(f_int N, f_int NCLIN) const
{	return NCLIN > 0 ? 2*N*N + 8*N + 5*NCLIN : N*N + 8 *N; }

QPSolverRelaxedQPOPTSOL::EInform QPSolverRelaxedQPOPT::call_qp_solver(bool warm_start)
{
	ITMAX_ = max_qp_iter_frac() * N_;

	// Set option parameters
	{
		namespace ns = QPOPT_CppDecl;
		namespace ft = FortranTypes;
		ns::reset_defaults();
		ns::set_logical_option( ns::WARM_START				, warm_start ? ft::FALSE : ft::TRUE		);
		ns::set_real_option(    ns::FEASIBILITY_TOLERANCE	, FEATOL_								);
		ns::set_real_option(    ns::INFINITE_BOUND_SIZE		, BIGBND_								);
		ns::set_int_option(     ns::ITERATION_LIMIT			, ITMAX_								);
		ns::set_int_option(     ns::PRINT_FILE				, 0										);
		ns::set_int_option(     ns::SUMMARY_FILE			, 0										);
		ns::set_int_option(     ns::PRINT_LEVEL				, 0										);
		ns::set_int_option(     ns::PROBLEM_TYPE			, ns::QP2								);
	}
	
	// Set the rest of the QPOPT inputs that could not be set in the constructor.

	LDA_	= ( A_.rows() > 0 ? A_.rows() : 1 );
	H_		= reinterpret_cast<f_dbl_prec*>(this);	// used to implement QPHESS
	AX_.resize( NCLIN_ > 0 ? NCLIN_ : 1 );

	QPOPT_CppDecl::qpopt(
		N_, NCLIN_, LDA_, LDH_, &A_(1,1), &BL_(1), &BU_(1)
		, &CVEC_(1), H_
		, FORTRAN_NAME_UL(QPHESS_SERVER_RELAX2,qphess_server_relax2)
		, &ISTATE_[0], &X_(1), INFORM_, ITER_, OBJ_, &AX_(1)
		, &CLAMDA_(1), &IWORK_[0], LIWORK_, &WORK_[0], LWORK_ );

	EInform return_inform;
	typedef QPSolverRelaxedQPOPTSOL bc;
	switch(INFORM_) {
	    case STRONG_LOCAL_MIN:
			return_inform = bc::STRONG_LOCAL_MIN;
			break;
	    case WEAK_LOCAL_MIN:
			return_inform = bc::WEAK_LOCAL_MIN;
			break;
	    case UNBOUNDED:
			throw Unbounded(
				"QPSolverRelaxedQPOPT::call_qp_solver() : Error,"
				" QPOPT returned that the QP is unbounded!" );
	    case INFEASIBLE:
			throw Infeasible(
				"QPSolverRelaxedQPOPT::call_qp_solver() : Error,"
				" QPOPT returned that the QP is infeasible!" );
	    case ITMAX_EXCEEDED:
			return_inform = bc::MAX_ITER_EXCEEDED;
			break;
	    case MAX_DOF_TOO_SMALL:
			throw std::runtime_error(
				"QPSolverRelaxedQPOPT::call_qp_solver() : Error,"
				" QPOPT says that the max dof is too small" );
	    case INVALID_INPUT:
			throw InvalidInput(
				"QPSolverRelaxedQPOPT::call_qp_solver() : Error,"
				" QPOPT returned that the input is invalid" );
	    case PROB_TYPE_NOT_REGOG:
			throw std::logic_error(
				"QPSolverRelaxedQPOPT::call_qp_solver() : Error,"
				" QPOPT says that the problem type is not recognized" );
			break;
	    default:
			assert(0); // Should not happen
	}

	return return_inform;
}

} // end namespace ConstrainedOptimizationPack

#endif // CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT
