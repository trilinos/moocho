// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproachOrthogonal_Step.cpp

#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproachOrthogonal_Step.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/GenMatrixOp.h"
#include "LinAlgPack/include/GenMatrixOut.h"
#include "LinAlgPack/include/assert_print_nan_inf.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgLAPack/include/LinAlgLAPack.h"

namespace ReducedSpaceSQPPack {

EvalNewPointTailoredApproachOrthogonal_Step::EvalNewPointTailoredApproachOrthogonal_Step(
		  const deriv_tester_ptr_t& 	deriv_tester
		, EFDDerivTesting				fd_deriv_testing
		)
	:
		EvalNewPointTailoredApproach_Step(deriv_tester,fd_deriv_testing)
{}
// protected

void EvalNewPointTailoredApproachOrthogonal_Step::calc_py_Ypy(
	  const GenMatrixSlice& D, VectorSlice* py, Vector* Ypy
	, EJournalOutputLevel olevel, std::ostream& out
	)
{
	const size_type
		n = D.rows()+D.cols(),
		r = D.cols(),
		dof = n-r;

	//
	// First let's define:
	//
	// Gc(decomp) = [ C, N ]
	// D = -inv(C)*N
	// Y = [ I ; -D' ]
	//
	// R = Gc(decomp)'*Y = C*(I + D*D')
	// 
	// Now we want to compute:
	// 
	// py = - inv(R) * c
	// 
	// To do so let's see what inv(R) is:
	// 
	// inv(R) = inv(I + D*D') * inv(C)
	// 
	// From the Sherman-Morrison-Woodbury Formula (Nocedal&Wright 1999, p. 605)
	// 
	// inv(I+D*D') = inv(I) - inv(I) * D * inv( I + D' * inv(I) * D ) * D' * inv(I)
	//             = I - D * inv(S) * D'
	//     where: S = I + D'*D
	//
	// This gives us:
	// 
	// inv(R) = ( I - D * inv(S) * D' ) * inv(C)
	// 
	// Therefore we need to compute and factorize S in order to factorize
	// R.
	// 
	// Therefore, we can compute py as:
	// 
	// py = -inv(R)*c
	//    = (I - D * inv(S) * D' ) * inv(C) * (-c)
	//                              \___________/
	//                               py on input
	//                               
	//    = py - D * inv(S) * D' * py
	//    
	// And  we can compute this as:
	// 
	// t1 = D'*py
	// 
	// t2 = inv(S) * t1
	// 
	// py = py - D * t2
	//

	// Form the lower triangular part of S = I + D'*D
	using BLAS_Cpp::lower;
	using BLAS_Cpp::nonunit;
	using LinAlgPack::nonconst_sym;
	using LinAlgPack::nonconst_tri_ele;
	GenMatrix S(dof,dof);
	S = 0.0;
	S.diag() = 1.0;
	LinAlgPack::syrk( BLAS_Cpp::no_trans, 1.0, D, 1.0
		, &nonconst_sym( S(), lower ) );
	// Factor the lower triangular part of S = I + D*D'
	LinAlgLAPack::potrf( &nonconst_tri_ele( S(), lower ) );

	// Now we have the cholesky factors fo S = L * L'
	tri_gms L( S(), lower, nonunit );

	// Now solve for py = py - D * inv(S) * D' * py

	// t = D'*py
	Vector t;
	LinAlgOpPack::V_MtV( &t, D, BLAS_Cpp::trans, *py );

	// t = inv(S) * t
	//    = inv(L*L') * t
	//    = inv(L')*inv(L)*t
	//    
	// Note that another temporary is not needed since
	// a direct call to the BLAS is used and aliaseing the
	// rhs and lhs is fine.
	// 
	LinAlgPack::V_invMtV( &t(), L, BLAS_Cpp::no_trans, t() );
	LinAlgPack::V_invMtV( &t(), L, BLAS_Cpp::trans, t() ); 	
	
	// py = py - D * t
	LinAlgPack::Vp_StMtV( py, -1.0, D, BLAS_Cpp::no_trans, t() );
	
	// Now all we have to do is to compute Ypy
	// 
	// Ypy = Y * py = [  I  ] py = [  py    ]
	//                [ -D' ]      [ -D'*py ]
	//                
	Ypy->resize(n);
	(*Ypy)(1,r) = *py;
	LinAlgOpPack::V_StMtV( &(*Ypy)(r+1,n), -1.0, D, BLAS_Cpp::trans, *py );

	// That's it man!
}

void EvalNewPointTailoredApproachOrthogonal_Step::print_calc_Y_py_Ypy(
	std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Orthogonal decomposition\n"
		<< L << "py_k = inv(I + D*D') * py_k\n"
		<< L << "Y = [ I ; -D' ] <: R^(n x m)   [Not computed explicity]\n"
		<< L << "Ypy_k = Y * py_k\n"
		;
}

}	// end namespace ReducedSpaceSQPPack 