// //////////////////////////////////////
// MatrixKKTFullSpaceRelaxed.hpp
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

#ifndef MATRIX_KKT_FULL_SPACE_RELAXED_H
#define MATRIX_KKT_FULL_SPACE_RELAXED_H

#include "ConstrainedOptimizationPackTypes.hpp"
#include "SparseLinAlgPack/src/MatrixWithOpFactorized.hpp"
#include "SparseLinAlgPack/src/MatrixConvertToSparseFortranCompatible.hpp"
#include "StandardCompositionMacros.hpp"

namespace ConstrainedOptimizationPack {

///
/** Implementation of a KKT matrix factorized in the full space.
  *
  * This class is used to represent the KKT matrix of the following
  * relaxed QP:
  *
  \begin{verbatim}
	min     [ g'  M ] * [  d  ] + 1/2 * [ d'  eta ] * [ G      ] * [  d  ]
	                    [ eta ]                       [     M  ]   [ eta ]

	s.t.    [ A'  -c ] * [  d  ] + c = 0
	                     [ eta ]
  \end{verbatim}
  *
  * The only matrix actually factorized is:
  *
  \begin{verbatim}
	K_bar = [ G  A ]
	        [ A'   ]
  \end{verbatim}
  *
  * The class has two modes.
  *
  * First mode is to not include the relaxation
  * term and therefore the KKT matrix is:
  *
  \begin{verbatim}
	K = [ G  A ]
	    [ A'   ]
  \end{verbatim}
  *
  * The second mode is the use the relaxation and he represented matrix is:
  *
  \begin{verbatim}
	    [ G       A   ]
	K = [     M   -c' ]
	    [ A'  -c      ]
  \end{verbatim}
  *
  * This class uses an aggregate DirectSparseFortranCompatibleSolver (DSFCS) 
  * object to factorize K above and then to solve for the linear systems
  * involving K.
  */
class MatrixKKTFullSpaceRelaxed
	: public MatrixWithOpFactorized
	, public MatrixConvertToSparseFortranCompatible
{
public:
	
	///
	typedef SparseSolverPack::DirectSparseFortranCompatibleSolver
		DirectSparseFortranCompatibleSolver;

	///
	class NotInitializedException : public std::logic_error
	{public: NotInitializedException (const std::string& what_arg) : std::logic_error(what_arg) {}};

	///
	class SingularMatrixException : public std::logic_error
	{public: SingularMatrixException (const std::string& what_arg) : std::logic_error(what_arg) {}};

	///
	class InvalidMatrixType : public std::logic_error
	{public: InvalidMatrixType (const std::string& what_arg) : std::logic_error(what_arg) {}};

	///
	enum ERunTests { RUN_TESTS, NO_TESTS };

	///
	enum EPrintMoreOrLess { PRINT_MORE, PRINT_LESS };

	/// <<std comp>> members for the direct sparse linear solver
	STANDARD_COMPOSITION_MEMBERS( DirectSparseFortranCompatibleSolver, direct_solver )

	///
	MatrixKKTFullSpaceRelaxed( const direct_solver_ptr_t& direct_solver = 0 );

	/** @name Initialize the relaxed or unrelaxed KKT matrix.
	  *
	  * These operations will factorize the matrix K.  If the matrix K
	  * is not full rank then a SingularMatrixException exception will
	  * be thrown.  The objects G and A must support the
	  * MatrixConvertToSparseFortranCompatible (MCTSFC) interface or
	  * the exception InvalidMatrixType will be thrown.
	  *
	  * Some of the common arguments that these initialization methods share
	  * are:
	  *
	  * \begin{itemize}
	  *	\item G [I] Hessian matrix ( must support MESFCE interface ).
	  *	\item A [I] Gradient of constraints matrix ( must support MESFCE interface ).
	  *	\item out [O] Output stream to print to.  This stream may be used for
	  *		output after initialization also so make sure that it remains valid
	  *		as long as this matrix object is is use.  For no output set out=NULL.
	  *	\item run_test [I] If set the true then many (expensive) tests will be
	  *		preformed to ensure that everything is working properly. 
	  *	\item print_more [I] If set the true then a lot more output may be produced
	  *		expecially if some error occurs. 
	  * \end{itemize}
	  *
	  * Important: It is vital that the definitions of G and A do not change
	  * externally while this object is being used.  To do so may invalidate
	  * the behavior of this object (especially the MatrixWithOp functions).
	  *
	  * This class will try to reuse the factorization structure from the
	  * last call to initialze(...) or initialize_relaxed(...) when possible.
	  * Namely if G and A have the same dimensions and same number of nonzeros
	  * of the matrices previously factorized, it will be assumed that the
	  * structure will be the same.  If this is not the case then the
	  * client should call release_memory(...) to wipe the slate clean and
	  * start over before calling initialize...(...) again.
	  */
	//@{

	///
	/** Initialize the nonrelaxed matrix.
	  *
	  */
	void initialize( const MatrixWithOp& G, const MatrixWithOp& A
		, std::ostream* out = 0, EPrintMoreOrLess print_what = PRINT_LESS
		, ERunTests test_what = NO_TESTS );

	///
	/** Initialize the relaxed matrix.
	  *
	  * If the unrelaxed QP is well scaled (near 1.0) then a reasonable
	  * value for bigM = M might be 1e+10 however this is problem specific.
	  */
	void initialize_relaxed( const MatrixWithOp& G, const MatrixWithOp& A
		, const VectorSlice& c, value_type bigM = 1e+10
		, std::ostream* out = 0, EPrintMoreOrLess print_what = PRINT_LESS
		, ERunTests test_what = NO_TESTS );

	///
	/** Set the matrix to uninitialized.
	  *
	  * The purpose of this method is for the client to specifically state that
	  * it is done using this object for now.  This is to avoid problems where
	  * the definitions of G and A might change and then another client unknowingly
	  * trys to use this object.
	  *
	  * Note that this does not erase storage of the factorization structure
	  * for example.
	  */
	void set_uninitialized();

	///
	/** Clear all allocated storage.
	  *
	  * The client should call this routine if he wants the new KKT matrix
	  * to be reanalyze and factorized the next time initialize...(...) is
	  * called.
	  */
	void release_memory();

	//@}

	// /////////////////////////////////////////////////////
	// Overridden from Matrix

	///
	size_type rows() const;

	///
	size_type cols() const;

	// /////////////////////////////////////////////////////////
	// Overridden from MatrixWithOp

	///
	std::ostream& output(std::ostream& out) const;

	///
	MatrixWithOp& operator=(const MatrixWithOp& m);

	/// (2) vs_lhs = alpha * op(M_rhs1) * vs_rhs2 + beta * vs_lhs (BLAS xGEMV)
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;

	// ////////////////////////////////////////////////////////////
	// Overridden from MatrixFactorized

	/// (1) v_lhs	= inv(op(M_rhs1)) * vs_rhs2
	void V_InvMtV( VectorSlice* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2) const;

	// ////////////////////////////////////////////////////////////
	// Overridden from MatrixConvertToSparseFortranCompatible

	///
	FortranTypes::f_int num_nonzeros( EExtractRegion extract_region ) const;

	///
	void coor_extract_nonzeros(
		  EExtractRegion extract_region
		, const FortranTypes::f_int len_Aval
			, FortranTypes::f_dbl_prec Aval[]
		, const FortranTypes::f_int len_Aij
			, FortranTypes::f_int Arow[]
			, FortranTypes::f_int Acol[]
			, const FortranTypes::f_int row_offset
			, const FortranTypes::f_int col_offset
		 ) const;

private:

	// //////////////////////////////
	// Private data members

	bool				initialized_;
	size_type			n_;		// Number of rows and columns in G and number of rows in A.
	size_type			m_;		// Number of columns in A
	bool				use_relaxation_;
	value_type			bigM_;
	EPrintMoreOrLess	print_what_;
	ERunTests			test_what_;
	std::ostream		*out_;
	const MatrixWithOp	*G_;
	const MatrixConvertToSparseFortranCompatible
						*convG_;
	size_type			G_nz_;	// Remember the number of nonzeros of G
	const MatrixWithOp	*A_;
	const MatrixConvertToSparseFortranCompatible
						*convA_;
	size_type			A_nz_;	// Remember the number of nonzeros of A

	// //////////////////////////////
	// Private member functions

	///
	void assert_matrices_set() const;

	///
	void assert_initialized() const;

	///
	/** Validate the types and sizes of G and A, set the member pointers G_ and A_
	  * and return the conversion interfaces convG and convA.
	  */
	void validate_and_set_matrices( const MatrixWithOp& G, const MatrixWithOp& A );

};	// end class MatrixKKTFullSpaceRelaxed

}	// end namespace ConstrainedOptimizationPack

#endif	// MATRIX_KKT_FULL_SPACE_RELAXED_H
