// ///////////////////////////////////////////////////////////////////
// SuperLUSolver.hpp
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

#ifdef SPARSE_SOLVER_PACK_USE_SUPERLU

#ifndef SSP_SUPERLU_SOLVER_H
#define SSP_SUPERLU_SOLVER_H

#include "ref_count_ptr.hpp"

namespace SuperLUPack {

///
/** Abstract interface to SuperLU.
 *
 * This class is abstract so that none of the SuperLU header files need
 * be presented to the client.  SuperLU uses some very bad software engineering
 * by injecting into the global namespace lots of inappropriate names.  All
 * SuperLU global names should be prefixed by SLU_ or SuperLU_ but they are not
 */
class SuperLUSolver {
public:

	/** @name Public Types */
	//@{

	///
	/** Abstract class for objects that represent the factorization
	  * structure of a particular class of matrices.
	  *
	  * This structure can be reused over and over again for factorizing matrices
	  * with the same structure but different nonzero elements.
	  */
	class FactorizationStructure {
	public:
		///
		virtual ~FactorizationStructure() {}
	};

	///
	/** Abstract class for objects that represent the factorization
	  * nonzeos of a particular matrix.
	  */
	class FactorizationNonzeros {
	public:
		///
		virtual ~FactorizationNonzeros() {}
	};

	//@}

	///
	virtual ~SuperLUSolver() {}

	/** @name Static members */
	//@{

	///
	static MemMngPack::ref_count_ptr<SuperLUSolver>                         create_solver();
	///
	static MemMngPack::ref_count_ptr<SuperLUSolver::FactorizationStructure> create_fact_struct();
	///
	static MemMngPack::ref_count_ptr<SuperLUSolver::FactorizationNonzeros>  create_fact_nonzeros();

	//@}

	/** @name Basis selection and factorization */
	//@{

	///
	/** Analyze and factor the matrix, finding a basis in the process
	 *
	 * param  m     [in] Number of rows in A
	 * param  n     [in] Number of columns in A
	 * param  nz    [in] Number of nonzeros in A
	 * param  a_val [in] Array (length nz) of the values of A
	 * param  a_row_i
	 *              [in] Array (length nz) of the rows indexes of A (zero-based
	 *              for SuperLU)
	 * param  a_col_ptr
	 *              [in] Array (length n+1) Column pointers for starts
	 *              of columns for A in a_val[] and a_col_ptr[]
	 * param  fact_struc
	 *              [out] Factorization structure for basis of A
	 * param  fact_nonzeros
	 *              [out] Factorization nonzeros for basis of A
	 * param  perm_r
	 *              [out] Array (length m) of row permutations for
	 *              basis selection of A (zero-based).
	 * param  perm_c
	 *              [out] Array (length n) of column permutations
	 *              for basis selection of A (zero-based).
	 * param  rank  [out] Rank of the basis of A selected.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>m >= n</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> Lots of others???
	 * </ul>
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void analyze_and_factor(
		int                         m
		,int                        n
		,int                        nz
		,const double               a_val[]
		,const int                  a_row_i[]
		,const int                  a_col_ptr[]
		,FactorizationStructure     *fact_struct
		,FactorizationNonzeros      *fact_nonzeros
		,int                        perm_r[]
		,int                        perm_c[]
		,int                        *rank
		) = 0;
	
	///
	/** Refactor the same basis.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void factor(
		int                             m
		,int                            n
		,int                            nz
		,const double                   a_val[]
		,const int                      a_row_i[]
		,const int                      a_col_ptr[]
		,const FactorizationStructure   &fact_struct
		,FactorizationNonzeros          *fact_nonzeros
		) = 0;

	//@}
	
	/** @name Solve for linear systems */
	//@{

	///
	/** Solve a set of linear systems.
	 *
	 * ToDo: Finish documentation!
	 */
	virtual void solve(
		const FactorizationStructure    &fact_struct
		,const FactorizationNonzeros    &fact_nonzeros
		,bool                           transp
		,int                            n
		,int                            nrhs
		,double                         rhs[]
		,int                            ldrhs
		) const = 0;

	//@}

}; // end class SuperLUSolver

} // end namespace SuperLUPack

#endif // SSP_SUPERLU_SOLVER_H

#endif // SPARSE_SOLVER_PACK_USE_SUPERLU
