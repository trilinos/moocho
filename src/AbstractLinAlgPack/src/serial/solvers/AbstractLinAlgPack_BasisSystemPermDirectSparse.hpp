// ////////////////////////////////////////////////////////////
// BasisSystemPermDirectSparse.hpp
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

#ifndef ABSTRACT_LIN_ALG_PACK_BASIS_PERM_DIRECT_SPARSE_SYSTEM_H
#define ABSTRACT_LIN_ALG_PACK_BASIS_PERM_DIRECT_SPARSE_SYSTEM_H

#include "DirectSparseSolver.hpp"
#include "AbstractLinAlgPack/src/BasisSystemPerm.hpp"
#include "DenseLinAlgPack/src/IVector.hpp"

namespace SparseSolverPack {

///
/** Permutatble basis system subclass that uses a direct sparse solver.
 *
 * This current implementation only allows undecomposed general inequality
 * constraints.
 *
 * ToDo: Finish documentation!
 */
class BasisSystemPermDirectSparse
	: public AbstractLinAlgPack::BasisSystemPerm
{
public:

	/** @name Public types */
	//@{

	///
	typedef MemMngPack::ref_count_ptr<DirectSparseSolver>   direct_solver_ptr_t;

	//@}

	/** @name Constructors / initializers */
	//@{

	/// Calls <tt>this->initialize()</tt>
	BasisSystemPermDirectSparse(
		const direct_solver_ptr_t&   direct_solver = MemMngPack::null
		);
	
	/// Initialize given a direct sparse solver object.
	void initialize(
		const direct_solver_ptr_t&   direct_solver
		);

	//@}

	/** @name Overridden from BasisSystem */
	//@{

	///
	const mat_nonsing_fcty_ptr_t factory_C() const;
	///
	const mat_fcty_ptr_t factory_D() const;
	///
	const mat_fcty_ptr_t factory_GcUP() const;
	///
	const mat_fcty_ptr_t factory_GhUP() const;
	///
	Range1D var_dep() const;
	///
	Range1D var_indep() const;
	///
	Range1D equ_decomp() const;
	///
	Range1D equ_undecomp() const;
	///
	void update_basis(
		const MatrixOp*         Gc
		,const MatrixOp*        Gh
		,MatrixOpNonsing*   C
		,MatrixOp*              D
		,MatrixOp*              GcUP
		,MatrixOp*              GhUP
		,EMatRelations              mat_rel
		,std::ostream               *out
		) const;

	//@}

	/** @name Overridded from BasisSystemPerm */
	//@{

	///
	const perm_fcty_ptr_t   factory_P_var() const;
	///
	const perm_fcty_ptr_t   factory_P_equ() const;
	///
	const perm_fcty_ptr_t   factory_P_inequ() const;
	///
	void set_basis(
		const Permutation          &P_var
		,const Range1D             &var_dep
		,const Permutation         *P_equ
		,const Range1D             *equ_decomp
		,const Permutation         *P_inequ
		,const Range1D             *inequ_decomp
		,const MatrixOp        *Gc
		,const MatrixOp        *Gh
		,MatrixOpNonsing   *C
		,MatrixOp              *D
		,MatrixOp              *GcUP
		,MatrixOp              *GhUP
		,EMatRelations              mat_rel
		,std::ostream               *out
		);
	///
	void select_basis(
		const Vector          *nu
		,const Vector         *lambdaI
		,MatrixOp               *Gc
		,MatrixOp               *Gh
		,Permutation                *P_var
		,Range1D                    *var_dep
		,Permutation                *P_equ
		,Range1D                    *equ_decomp
		,Permutation                *P_inequ
		,Range1D                    *inequ_decomp
		,MatrixOpNonsing    *C
		,MatrixOp               *D
		,MatrixOp               *GcUP
		,MatrixOp               *GhUP
		,EMatRelations              mat_rel
		,std::ostream               *out
		);
	
	//@}

private:

	// ///////////////////////////////
	// Private data members

	direct_solver_ptr_t   direct_solver_;
	size_type             n_;
	size_type	          m_;
	size_type             mI_;
	size_type             r_;
	size_type             Gc_nz_;
	Range1D               init_var_rng_;
	IVector               init_var_inv_perm_;  // If init_var_rng is full range then this is ignored
    Range1D               init_equ_rng_;
	IVector               init_equ_inv_perm_;  // If init_equ_rng is full range then this is ignored
	Range1D               var_dep_;       // used by factor()
	Range1D               var_indep_;     // used by factor()
    Range1D               equ_decomp_;    // used by factor()
    Range1D               equ_undecomp_;  // used by factor()

	// ///////////////////////////////
	// Private member functions

	///
	MemMngPack::ref_count_ptr<DirectSparseSolver::BasisMatrix>
	get_basis_matrix( MatrixOpNonsingAggr &C_aggr ) const;

	///
	void set_A_mctse(
		size_type                    n
		,size_type                   m
		,const MatrixPermAggr        &Gc_pa
		,MatrixConvertToSparseEncap  *A_mctse
		) const;

	///
	void update_basis_and_auxiliary_matrices(
		const MatrixOp& Gc
		,const MemMngPack::ref_count_ptr<DirectSparseSolver::BasisMatrix>& C_bm
		,MatrixOpNonsingAggr *C_aggr
		,MatrixOp* D, MatrixOp* GcUP, MatrixOp* GhUp
		) const;

	///
	void do_some_basis_stuff(
		const MatrixOp& Gc, const MatrixOp* Gh
		,const Range1D& var_dep, const Range1D& equ_decomp
		,const MemMngPack::ref_count_ptr<DirectSparseSolver::BasisMatrix>& C_bm
		,MatrixOpNonsingAggr *C_aggr
		,MatrixOp* D, MatrixOp* GcUP, MatrixOp* GhUP
		);


}; // end class BasisSystemPermDirectSparse

}  // end namespace SparseSolverPack

#endif // ABSTRACT_LIN_ALG_PACK_BASIS_PERM_DIRECT_SPARSE_SYSTEM_H
