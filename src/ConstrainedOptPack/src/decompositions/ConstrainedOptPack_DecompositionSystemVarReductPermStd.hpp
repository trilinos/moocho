// /////////////////////////////////////////////////////////////////////////////
// DecompositionSystemVarReductPermStd.h
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

#ifndef DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_STD_H
#define DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_STD_H

#include <stdexcept>

#include "DecompositionSystemVarReductPerm.h"
#include "AbstractLinAlgPack/src/VectorSpace.h"

namespace ConstrainedOptimizationPack {

///
/** Concreate subclass of \c DecompositionSystemVarReductPerm that uses an
 * aggregate \c DecompostionSystemVarReductImp object.
 *
 * ToDo: Finish documentation!
 */
class DecompositionSystemVarReductPermStd : public DecompositionSystemVarReductPerm {
public:

	/** @name Public types */
	//@{

	///
	typedef MemMngPack::ref_count_ptr<DecompositionSystemVarReductImp>    decomp_sys_imp_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<BasisSystemPerm>                    basis_sys_ptr_t;

	//@}

	/** @name Constructors / initializers */
	//@{

	/// Calls \c this->initialize().
	DecompositionSystemVarReductPermStd(
		const decomp_sys_imp_ptr_t&        decomp_sys_imp  = MemMngPack::null
		,const basis_sys_ptr_t&            basis_sys       = MemMngPack::null
		,bool                              basis_selected  = false
		,EExplicitImplicit                 D_imp           = MAT_IMP_AUTO
		,EExplicitImplicit                 Uz_imp          = MAT_IMP_AUTO
		,EExplicitImplicit                 Vz_imp          = MAT_IMP_AUTO
		);

	/// Initialize given decomposition system and basis system objects.
	void initialize(
		const decomp_sys_imp_ptr_t&        decomp_sys_imp
		,const basis_sys_ptr_t&            basis_sys
		,bool                              basis_selected  = false
		,EExplicitImplicit                 D_imp           = MAT_IMP_AUTO
		,EExplicitImplicit                 Uz_imp          = MAT_IMP_AUTO
		,EExplicitImplicit                 Vz_imp          = MAT_IMP_AUTO
		);

	//@}

	/** @name Access */
	//@{

	///
	const decomp_sys_imp_ptr_t& decomp_sys_imp() const;
	///
	const basis_sys_ptr_t& basis_sys() const;

	//@}

	/** @name Overridden from DecompositionSystem */
	//@{

	///
	size_type n() const;
	///
	size_type m() const;
	///
	size_type r() const;
	///
	Range1D equ_decomp() const;
	///
	Range1D equ_undecomp() const;
	///
	const VectorSpace::space_ptr_t space_range() const;
	///
	const VectorSpace::space_ptr_t space_null() const;
	///
	const mat_fcty_ptr_t factory_Z() const;
	///
	const mat_fcty_ptr_t factory_Y() const;
	///
	const mat_nonsing_fcty_ptr_t factory_R() const;
	///
	const mat_fcty_ptr_t factory_Uz() const;
	///
	const mat_fcty_ptr_t factory_Uy() const;
	///
	const mat_fcty_ptr_t factory_Vz() const;
	///
	const mat_fcty_ptr_t factory_Vy() const;
	///
	void update_decomp(
		std::ostream              *out
		,EOutputLevel             olevel
		,ERunTests                test_what
		,const MatrixWithOp       &Gc
		,const MatrixWithOp       *Gh
		,MatrixWithOp             *Z
		,MatrixWithOp             *Y
		,MatrixWithOpNonsingular  *R
		,MatrixWithOp             *Uz
		,MatrixWithOp             *Uy
		,MatrixWithOp             *Vz
		,MatrixWithOp             *Vy
		,EMatRelations            mat_rel
		) const;
	///
	void print_update_decomp(
		std::ostream& out, const std::string& leading_str ) const;

	//@}

	/** @name Overridden from DecompositionSystemVarReduct */
	//@{

	///
	Range1D var_indep() const;
	///
	Range1D var_dep() const;

	//@}

	/** @name Overridden from DecompositionSystemVarReductPerm */
	//@{

	///
	const perm_fcty_ptr_t   factory_P_var() const;
	///
	const perm_fcty_ptr_t   factory_P_equ() const;
	///
	bool has_basis() const;
	///
	void set_decomp(
		std::ostream              *out
		,EOutputLevel             olevel
		,ERunTests                test_what
		,const Permutation        &P_var
		,const Range1D            &var_dep
		,const Permutation        *P_equ
		,const Range1D            *equ_decomp
		,const MatrixWithOp       &Gc
		,const MatrixWithOp       *Gh
		,MatrixWithOp             *Z
		,MatrixWithOp             *Y
		,MatrixWithOpNonsingular  *R
		,MatrixWithOp             *Uz
		,MatrixWithOp             *Uy
		,MatrixWithOp             *Vz
		,MatrixWithOp             *Vy
		,EMatRelations            mat_rel
		);
	///
	void select_decomp(
		std::ostream              *out
		,EOutputLevel             olevel
		,ERunTests                test_what
		,const VectorWithOp       *nu
		,MatrixWithOp             *Gc
		,MatrixWithOp             *Gh
		,Permutation              *P_var
		,Range1D                  *var_dep
		,Permutation              *P_equ
		,Range1D                  *equ_decomp
		,MatrixWithOp             *Z
		,MatrixWithOp             *Y
		,MatrixWithOpNonsingular  *R
		,MatrixWithOp             *Uz
		,MatrixWithOp             *Uy
		,MatrixWithOp             *Vz
		,MatrixWithOp             *Vy
		,EMatRelations            mat_rel
		);

	//@}

private:

	// /////////////////////////
	// Private data members

	bool                        basis_selected_;  // True if a basis is currently selected
	decomp_sys_imp_ptr_t        decomp_sys_imp_;
	basis_sys_ptr_t             basis_sys_;

	// /////////////////////////
	// Private member functions

	///
	void assert_basis_selected() const;

	// Not defined and not to be called!
	DecompositionSystemVarReductPermStd();
	DecompositionSystemVarReductPermStd(const DecompositionSystemVarReductPermStd&);
	DecompositionSystemVarReductPermStd& operator=(const DecompositionSystemVarReductPermStd&);
	
};	// end class DecompositionSystemVarReductPermStd

// ///////////////////////////////////////
// Inline members

inline
const DecompositionSystemVarReductPermStd::decomp_sys_imp_ptr_t&
DecompositionSystemVarReductPermStd::decomp_sys_imp() const
{
	return decomp_sys_imp_;
}

inline
const DecompositionSystemVarReductPermStd::basis_sys_ptr_t&
DecompositionSystemVarReductPermStd::basis_sys() const
{
	return basis_sys_;
}

}	// end namespace ConstrainedOptimizationPack

#endif // DECOMPOSITION_SYSTEM_VAR_REDUCT_PERM_STD_H
