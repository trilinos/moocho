// ////////////////////////////////////////////////////////////
// BasisSystemPermDirectSparse.cpp
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

#include "SparseSolverPack/src/BasisSystemPermDirectSparse.hpp"
#include "SparseLinAlgPack/src/PermutationSerial.hpp"
#include "SparseLinAlgPack/src/MatrixConvertToSparseEncap.hpp"
#include "SparseLinAlgPack/src/MatrixExtractSparseElements.hpp"
#include "SparseLinAlgPack/src/MatrixSymPosDefCholFactor.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/MatrixOp.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/MatrixOpNonsingAggr.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/MatrixPermAggr.hpp"
#include "SparseLinAlgPack/src/MultiVectorMutableDense.hpp"
#include "AbstractFactoryStd.hpp"
#include "ThrowException.hpp"
#include "dynamic_cast_verbose.hpp"

namespace SparseSolverPack {

BasisSystemPermDirectSparse::BasisSystemPermDirectSparse(
	const direct_solver_ptr_t&   direct_solver
	)
	:BasisSystemPerm(
		MemMngPack::rcp(
			new MemMngPack::AbstractFactoryStd<MatrixSymOp,MatrixSymPosDefCholFactor>())               // D'*D
		,MemMngPack::rcp(
			new MemMngPack::AbstractFactoryStd<MatrixSymOpNonsing,MatrixSymPosDefCholFactor>())    // S
		)
{
	this->initialize(direct_solver);
}
	
void BasisSystemPermDirectSparse::initialize(
	const direct_solver_ptr_t&   direct_solver
	)
{
	direct_solver_ = direct_solver;
	n_ = m_ = mI_ = r_ = Gc_nz_ = 0;
	init_var_inv_perm_.resize(0);
	init_equ_inv_perm_.resize(0);
	init_var_rng_    = Range1D::Invalid;
    init_equ_rng_    = Range1D::Invalid;
	var_dep_         = Range1D::Invalid;
	var_indep_       = Range1D::Invalid;
    equ_decomp_      = Range1D::Invalid;
    equ_undecomp_    = Range1D::Invalid;
}

// Overridden from BasisSystem

const BasisSystem::mat_nonsing_fcty_ptr_t
BasisSystemPermDirectSparse::factory_C() const
{
	return MemMngPack::rcp(
		new MemMngPack::AbstractFactoryStd<MatrixOpNonsing,MatrixOpNonsingAggr>()
		);
}

const BasisSystem::mat_fcty_ptr_t
BasisSystemPermDirectSparse::factory_D() const
{
	return MemMngPack::rcp(
		new MemMngPack::AbstractFactoryStd<MatrixOp,MultiVectorMutableDense>()
		);
}

const BasisSystem::mat_fcty_ptr_t
BasisSystemPermDirectSparse::factory_GcUP() const
{
	return MemMngPack::rcp(
		new MemMngPack::AbstractFactoryStd<MatrixOp,MultiVectorMutableDense>()
		);
}

const BasisSystem::mat_fcty_ptr_t
BasisSystemPermDirectSparse::factory_GhUP() const
{
	return MemMngPack::rcp(
		new MemMngPack::AbstractFactoryStd<MatrixOp,MultiVectorMutableDense>()
		);
}

Range1D BasisSystemPermDirectSparse::var_dep() const
{
	return var_dep_;
}

Range1D BasisSystemPermDirectSparse::var_indep() const
{
	return var_indep_;
}

Range1D BasisSystemPermDirectSparse::equ_decomp() const
{
	return equ_decomp_;
}

Range1D BasisSystemPermDirectSparse::equ_undecomp() const
{
	return equ_undecomp_;
}

void BasisSystemPermDirectSparse::update_basis(
	const MatrixOp*         Gc
	,const MatrixOp*        Gh
	,MatrixOpNonsing*   C
	,MatrixOp*              D
	,MatrixOp*              GcUP
	,MatrixOp*              GhUP
	,EMatRelations              mat_rel
	,std::ostream               *out
	) const
{
	namespace mmp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	if(out)
		*out << "\nUsing a direct sparse solver to update basis ...\n";
#ifdef _DEBUG
	// Validate input
	THROW_EXCEPTION(
		Gc == NULL, std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error, "
		"Must have equality constriants in this current implementation! " );
#endif
	const size_type
		n  = Gc->rows(),
		m  = Gc->cols(),
		mI = Gh ? Gh->cols() : 0;
#ifdef _DEBUG
	const size_type Gc_rows = n, Gc_cols = m, Gc_nz = Gc->nz()
		,Gh_rows = Gh ? Gh->rows() : 0, Gh_cols = mI;
	THROW_EXCEPTION(
		Gc_rows != n_ || Gc_cols != m_ || Gc_nz != Gc_nz_, std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error, "
		"This matrix object is not compatible with last call to set_basis() or select_basis()!" );
	THROW_EXCEPTION(
		Gh != NULL && (Gh_rows != Gc_rows), std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error, "
		"Gc and Gh are not compatible!" )
	THROW_EXCEPTION(
		C == NULL, std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error!" );
#endif
	// Get the aggregate matrix object for Gc
	const MatrixPermAggr	
		&Gc_pa = dyn_cast<const MatrixPermAggr>(*Gc);
	// Get the basis matrix object from the aggregate or allocate one
	MatrixOpNonsingAggr
		&C_aggr = dyn_cast<MatrixOpNonsingAggr>(*C);
	mmp::ref_count_ptr<DirectSparseSolver::BasisMatrix>
		C_bm = get_basis_matrix(C_aggr);
	// Setup the encapulated convert-to-sparse matrix object
	MatrixConvertToSparseEncap A_mctse;
	set_A_mctse( n, m, Gc_pa, &A_mctse );
	// Refactor this basis (it had better be full rank)!
	try {
		direct_solver_->factor(
			A_mctse
			,C_bm.get()
			,mmp::null // Same factorization structure as before
			,out
			);
	}
	catch(const DirectSparseSolver::FactorizationFailure& excpt) {
		if(out)
			*out << "\nCurrent basis is singular : " << excpt.what() << std::endl
				 << "Throwing SingularBasis exception to client ...\n";
		THROW_EXCEPTION(
			true, SingularBasis
			,"BasisSystemPermDirectSparse::update_basis(...) : Error, the current basis "
			"is singular : " << excpt.what() );
	}
	// Update the aggregate basis matrix and compute the auxiliary projected matrices
	update_basis_and_auxiliary_matrices( *Gc, C_bm, &C_aggr, D, GcUP, GhUP );
}

// Overridded from BasisSystemPerm

const AbstractLinAlgPack::BasisSystemPerm::perm_fcty_ptr_t
BasisSystemPermDirectSparse::factory_P_var() const
{
	assert(0); // ToDo: Implement using PermutationSerial
	return MemMngPack::null;
}

const AbstractLinAlgPack::BasisSystemPerm::perm_fcty_ptr_t
BasisSystemPermDirectSparse::factory_P_equ() const
{
	assert(0); // ToDo: Implement using PermutationSerial
	return MemMngPack::null;
}

const AbstractLinAlgPack::BasisSystemPerm::perm_fcty_ptr_t
BasisSystemPermDirectSparse::factory_P_inequ() const
{
	assert(0); // ToDo: Implement using PermutationSerial
	return MemMngPack::null;
}

void BasisSystemPermDirectSparse::set_basis(
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
	)
{
	namespace mmp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	if(out)
		*out << "\nUsing a direct sparse solver to set a new basis ...\n";
#ifdef _DEBUG
	// Validate input
	THROW_EXCEPTION(
		Gc == NULL, std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error, "
		"Must have equality constriants in this current implementation! " );
#endif
	const size_type
		n  = Gc->rows(),
		m  = Gc->cols(),
		mI = Gh ? Gh->cols() : 0;
#ifdef _DEBUG
	const size_type Gc_rows = n, Gc_cols = m, Gc_nz = Gc->nz()
		,Gh_rows = Gh ? Gh->rows() : 0, Gh_cols = mI;
	THROW_EXCEPTION(
		P_equ == NULL || equ_decomp == NULL, std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error!" );
	THROW_EXCEPTION(
		Gh != NULL && (Gh_rows != Gc_rows), std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error, "
		"Gc and Gh are not compatible!" )
	THROW_EXCEPTION(
		Gh != NULL && (P_inequ != NULL || inequ_decomp != NULL)
		,std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error, "
		"Can not handle decomposed inequalities yet!" );
	THROW_EXCEPTION(
		C == NULL, std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error!" );
	if(Gh) THROW_EXCEPTION(
		var_dep.size() != equ_decomp->size(), std::invalid_argument
		,"BasisSystemPermDirectSparse::set_basis(...) : Error!" );
#endif
	// Get the aggreate matrix object for Gc
	const MatrixPermAggr	
		&Gc_pa = dyn_cast<const MatrixPermAggr>(*Gc);
	// Get the basis matrix object from the aggregate or allocate one
	MatrixOpNonsingAggr
		&C_aggr = dyn_cast<MatrixOpNonsingAggr>(*C);
	mmp::ref_count_ptr<DirectSparseSolver::BasisMatrix>
		C_bm = get_basis_matrix(C_aggr);
	// Get at the concreate permutation vectors
	const PermutationSerial
		&P_var_s = dyn_cast<const PermutationSerial>(P_var),
		&P_equ_s = dyn_cast<const PermutationSerial>(*P_equ);
	// Setup the encapulated convert-to-sparse matrix object
	init_var_inv_perm_  = *P_var_s.inv_perm();
	init_var_rng_       = var_dep;
	init_equ_inv_perm_  = *P_equ_s.inv_perm();
    init_equ_rng_       = *equ_decomp;
	MatrixConvertToSparseEncap A_mctse;
	set_A_mctse( n, m, Gc_pa, &A_mctse );
	// Analyze and factor this basis (it had better be full rank)!
	IVector row_perm_ds, col_perm_ds; // Must store these even though we don't want them!
	size_type rank = 0;
	direct_solver_->analyze_and_factor(
		A_mctse
		,&row_perm_ds
		,&col_perm_ds
		,&rank
		,C_bm.get()
		,out
		);
	if( rank < var_dep.size() ) {
		assert(0); // ToDo: Throw an exception with a good error message!
	}
	// Update the rest of the basis stuff
	do_some_basis_stuff(*Gc,Gh,var_dep,*equ_decomp,C_bm,&C_aggr,D,GcUP,GhUP);
}

void BasisSystemPermDirectSparse::select_basis(
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
	)
{
	namespace mmp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	if(out)
		*out << "\nUsing a direct sparse solver to select a new basis ...\n";
#ifdef _DEBUG
	// Validate input
	const char msg_err_head[] = "BasisSystemPermDirectSparse::set_basis(...) : Error!";
	THROW_EXCEPTION(
		Gc == NULL, std::invalid_argument
		,msg_err_head << " Must have equality constriants in this current implementation! " );
#endif
	const size_type
		n  = Gc->rows(),
		m  = Gc->cols(),
		mI = Gh ? Gh->cols() : 0;
#ifdef _DEBUG
	// Validate input
	const size_type Gc_rows = Gc->rows(), Gc_cols = Gc->cols(), Gc_nz = Gc->nz()
		,Gh_rows = Gh ? Gh->rows() : 0, Gh_cols = Gh ? Gh->cols() : 0;
	THROW_EXCEPTION(
		P_var == NULL || var_dep == NULL, std::invalid_argument, msg_err_head );
	THROW_EXCEPTION(
		P_equ == NULL || equ_decomp == NULL, std::invalid_argument, msg_err_head );
	THROW_EXCEPTION(
		Gh != NULL && (Gh_rows != Gc_rows), std::invalid_argument
		,msg_err_head << " Gc and Gh are not compatible!" )
	THROW_EXCEPTION(
		Gh != NULL && (P_inequ != NULL || inequ_decomp != NULL), std::invalid_argument
		,msg_err_head << "Can not handle decomposed inequalities yet!" );
	THROW_EXCEPTION(
		C == NULL, std::invalid_argument, msg_err_head );
#endif
	// Get the aggreate matrix object for Gc
	MatrixPermAggr	
		&Gc_pa = dyn_cast<MatrixPermAggr>(*Gc);
	// Get the basis matrix object from the aggregate or allocate one
	MatrixOpNonsingAggr
		&C_aggr = dyn_cast<MatrixOpNonsingAggr>(*C);
	mmp::ref_count_ptr<DirectSparseSolver::BasisMatrix>
		C_bm = get_basis_matrix(C_aggr);
	// Setup the encapulated convert-to-sparse matrix object
	// ToDo: Use nu to exclude variables that are at a bound!
	init_var_rng_       = Range1D(1,n);
	init_var_inv_perm_.resize(0); // Not used since above is full range
	init_equ_rng_       = Range1D(1,m);
	init_equ_inv_perm_.resize(0); // Not used since above is full range
	MatrixConvertToSparseEncap A_mctse;
	set_A_mctse( n, m, Gc_pa, &A_mctse );
	// Analyze and factor this basis (it had better be full rank)!
	mmp::ref_count_ptr<IVector>
		var_perm_ds = mmp::rcp(new IVector),
		equ_perm_ds = mmp::rcp(new IVector);
	size_type rank = 0;
	direct_solver_->analyze_and_factor(
		A_mctse
		,equ_perm_ds.get()
		,var_perm_ds.get()
		,&rank
		,C_bm.get()
		,out
		);
	if( rank == 0 ) {
		assert( rank == 0 ); // ToDo: Throw exception with good error message!
	}
	// Return the selected basis
	// ToDo: Use var_perm_ds and equ_perm_ds together with nu to
	// get the overall permuations for all of the variables.
	PermutationSerial
		&P_var_s = dyn_cast<PermutationSerial>(*P_var),
		&P_equ_s = dyn_cast<PermutationSerial>(*P_equ);
	// Create the overall permutations to set to the permutation matrices!
	*var_dep = Range1D(1,rank);
	P_var_s.initialize( var_perm_ds, mmp::null ); 
	*equ_decomp = Range1D(1,rank);
	P_equ_s.initialize( equ_perm_ds, mmp::null );
	// Setup Gc_aggr with Gc_perm
	const int          num_row_part = 2;
	const int          num_col_part = 2;
	const index_type   row_part[3]  = { 1, rank, n+1 };
	const index_type   col_part[3]  = { 1, rank, m+1 };
	Gc_pa.initialize(
		Gc_pa.mat_orig()
		,mmp::rcp(new PermutationSerial(var_perm_ds,mmp::null)) // var_perm_ds reuse is okay!
		,mmp::rcp(new PermutationSerial(equ_perm_ds,mmp::null)) // equ_perm_ds resue is okay!
		,Gc_pa.mat_orig()->perm_view(
			P_var,row_part,num_row_part
			,P_equ,col_part,num_col_part
			)
		);
	// Update the rest of the basis stuff
	do_some_basis_stuff(*Gc,Gh,*var_dep,*equ_decomp,C_bm,&C_aggr,D,GcUP,GhUP);
}

// private

MemMngPack::ref_count_ptr<DirectSparseSolver::BasisMatrix>
BasisSystemPermDirectSparse::get_basis_matrix( MatrixOpNonsingAggr &C_aggr ) const
{
	namespace mmp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	mmp::ref_count_ptr<DirectSparseSolver::BasisMatrix> C_bm;
	if( C_aggr.mns().get() ) {
		C_bm = mmp::rcp_dynamic_cast<DirectSparseSolver::BasisMatrix>(
			mmp::rcp_const_cast<MatrixNonsing>(C_aggr.mns() ) );
		if(C_bm.get() == NULL)
			dyn_cast<const DirectSparseSolver::BasisMatrix>(*C_aggr.mns()); // Throws exception!
	}
	else {
		C_bm = direct_solver_->basis_matrix_factory()->create();
	}
	return C_bm;
}

void BasisSystemPermDirectSparse::set_A_mctse(
	size_type                    n
	,size_type                   m
	,const MatrixPermAggr        &Gc_pa
	,MatrixConvertToSparseEncap  *A_mctse
	) const
{
	namespace mmp = MemMngPack;
	A_mctse->initialize(
		mmp::rcp_dynamic_cast<const MatrixExtractSparseElements>(Gc_pa.mat_orig())
		,mmp::rcp( init_var_rng_.size()    < n ? &init_var_inv_perm_ : NULL, false )
		,init_var_rng_
		,mmp::rcp( init_equ_rng_.size() < m ? &init_equ_inv_perm_ : NULL, false )
		,init_equ_rng_
		,BLAS_Cpp::trans
		);
}

void BasisSystemPermDirectSparse::update_basis_and_auxiliary_matrices(
	const MatrixOp& Gc
	,const MemMngPack::ref_count_ptr<DirectSparseSolver::BasisMatrix>& C_bm
	,MatrixOpNonsingAggr *C_aggr
	,MatrixOp* D, MatrixOp* GcUP, MatrixOp* GhUP
	) const
{
	namespace mmp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	// Initialize the aggregate basis matrix object.
	C_aggr->initialize(
		Gc.sub_view(var_dep_,equ_decomp_)
		,BLAS_Cpp::trans
		,C_bm
		,BLAS_Cpp::no_trans
		);
	// Compute the auxiliary projected matrices
	// Get the concreate type of the direct sensitivity matrix (if one was passed in)
	if( D ) {
		MultiVectorMutableDense *D_mvd = &dyn_cast<MultiVectorMutableDense>(*D);
		assert( D ); // ToDo: Throw exception!
		// D = -inv(C) * N
		D_mvd->initialize(var_dep_.size(),var_indep_.size());
		AbstractLinAlgPack::M_StInvMtM(
			D_mvd, -1.0, *C_bm, BLAS_Cpp::no_trans
			,*Gc.sub_view(var_indep_,equ_decomp_),BLAS_Cpp::trans // N = Gc(var_indep,equ_decomp)'
			);
	}
	if( GcUP ) {
		assert(0); // ToDo: Implement!
	}
	if( GhUP ) {
		assert(0); // ToDo: Implement!
	}
}

void BasisSystemPermDirectSparse::do_some_basis_stuff(
	const MatrixOp& Gc, const MatrixOp* Gh
	,const Range1D& var_dep, const Range1D& equ_decomp
	,const MemMngPack::ref_count_ptr<DirectSparseSolver::BasisMatrix>& C_bm
	,MatrixOpNonsingAggr *C_aggr
	,MatrixOp* D, MatrixOp* GcUP, MatrixOp* GhUP
	)
{
	const size_type
		n  = Gc.rows(),
		m  = Gc.cols(),
		mI = Gh ? Gh->cols() : 0;
	// Set the ranges
	var_dep_      = var_dep;
	var_indep_   = ( var_dep.size() == n
					  ? Range1D::Invalid
					  : ( var_dep.lbound() == 1
						  ? Range1D(var_dep.size()+1,n)
						  : Range1D(1,n-var_dep.size()) ) );
	equ_decomp_   = equ_decomp;
	equ_undecomp_ = ( equ_decomp.size() == m
					  ? Range1D::Invalid
					  : ( equ_decomp.lbound() == 1
						  ? Range1D(equ_decomp.size()+1,m)
						  : Range1D(1,m-equ_decomp.size()) ) );
	// Set the basis system dimensions
	n_     = n;
	m_     = m;
	mI_    = mI;
	r_     = var_dep.size();
	Gc_nz_ = Gc.nz();
	// Update the aggregate basis matrix and compute the auxiliary projected matrices
	update_basis_and_auxiliary_matrices( Gc, C_bm, C_aggr, D, GcUP, GhUP );
}
	
} // end namespace SparseSolverPack
