// /////////////////////////////////////////////////////////////
// DirectSparseSolverMA28.cpp
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

#ifdef SPARSE_SOLVER_PACK_USE_MA28

#include <assert.h>

#include <ostream>

#include "SparseSolverPack/include/DirectSparseSolverMA28.h"
#include "SparseSolverPack/include/MatrixScaling_Strategy.h"
#include "SparseLinAlgPack/include/VectorDenseEncap.h"
#include "LinAlgPack/include/PermVecMat.h"
#include "AbstractFactoryStd.h"
#include "ThrowException.h"
#include "WorkspacePack.h"
#include "dynamic_cast_verbose.h"

namespace SparseSolverPack {

//
// Implementation of DirectSparseSolver(Imp) interface using MA28.
//
// Here are some little bits of knowledge about MA28 that I need
// to record after many hours of hard work.
//
// * The 1979 paper in ACM TOMS (Vol. 5, No. 1, pages 27), seems
// to suggest that MA28 pivots by column for numerical stability
// but I am not sure about this.
//
// * When factoring a rectangular matrix, you must set 
// LBLOCK = .FALSE. or the row and column permutations
// extracted from IKEEP(:,2) and IKEEP(:,3) respectivly
// are meaningless.
//
// ToDo: Finish this discussion!
//

// ToDo:
// a) Add an option for printing the values of the common
//    block parameters to out or to a file.  This can
//    be used to get a feel for the performance of
//    ma28
// b) Add provisions for an external client to change
//    the control options of MA28.  Most of these are
//    stored as common block variables.

// //////////////////////////////////////////////////
// DirectSparseSolverMA28::FactorizationStructureMA28

DirectSparseSolverMA28::FactorizationStructureMA28::FactorizationStructureMA28()
	:m_(0),n_(0),max_n_(0),nz_(0),licn_(0),lirn_(0)
	 ,u_(0.1),scaling_(NO_SCALING)
	 ,cntr_var_changed_(false)
{}

// //////////////////////////////////////////////////
// DirectSparseSolverMA28::BasisMatrixMA28

// Overridden from BasisMatrixImp

MemMngPack::ref_count_ptr<DirectSparseSolverImp::BasisMatrixImp>
DirectSparseSolverMA28::BasisMatrixMA28::create_matrix() const
{
	return MemMngPack::rcp(new BasisMatrixMA28);
}

void DirectSparseSolverMA28::BasisMatrixMA28::V_InvMtV(
	VectorWithOpMutable* y, BLAS_Cpp::Transp M_trans, const VectorWithOp& x
	) const 
{
	using DynamicCastHelperPack::dyn_cast;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();
	size_type k;

	// Get concrete objects
	const FactorizationStructureMA28
		&fs = dyn_cast<const FactorizationStructureMA28>(*this->get_fact_struc());
	const FactorizationNonzerosMA28
		&fn = dyn_cast<const FactorizationNonzerosMA28>(*this->get_fact_nonzeros());

	// Validate input
#ifdef _DEBUG
	THROW_EXCEPTION(
		y == NULL, std::invalid_argument
		,"DirectSparseSolverMA28::BasisMatrixMA28::V_InvMtV(...) : Error! " );
	THROW_EXCEPTION(
		fs.cntr_var_changed_, std::logic_error
		,"DirectSparseSolverMA28::BasisMatrixMA28::V_InvMtV(...), Error, "
		"You can't change control variables between calls of \'analyze_and factorize\' and \'solve\'");	
#endif
	const size_type y_dim = y->dim(), x_dim = x.dim();
#ifdef _DEBUG
	THROW_EXCEPTION(
		fs.rank_ != y_dim, std::invalid_argument
		,"DirectSparseSolverMA28::BasisMatrixMA28::V_InvMtV(...) : Error! " );
	THROW_EXCEPTION(
		fs.rank_ != x_dim, std::invalid_argument
		,"DirectSparseSolverMA28::BasisMatrixMA28::V_InvMtV(...) : Error! " );
#endif

	VectorDenseMutableEncap  yd(*y);
	VectorDenseEncap         xd(x);

	// Allocate workspace memory
	wsp::Workspace<value_type>  xfull_s(wss,fs.max_n_,false);
	VectorSlice                 xfull(&xfull_s[0],xfull_s.size());
	wsp::Workspace<value_type>  ws(wss,fs.max_n_,false);
	VectorSlice                 w(&ws[0],ws.size());

	// Get a context for transpose or no transpose
	const IVector
		&row_perm = (M_trans == BLAS_Cpp::no_trans) ? fs.row_perm_ : fs.col_perm_,
		&col_perm = (M_trans == BLAS_Cpp::no_trans) ? fs.col_perm_ : fs.row_perm_;

	// Copy x into positions in full w
	// Here, the rhs vector is set with only those equations that
	// are part of the nonsingular set.  It is important that the
	// ordering be the same as the original ordering sent to
	// MA28AD().
	xfull = 0.0;
	for( k = 1; k <= x_dim; ++k ) 
		xfull(row_perm(k)) = xd()(k);
	
	// Scale the rhs
	if( fs.matrix_scaling_.get() )
		fs.matrix_scaling_->scale_rhs( M_trans, xfull.raw_ptr() );

	// Solve for the rhs
	FortranTypes::f_int mtype = ( (M_trans == BLAS_Cpp::no_trans) ? 1 : 0 );
	fs.ma28_.ma28cd(
		fs.max_n_, &fn.a_[0], fs.licn_, &fs.icn_[0], &fs.ikeep_[0]
		,xfull.raw_ptr(), w.raw_ptr(), mtype );

	// Scale the lhs
	if( fs.matrix_scaling_.get() )
		fs.matrix_scaling_->scale_rhs( M_trans, xfull.raw_ptr() );

	// Copy the solution into y
	// Here, the solution vector is set with only those variables that
	// are in the basis.  It is important that the
	// ordering be the same as the original ordering sent to
	// MA28AD().
	for( k = 1; k <= y_dim; ++k )
		yd()(k) = xfull(col_perm(k));
	
}

// //////////////////////////////////////////////////
// DirectSparseSolverMA28

// Constructors/initializers

DirectSparseSolverMA28::DirectSparseSolverMA28()
	: estimated_fillin_ratio_(10.0)
{}

// Overridden from DirectSparseSolver

const DirectSparseSolver::basis_matrix_factory_ptr_t
DirectSparseSolverMA28::basis_matrix_factory() const
{
	namespace mmp = MemMngPack;
	return mmp::rcp(new mmp::AbstractFactoryStd<BasisMatrix,BasisMatrixMA28>());
}

void DirectSparseSolverMA28::estimated_fillin_ratio(
	value_type estimated_fillin_ratio
	)
{
	estimated_fillin_ratio_ = estimated_fillin_ratio;
}

// Overridden from DirectSparseSolverImp

const MemMngPack::ref_count_ptr<DirectSparseSolver::FactorizationStructure>
DirectSparseSolverMA28::create_fact_struc() const
{
	return MemMngPack::rcp(new FactorizationStructureMA28);
}

const MemMngPack::ref_count_ptr<DirectSparseSolverImp::FactorizationNonzeros>
DirectSparseSolverMA28::create_fact_nonzeros() const
{
	return MemMngPack::rcp(new FactorizationNonzerosMA28);
}

void DirectSparseSolverMA28::imp_analyze_and_factor(
	const SparseLinAlgPack::MatrixConvertToSparse   &A
	,const BasisMatrix::fact_struc_ptr_t            &fact_struc
	,const BasisMatrixImp::fact_nonzeros_ptr_t      &fact_nonzeros
	,LinAlgPack::IVector                            *row_perm
	,LinAlgPack::IVector                            *col_perm
	,size_type                                      *rank
	,std::ostream                                   *out
	)
{
	using DynamicCastHelperPack::dyn_cast;
	typedef MatrixConvertToSparse MCTS;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	if(out)
		*out << "\nUsing MA28 to analyze and factor a new matrix ...\n";

	// Get the concrete factorization and nonzeros objects
	FactorizationStructureMA28
		&fs = dyn_cast<FactorizationStructureMA28>(*fact_struc);
	FactorizationNonzerosMA28
		&fn = dyn_cast<FactorizationNonzerosMA28>(*fact_nonzeros);

	// Set MA28 parameters
	fs.ma28_.nsrch(4);  // This is a well known trick to increase speed!
	fs.ma28_.lblock(0); // Do not permute to block triangular form (*** This is critical!)

	// Get the dimensions of things.
	const index_type
		m = A.rows(),
		n = A.cols(),
		nz = A.num_nonzeros( MCTS::EXTRACT_FULL_MATRIX ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM );

	// Validate input
	THROW_EXCEPTION(
		n <= 0 || m <= 0 || m > n, std::invalid_argument
		,"DirectSparseSolverMA28::imp_analyze_and_factor(...) : Error!" );

	// Sets control variable change flag to false.
	fs.cntr_var_changed_ = false;

	// Memorize the dimenstions for checks later
	fs.m_ = m; fs.n_ = n; fs.nz_ = nz;
	fs.max_n_ = std::_MAX(fs.m_,fs.n_);

	// By default set licn and ircn equal to fillin_ratio * nz.
	if( fs.licn_ < fs.nz_ ) fs.licn_ = estimated_fillin_ratio_ * fs.nz_;
	if( fs.lirn_ < fs.nz_ ) fs.lirn_ = estimated_fillin_ratio_ * fs.nz_;

	// Initialize matrix factorization storage and temporary storage
	fs.ivect_.resize(fs.nz_); // perminatly stores nz row indexes
	fs.jvect_.resize(fs.nz_); // perminatly stores nz column indexes
	fs.icn_.resize(fs.licn_); // first nz entries stores the column indexes
	fn.a_.resize(fs.licn_);
	fs.ikeep_.resize( fs.ma28_.lblock() ? 5*fs.max_n_ :  4*fs.max_n_ + 1 );
	wsp::Workspace<index_type>  irn_tmp_(wss,fs.lirn_), iw(wss,8*fs.max_n_);
	wsp::Workspace<value_type>  w(wss,fs.max_n_);

	// Fill in the matrix information
	A.coor_extract_nonzeros(
		MCTS::EXTRACT_FULL_MATRIX
		,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM
		,fs.nz_
		,&fn.a_[0]
		,fs.nz_
		,&fs.ivect_[0]
		,&fs.jvect_[0]
		);
	std::copy( &fs.ivect_[0], &fs.ivect_[0] + fs.nz_, &irn_tmp_[0] );
	std::copy( &fs.jvect_[0], &fs.jvect_[0] + fs.nz_, &fs.icn_[0] );

	// Scale the matrix
	if( fs.matrix_scaling_.get() )
		fs.matrix_scaling_->scale_matrix(
			fs.m_, fs.n_, fs.nz_, &fs.ivect_[0], &fs.jvect_[0], true
			,&fn.a_[0]
			);

	// Analyze and factor the matrix
	index_type iflag = 0;
	fs.ma28_.ma28ad(
		fs.max_n_, fs.nz_, &fn.a_[0], fs.licn_, &irn_tmp_[0], fs.lirn_, &fs.icn_[0], fs.u_
		,&fs.ikeep_[0], &iw[0], &w[0], &iflag
		);

	// Todo : deal with resizing for insufficient storage when needed.  Write a loop!

	if(iflag != 0 && out)
		*out << "\nMA28AD returned iflag = " << iflag << " != 0!\n";

	// Check for errors and throw exception if you have to.
	ThrowIFlagException(iflag);

	// Extract the basis matrix selection

	*rank = fs.ma28_.irank();

	row_perm->resize(fs.m_);
	if( *rank < fs.m_ ) {
		index_type
			*row_perm_ikeep = &fs.ikeep_[fs.max_n_],
			*row_perm_itr   = &(*row_perm)(1),
			*row_perm_last  = row_perm_itr + fs.m_;
		for(; row_perm_itr != row_perm_last;)
			*row_perm_itr++ = abs(*row_perm_ikeep++);
		// Sort partitions in assending order (required!)
		std::sort(&(*row_perm)[0]           , &(*row_perm)[0] + (*rank) );
		std::sort(&(*row_perm)[0] + (*rank)	, &(*row_perm)[0] + m       );
	}
	else {
		LinAlgPack::identity_perm( row_perm );
	}

	col_perm->resize(fs.n_);
	if( *rank < fs.n_ ) {
		index_type
			*col_perm_ikeep = &fs.ikeep_[2*fs.max_n_],
			*col_perm_itr   = &(*col_perm)(1),
			*col_perm_last  = col_perm_itr + fs.n_;
		for(; col_perm_itr != col_perm_last;)
			*col_perm_itr++ = abs(*col_perm_ikeep++);
		// Sort partitions in assending order (required!)
		std::sort(&(*col_perm)[0]           , &(*col_perm)[0] + (*rank) );
		std::sort(&(*col_perm)[0] + (*rank)	, &(*col_perm)[0] + n       );
	}
	else {
		LinAlgPack::identity_perm( col_perm );
	}

	// Set internal copy of basis selection
	fs.row_perm_ = *row_perm;
	fs.col_perm_ = *col_perm;
	fs.rank_     = *rank;

}

void DirectSparseSolverMA28::imp_factor(
	const SparseLinAlgPack::MatrixConvertToSparse   &A
	,const BasisMatrix::fact_struc_ptr_t            &fact_struc
	,const BasisMatrixImp::fact_nonzeros_ptr_t      &fact_nonzeros
	,std::ostream                                   *out
	)
{
	using DynamicCastHelperPack::dyn_cast;
	typedef MatrixConvertToSparse MCTS;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

	if(out)
		*out << "\nUsing MA28 to factor a new matrix ...\n";

	// Get the concrete factorization and nonzeros objects
	const FactorizationStructureMA28
		&fs = dyn_cast<FactorizationStructureMA28>(*fact_struc);
	FactorizationNonzerosMA28
		&fn = dyn_cast<FactorizationNonzerosMA28>(*fact_nonzeros);

	// Get the dimensions of things.
	const index_type
		m = A.rows(),
		n = A.cols(),
		nz = A.num_nonzeros( MCTS::EXTRACT_FULL_MATRIX ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM );

	// Validate input
#ifdef _DEBUG
	THROW_EXCEPTION(
		m != fs.m_ || n != fs.n_ || fs.nz_ != nz, std::invalid_argument
		,"DirectSparseSolverMA28::imp_factor(...) : Error, "
		"A is not compatible with matrix passed to imp_analyze_and_factor()!" );
#endif

	// Initialize matrix factorization storage and temporary storage
	if(fn.a_.size() < fs.licn_)  fn.a_.resize(fs.licn_);
	wsp::Workspace<index_type>   iw(wss,5*fs.max_n_);
	wsp::Workspace<value_type>   w(wss,fs.max_n_);

	// Fill in the matrix nonzeros (we already have the structure)
	A.coor_extract_nonzeros(
		MCTS::EXTRACT_FULL_MATRIX
		,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM
		,fs.nz_
		,&fn.a_[0]
		,0
		,NULL
		,NULL
		);

	// Scale the matrix
	if( fs.matrix_scaling_.get() )
		fs.matrix_scaling_->scale_matrix(
			fs.m_, fs.n_, fs.nz_, &fs.ivect_[0], &fs.jvect_[0], false
			,&fn.a_[0]
			);

	// Factor the matrix
	index_type iflag = 0;
	fs.ma28_.ma28bd(
		fs.max_n_, fs.nz_, &fn.a_[0], fs.licn_, &fs.ivect_[0], &fs.jvect_[0], &fs.icn_[0]
		,&fs.ikeep_[0], &iw[0], &w[0], &iflag
		);

	if(iflag != 0 && out)
		*out << "\nMA28BD returned iflag = " << iflag << " != 0!\n";

	// Check for errors and throw exception if you have to.
	ThrowIFlagException(iflag);

}

// private

void DirectSparseSolverMA28::ThrowIFlagException(index_type iflag)
{
	E_IFlag e_iflag = static_cast<E_IFlag>(iflag);
	const char msg_err_head[] = "DirectSparseSolverMA28::ThrowIFlagException(iflag) : Error";
	switch(e_iflag) {
		case SLOW_ITER_CONV :
			THROW_EXCEPTION(
				true, std::runtime_error
				,msg_err_head << ", Convergence to slow" );
		case MAXIT_REACHED :
			THROW_EXCEPTION(
				true, std::runtime_error
				,msg_err_head << ", Maximum iterations exceeded");
		case MA28BD_CALLED_WITH_DROPPED :
			THROW_EXCEPTION(
				true, std::logic_error
				,msg_err_head << ", ma28bd called with elements dropped in ma28ad");
		case DUPLICATE_ELEMENTS :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", Duplicate elements have been detected");
		case NEW_NONZERO_ELEMENT :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", A new non-zero element has be passed to ma28bd that was not ot ma28ad");
		case N_OUT_OF_RANGE :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", 1 <=max(n,m) <= 32767 has been violated");
		case NZ_LE_ZERO :
			THROW_EXCEPTION(
				true, std::logic_error
				,msg_err_head << ", nz <= 0 has been violated");
		case LICN_LE_NZ :
			THROW_EXCEPTION(
				true, std::logic_error
				,msg_err_head << ", licn <= nz has been violated");
		case LIRN_LE_NZ :
			THROW_EXCEPTION(
				true, std::logic_error
				,msg_err_head << ", lirn <= nz has been violated");
		case ERROR_DURRING_BLOCK_TRI :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", An error has occured durring block triangularization");
		case LICN_AND_LIRN_TOO_SMALL :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", licn and lirn are to small to hold matrix factorization");
		case LICN_TOO_SMALL :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", licn is to small to hold matrix factorization");
		case LICN_FAR_TOO_SMALL :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", licn is to far small to hold matrix factorization");
		case LIRN_TOO_SMALL :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", lirn is to small to hold matrix factorization");
		case NUMERICALLY_SINGULAR :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", matrix is numerically singular, see \'abort2\'");
		case STRUCTURALLY_SINGULAR :
			THROW_EXCEPTION(
				true, FactorizationFailure
				,msg_err_head << ", matrix is structurally singular, see \'abort1\'");
		default:
			return; // We don't throw exceptions for other values of iflag.
	}
}

}	// end namespace SparseSolverPack 

#endif // SPARSE_SOLVER_PACK_USE_MA28
