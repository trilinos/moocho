// /////////////////////////////////////////////////////////
// QPSchurInitKKTSystemHessianFixedFree.cpp

#include "ConstrainedOptimizationPack/include/QPSchurInitKKTSystemHessianFixedFree.h"
#include "ConstrainedOptimizationPack/include/initialize_Q_R_Q_X.h"
#include "SparseLinAlgPack/include/MatrixSymWithOp.h"
#include "SparseLinAlgPack/include/sparse_bounds.h"
#include "SparseLinAlgPack/include/GenPermMatrixSlice.h"
#include "SparseLinAlgPack/include/GenPermMatrixSliceOp.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "Misc/include/dynamic_cast_verbose.h"
#include "Misc/include/WorkspacePack.h"
#include "Misc/include/profile_hack.h"

namespace LinAlgOpPack {
    using SparseLinAlgPack::Vp_StMtV;
}

namespace ConstrainedOptimizationPack {

void QPSchurInitKKTSystemHessianFixedFree::initialize_kkt_system(
	const VectorSlice&    g
	,const MatrixWithOp&  G
	,value_type           etaL
	,const SpVectorSlice& dL
	,const SpVectorSlice& dU
	,const MatrixWithOp*  F
	,BLAS_Cpp::Transp     trans_F
	,const VectorSlice*   f
	,const VectorSlice&   d
	,const SpVectorSlice& nu
	,size_type*           n_R
	,i_x_free_t*          i_x_free
	,i_x_fixed_t*         i_x_fixed
	,bnd_fixed_t*         bnd_fixed
	,j_f_decomp_t*        j_f_decomp
	,Vector*              b_X
	,Ko_ptr_t*            Ko
	,Vector*              fo
	) const
{
	using DynamicCastHelperPack::dyn_cast;
	using LinAlgOpPack::V_mV;
	namespace rcp = ReferenceCountingPack;
	namespace wsp = WorkspacePack;
	wsp::WorkspaceStore* wss = WorkspacePack::default_workspace_store.get();

#ifdef PROFILE_HACK_ENABLED
	ProfileHackPack::ProfileTiming profile_timing( "QPSchurInitKKTSystemHessianFixedFree::initialize_kkt_system(...)" );
#endif

	// Validate type of and convert G
#ifdef _WINDOWS
	const MatrixSymWithOp&
		G_sym = dynamic_cast<const MatrixSymWithOp&>(G);
#else
	const MatrixSymWithOp&
		G_sym = dyn_cast<const MatrixSymWithOp>(G);
#endif

	const size_type nd = g.size();

	// Determine the number of initially fixed variables
	wsp::Workspace<EBounds> x_frfx(wss,nd);
	std::fill_n( &x_frfx[0], nd, FREE ); // make all free initially
	size_type
		num_init_fixed = 0;
	{
		const value_type inf_bnd = std::numeric_limits<value_type>::max();
		SparseLinAlgPack::sparse_bounds_itr
			dLU_itr(
				dL.begin(), dL.end(), dL.offset(),
				dU.begin(), dU.end(), dU.offset(), inf_bnd );
		SpVectorSlice::const_iterator
			nu_itr = nu.begin(),
			nu_end = nu.end();
		const SpVector::difference_type o = nu.offset();
		while( !dLU_itr.at_end() || nu_itr != nu_end ) {
			if( dLU_itr.at_end() ) { // Add the rest of the elements in nu
				for( ; nu_itr != nu_end; ++num_init_fixed, ++nu_itr )
					x_frfx[nu_itr->indice() + o - 1] = ( nu_itr->value() > 0.0 ? UPPER : LOWER );
			}
			else { // Be carefull to add fixed dL(i) == dU(i)
				// Add elements in nu up to the current dLU_itr.indice()
				for( ; nu_itr != nu_end && nu_itr->indice() + o < dLU_itr.indice(); ++num_init_fixed, ++nu_itr )
					x_frfx[nu_itr->indice() + o - 1] = ( nu_itr->value() > 0.0 ? UPPER : LOWER );
				if( dLU_itr.lbound() == dLU_itr.ubound() ) {
					// This is a permanently fixed variable!
					x_frfx[dLU_itr.indice() - 1] = EQUALITY;
					++num_init_fixed;
					// Don't add a duplicate entry in nu
					if( nu_itr != nu_end && nu_itr->indice() + o == dLU_itr.indice() )
						++nu_itr;
				}
				++dLU_itr;
			}
		}
	}
	assert( nd >= num_init_fixed );

	// n_R
	*n_R = nd - num_init_fixed;
	
	// Set up i_x_free[], i_x_fixed[], bnd_fixed[], and b_X
	i_x_free->resize(*n_R);
	i_x_fixed->resize(num_init_fixed+1);
	bnd_fixed->resize(num_init_fixed+1);
	b_X->resize(num_init_fixed+1);
	{
		const value_type inf_bnd = std::numeric_limits<value_type>::max();
		SparseLinAlgPack::sparse_bounds_itr
			dLU_itr(
				dL.begin(), dL.end(), dL.offset(),
				dU.begin(), dU.end(), dU.offset(), inf_bnd );
		size_type i_R = 0, i_X = 0;
		for( size_type i = 1; i <= nd; ++i ) {
			const EBounds
				bnd_i = x_frfx[i-1];
			if( bnd_i == FREE ) {
				(*i_x_free)[i_R] = i;
				++i_R;
			}
			else {
				(*i_x_fixed)[i_X] = i;
				(*bnd_fixed)[i_X] = bnd_i;
				assert( !dLU_itr.at_end() );    // find entry in b_X
				while( dLU_itr.indice() < i )
					++dLU_itr;
				assert( dLU_itr.indice() == i );
				value_type b_X_val = 0.0;
				switch( bnd_i ) {
					case EQUALITY:
					case LOWER:
						b_X_val = dLU_itr.lbound();
						break;
					case UPPER:
						b_X_val = dLU_itr.ubound();
						break;
					default:
						assert(0); // Local error only?
				}
				(*b_X)[i_X] = b_X_val;
				++i_X;
			}
		}
		(*i_x_fixed)[i_X] = nd+1;   // built-in relaxation variable
		(*bnd_fixed)[i_X] = LOWER;
		(*b_X)[i_X]       = etaL;
		++i_X;
	}
	
	// j_f_decomp[] = empty
	j_f_decomp->resize(0);

	// Initialize temporary Q_R and Q_X (not including extra relaxation variable)
	wsp::Workspace<size_type>
		Q_R_row_i(wss,*n_R),
		Q_R_col_j(wss,*n_R),
		Q_X_row_i(wss,num_init_fixed),
		Q_X_col_j(wss,num_init_fixed);
	GenPermMatrixSlice
		Q_R, Q_X;
	initialize_Q_R_Q_X(
		*n_R,num_init_fixed,&(*i_x_free)[0],&(*i_x_fixed)[0],false
		,&Q_R_row_i[0],&Q_R_col_j[0],&Q_R
		,&Q_X_row_i[0],&Q_X_col_j[0],&Q_X
		);

	//
	// Create and initialize object for Ko = G_RR = Q_R'*G*Q_R
	//

	// Compute the dense matrix G_RR
	GenMatrix G_RR_dense(*n_R,*n_R);
	sym_gms sym_G_RR_dense(G_RR_dense(),BLAS_Cpp::lower);
	SparseLinAlgPack::Mp_StPtMtP(
		&sym_G_RR_dense, 1.0, MatrixSymWithOp::DUMMY_ARG
		,G_sym, Q_R, BLAS_Cpp::no_trans, 0.0 );
	// Initialize a factorization object for this matrix
	typedef rcp::ref_count_ptr<MatrixSymPosDefCholFactor> G_RR_ptr_t;
	G_RR_ptr_t
		G_RR_ptr = new MatrixSymPosDefCholFactor();
	G_RR_ptr->initialize(sym_G_RR_dense);
	
	*Ko = rcp::rcp_implicit_cast<Ko_ptr_t::element_type>(G_RR_ptr); // Ko is initialized!

	// ToDo: (2001/07/05) We could be more carefull about how memory is initialized and reused
	// in the future but this implementation is just easier.

	// fo = - Q_R'*g - Q_R'*G*(Q_X*b_X)
	LinAlgOpPack::V_StMtV( fo, -1.0, Q_R, BLAS_Cpp::trans, g );
	if( num_init_fixed ) {
		SpVector b_XX;
		V_MtV( &b_XX, Q_X, BLAS_Cpp::no_trans, (*b_X)(1,num_init_fixed) );
		SparseLinAlgPack::Vp_StPtMtV( &(*fo)(), -1.0, Q_R, BLAS_Cpp::trans, G, BLAS_Cpp::no_trans, b_XX() );
	}

}

} // end namesapce ConstrainedOptimizationPack