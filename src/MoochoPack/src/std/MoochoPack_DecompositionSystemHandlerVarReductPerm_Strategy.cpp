// ////////////////////////////////////////////////////////////////////////////
// DecompositionSystemHandlerVarReductPerm_Strategy.cpp
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

#ifndef MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS

#include <ostream>
#include <typeinfo>

#include "ReducedSpaceSQPPack/src/std/DecompositionSystemHandlerVarReductPerm_Strategy.hpp"
#include "ReducedSpaceSQPPack/src/ReducedSpaceSQPPackExceptions.hpp"
#include "ReducedSpaceSQPPack/src/rsqp_algo_conversion.hpp"
#include "IterationPack/src/print_algorithm_step.hpp"
#include "ConstrainedOptimizationPack/src/DecompositionSystemVarReductPerm.hpp"
#include "NLPInterfacePack/src/NLPFirstOrderInfo.hpp"
#include "NLPInterfacePack/src/NLPVarReductPerm.hpp"
#include "AbstractLinAlgPack/src/PermutationOut.hpp"
#include "AbstractLinAlgPack/src/MatrixWithOpNonsingular.hpp"
#include "AbstractLinAlgPack/src/MatrixWithOpOut.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpMutable.hpp"
#include "AbstractLinAlgPack/src/VectorStdOps.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpOut.hpp"
#include "AbstractLinAlgPack/src/assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack/src/LinAlgOpPack.hpp"
#include "dynamic_cast_verbose.hpp"
#include "ThrowException.hpp"

namespace ReducedSpaceSQPPack {

// Constructors / initializers

DecompositionSystemHandlerVarReductPerm_Strategy::DecompositionSystemHandlerVarReductPerm_Strategy()
{}

// Overridden from DecompositionSystemHandler_Strategy

bool DecompositionSystemHandlerVarReductPerm_Strategy::update_decomposition(
	rSQPAlgo                                &algo
	,rSQPState                              &s
	,NLPFirstOrderInfo                      &nlp
	,EDecompSysTesting                      decomp_sys_testing
	,EDecompSysPrintLevel                   decomp_sys_testing_print_level
	,bool                                   *new_decomp_selected
	)
{
	using DynamicCastHelperPack::dyn_cast;

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out << "Updating the decomposition ...\n";
	}

	DecompositionSystemVarReductPerm   &decomp_sys_perm = dyn_cast<DecompositionSystemVarReductPerm>(s.decomp_sys());
	NLPVarReductPerm                   &nlp_vrp         = dyn_cast<NLPVarReductPerm>(nlp);

	const size_type
		n  = nlp.n(),
		nb = nlp.num_bounded_x(),
		m  = nlp.m(),
		mI = nlp.mI();
	size_type
		r  = s.decomp_sys().equ_decomp().size();

	bool decomp_updated     = false;
	bool get_new_basis      = false;
	bool new_basis_selected = false;

	if( select_new_decomposition_ ) {
		if( olevel >= PRINT_ALGORITHM_STEPS )
			out << "\nSome client called select_new_decomposition() so we will select a new basis ...\n";
		get_new_basis = true;
		select_new_decomposition_ = false;
	}
	else if( !decomp_sys_perm.has_basis() ) {
		if( olevel >= PRINT_ALGORITHM_STEPS )
			out << "\nDecompositionSystemVarReductPerm object currently does not have a basis so we must select one ...\n";
		get_new_basis = true;
	}

	// Get the iteration quantity container objects
	IterQuantityAccess<index_type>
		&num_basis_iq = s.num_basis();
	IterQuantityAccess<VectorWithOpMutable>
		&x_iq   = s.x(),
		&nu_iq  = s.nu();
	IterQuantityAccess<MatrixWithOp>
		*Gc_iq  = m  > 0                                   ? &s.Gc() : NULL,
		*Gh_iq  = mI > 0                                   ? &s.Gh() : NULL,
		*Z_iq   = ( n > m && r > 0 )    || get_new_basis   ? &s.Z()  : NULL,
		*Y_iq   = ( r > 0 )             || get_new_basis   ? &s.Y()  : NULL,
		*Uz_iq  = ( m  > 0 && m  > r )  || get_new_basis   ? &s.Uz() : NULL,
		*Uy_iq  = ( m  > 0 && m  > r )  || get_new_basis   ? &s.Uy() : NULL,
		*Vz_iq  = ( mI > 0 ) && ( m > 0 || get_new_basis ) ? &s.Vz() : NULL,
		*Vy_iq  = ( mI > 0 ) && ( m > 0 || get_new_basis ) ? &s.Vy() : NULL;
	IterQuantityAccess<MatrixWithOpNonsingular>
		*R_iq   = ( m > 0 )                                ? &s.R()  : NULL;
	
	//
	// Update (or select a new) range/null decomposition
	//

	// Determine if we will test the decomp_sys or not
	const DecompositionSystem::ERunTests
		ds_test_what = ( ( decomp_sys_testing == DST_TEST
						   || ( decomp_sys_testing == DST_DEFAULT
								&& algo.algo_cntr().check_results() ) )
						 ? DecompositionSystem::RUN_TESTS
						 : DecompositionSystem::NO_TESTS );
		
	// Determine the output level for decomp_sys				
	DecompositionSystem::EOutputLevel ds_olevel;
	switch(olevel) {
		case PRINT_NOTHING:
		case PRINT_BASIC_ALGORITHM_INFO:
			ds_olevel = DecompositionSystem::PRINT_NONE;
			break;
		case PRINT_ALGORITHM_STEPS:
		case PRINT_ACTIVE_SET:
			ds_olevel = DecompositionSystem::PRINT_BASIC_INFO;
			break;
		case PRINT_VECTORS:
			ds_olevel = DecompositionSystem::PRINT_VECTORS;
			break;
		case PRINT_ITERATION_QUANTITIES:
			ds_olevel = DecompositionSystem::PRINT_EVERY_THING;
			break;
		default:
			assert(0); // Should not get here!
	};

	if( !get_new_basis ) {
		// Form the decomposition of Gc and Gh and update the decomposition system matrices
		if( olevel >= PRINT_ALGORITHM_STEPS ) {
			out << "\nUpdating the range/null decompostion matrices ...\n";
		}
		try {
			s.decomp_sys().update_decomp(
				&out                               // out
				,ds_olevel                         // olevel
				,ds_test_what                      // test_what
				,Gc_iq->get_k(0)                   // Gc
				,Gh_iq ? &Gh_iq->get_k(0) : NULL   // Gh
				,&Z_iq->set_k(0)                   // Z
				,&Y_iq->set_k(0)                   // Y
				,&R_iq->set_k(0)                   // R
				,Uz_iq ? &Uz_iq->set_k(0) : NULL   // Uz
				,Uy_iq ? &Uy_iq->set_k(0) : NULL   // Uy
				,Vz_iq ? &Vz_iq->set_k(0) : NULL   // Vz
				,Vy_iq ? &Vy_iq->set_k(0) : NULL   // Vy
				,DecompositionSystem::MATRICES_ALLOW_DEP_IMPS // ToDo: Change this!
				);
			s.equ_decomp(   s.decomp_sys().equ_decomp()   );
			s.equ_undecomp( s.decomp_sys().equ_undecomp() );
			decomp_updated = true;
		}
		catch( const DecompositionSystem::SingularDecomposition& except) {
			if( olevel >= PRINT_BASIC_ALGORITHM_INFO )
				out
					<< "\nOops!  This decomposition was singular; must select a new basis!\n"
					<< except.what() << std::endl;
		}
	}

	if( !decomp_updated ) {
		if( s.get_P_var_current().get() == NULL ) {
			MemMngPack::ref_count_ptr<Permutation>
				P_var = nlp_vrp.factory_P_var()->create(),
				P_equ = nlp_vrp.factory_P_equ()->create();
			Range1D
				var_dep,
				equ_decomp;
			nlp_vrp.get_basis(
				P_var.get(), &var_dep, P_equ.get(), &equ_decomp, NULL, NULL );
			s.set_P_var_current( P_var );
			s.set_P_equ_current( P_equ );
		}
		MemMngPack::ref_count_ptr<Permutation>
			P_var = nlp_vrp.factory_P_var()->create(),
			P_equ = nlp_vrp.factory_P_equ()->create();
		Range1D
			var_dep,
			equ_decomp;
		bool nlp_selected_basis = false;
		bool do_permute_x = true;
		if( nlp_vrp.nlp_selects_basis() ) {
			// The nlp may select the new (or first) basis.
				
			// If this is the first basis, the NLPVarReductPerm interface specifies that it
			// will already be set for the nlp.  Check to see if this is the first basis
			// and if not, ask the nlp to give you the next basis.
			// I must form a loop here to deal with the
			// possibility that the basis the nlp selects will be singular.
			if( olevel >= PRINT_BASIC_ALGORITHM_INFO )
				out
					<< "\nThe NLP will attempt to select a basis "
					<< "(k = " << s.k() << ")...\n";
			// If decomp_sys_per->has_basis() == false, the first execution of the while()
			// statement will not execute get_next_basis(...).		
			bool very_first_basis = !decomp_sys_perm.has_basis();
			bool a_basis_was_singular = false;
			if(very_first_basis)
				nlp_vrp.get_basis(
					P_var.get(), &var_dep, P_equ.get(), &equ_decomp, NULL, NULL );
			while( very_first_basis
				   || nlp_vrp.get_next_basis(
					   P_var.get(), &var_dep, P_equ.get(), &equ_decomp, NULL, NULL )
				)
			{
				try {
					very_first_basis = false;
					decomp_sys_perm.set_decomp(
						&out                               // out
						,ds_olevel                         // olevel
						,ds_test_what                      // test_what
						,*P_var                            // P_var
						,var_dep                           // var_dep
						,P_equ.get()                       // P_equ
						,&equ_decomp                       // equ_decomp
						,Gc_iq->get_k(0)                   // Gc
						,Gh_iq ? &Gh_iq->get_k(0) : NULL   // Gh
						,&Z_iq->set_k(0)                   // Z
						,&Y_iq->set_k(0)                   // Y
						,&R_iq->set_k(0)                   // R
						,Uz_iq ? &Uz_iq->set_k(0) : NULL   // Uz
						,Uy_iq ? &Uy_iq->set_k(0) : NULL   // Uy
						,Vz_iq ? &Vz_iq->set_k(0) : NULL   // Vz
						,Vy_iq ? &Vy_iq->set_k(0) : NULL   // Vy
						,DecompositionSystem::MATRICES_ALLOW_DEP_IMPS // ToDo: Change this to MATRICES_INDEP_IMPS
						);
					// If you get here the basis was not singular.
					nlp_selected_basis = true;
					break; // break out of the while(...) loop
				}
				// Catch the singularity exceptions and loop around
				catch( const DecompositionSystem::SingularDecomposition& except )
				{
					if( olevel >= PRINT_BASIC_ALGORITHM_INFO )
						out
							<< "\nOops!  This decomposition was singular; ask the NLP for another basis!\n"
							<< except.what() << std::endl;
					a_basis_was_singular = true;
				}
				// Any other exception gets thrown clean out of here.
			}
			do_permute_x = 	!( very_first_basis && !a_basis_was_singular );
			if( olevel >= PRINT_BASIC_ALGORITHM_INFO && !nlp_selected_basis )
				out
					<< "\nThe NLP was unable to provide a nonsigular basis "
					<< "(k = " << s.k() << ")\n";
		}
		if(!nlp_selected_basis) {
			// If you get into here then the nlp could not select a nonsingular
			// basis so we will let the decomposition system select a basis.
			// and give it to the nlp.
			
			if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
				out
					<< "\nThe decomposition system object is selecting the basis "
					<< "(k = " << s.k() << ")...\n";
			}
			decomp_sys_perm.select_decomp(
				&out                                        // out
				,ds_olevel                                  // olevel
				,ds_test_what                               // test_what
				,nu_iq.updated_k(0)?&nu_iq.get_k(0):NULL    // nu
				,&Gc_iq->get_k(0)                           // Gc
				,Gh_iq ? &Gh_iq->get_k(0) : NULL            // Gh
				,P_var.get()                                // P_var
				,&var_dep                                   // var_dep
				,P_equ.get()                                // P_equ
				,&equ_decomp                                // equ_decomp
				,&Z_iq->set_k(0)                            // Z
				,&Y_iq->set_k(0)                            // Y
				,&R_iq->set_k(0)                            // R
				,Uz_iq ? &Uz_iq->set_k(0) : NULL            // Uz
				,Uy_iq ? &Uy_iq->set_k(0) : NULL            // Uy
				,Vz_iq ? &Vz_iq->set_k(0) : NULL            // Vz
				,Vy_iq ? &Vy_iq->set_k(0) : NULL            // Vy
				,DecompositionSystem::MATRICES_ALLOW_DEP_IMPS // ToDo: Change this to MATRICES_INDEP_IMPS
				);
			nlp_vrp.set_basis(	*P_var, var_dep, P_equ.get(), &equ_decomp, NULL, NULL );
		}
		
		// If you get here then no unexpected exceptions where thrown and a new
		// basis has been selected.
				
		new_basis_selected = true;
		r                  = s.decomp_sys().equ_decomp().size();

		// Record this basis change

		const int
			last_updated_k = num_basis_iq.last_updated();
		const index_type
			num_basis = ( last_updated_k != IterQuantity::NONE_UPDATED ? num_basis_iq.get_k(last_updated_k) : 0 ) + 1;
		num_basis_iq.set_k(0) = num_basis;

		s.var_dep(      decomp_sys_perm.var_dep()      );
		s.var_indep(    decomp_sys_perm.var_indep()    );
		s.equ_decomp(   decomp_sys_perm.equ_decomp()   );
		s.equ_undecomp( decomp_sys_perm.equ_undecomp() );
					
		s.set_P_var_last( s.get_P_var_current() );
		s.set_P_equ_last( s.get_P_equ_current() );

		s.set_P_var_current( P_var );
		s.set_P_equ_current( P_equ );

		if( olevel >= PRINT_VECTORS ) {
			out	<< "\nNew permutations:"
				<< "\nP_var_current() =\n" << s.P_var_current()
				<< "\nP_equ_current() =\n" << s.P_equ_current();
		}
					
		if( do_permute_x ) {
			// Sort x according to this new basis.
			VectorWithOpMutable &x = x_iq.get_k(0);
			s.P_var_last().permute( BLAS_Cpp::trans, &x ); // Permute back to original order
			if( olevel >= PRINT_VECTORS ) {
				out	<< "\nx resorted to the original order\n" << x;
			}
			s.P_var_current().permute( BLAS_Cpp::no_trans, &x ); // Permute to new (current) order
			if( olevel >= PRINT_VECTORS ) {
				out	<< "\nx resorted to new basis \n" << x;
			}
		}
					
		// Set the new range and null spaces (these will update all of the set vectors!)
		s.set_space_range( decomp_sys_perm.space_range() );
		s.set_space_null(  decomp_sys_perm.space_null()  );

	}

	*new_decomp_selected = new_basis_selected;

	return true;
}

void DecompositionSystemHandlerVarReductPerm_Strategy::print_update_decomposition(
	const rSQPAlgo                          &algo
	,const rSQPState                        &s
	,std::ostream                           &out
	,const std::string                      &L
	) const
{
	out
		<< L << "*** Updating or selecting a new decomposition using a variable reduction\n"
		<< L << "*** range/null decomposition object.\n"
		<< L << "if decomp_sys does not support the DecompositionSystemVarReductPerm interface then throw exception!\n"
		<< L << "if nlp does not support the NLPVarReductPerm interface then throw exception!\n"
		<< L << "decomp_updated     = false\n"
		<< L << "get_new_basis      = false\n"
		<< L << "new_basis_selected = false\n"
		<< L << "if( select_new_decomposition == true ) then\n"
		<< L << "  get_new_basis            = true\n"
		<< L << "  select_new_decomposition = false\n"
		<< L << "end\n"
		<< L << "if (decomp_sys does not have a basis) then\n"
		<< L << "  get_new_basis = true\n"
		<< L << "end\n"
		<< L << "if (get_new_basis == true) then\n"
		<< L << "  begin update decomposition\n"
		<< L << "  (class = \'" << typeid(s.decomp_sys()).name() << "\')\n"
		;
	s.decomp_sys().print_update_decomp( out, L + "    " );
	out
		<< L << "  end update decomposition\n"
		<< L << "if SingularDecomposition exception was not thrown then\n"
		<< L << "  decomp_updated = true\n"
		<< L << "end\n"
		<< L << "if (decomp_updated == false) then\n"
		<< L << "  nlp_selected_basis = false\n"
		<< L << "  if (nlp.selects_basis() == true) then\n"
		<< L << "    for each basis returned from nlp.get_basis(...) or nlp.get_next_basis()\n"
		<< L << "      decomp_sys.set_decomp(Gc_k,Gh_k...) -> Z_k,Y_k,R_k,Uz_k,Uy_k,Vz_k,Vy_k \n"
		<< L << "      if SingularDecompositon exception was not thrown then\n"
		<< L << "        nlp_selected_basis = true\n"
		<< L << "        exit loop\n"
		<< L << "      end\n"
		<< L << "    end\n"
		<< L << "  end\n"
		<< L << "  if (nlp_selected_basis == false) then\n"
		<< L << "    decomp_sys.select_decomp(Gc_k,Gh_k...) -> P_var,P_equ,Z,Y,R,Uz,Uy,Vz,Vy\n"
		<< L << "                                              and permute Gc,Gh\n"
		<< L << "  end\n"
		<< L << "  *** If you get here then no unexpected exceptions were thrown and a new\n"
		<< L << "  *** basis has been selected\n"
		<< L << "  num_basis_k = num_basis_k(last_updated) + 1\n"
		<< L << "  P_var_last = P_var_current\n"
		<< L << "  P_equ_last = P_equ_current\n"
		<< L << "  P_var_current = P_var\n"
		<< L << "  P_equ_current = P_equ\n"
		<< L << "  Resort x_k according to P_var_current\n"
		<< L << "end\n"
		;
}

// Overridden from DecompositionSystemHandlerSelectNew_Strategy

void DecompositionSystemHandlerVarReductPerm_Strategy::select_new_decomposition( bool select_new_decomposition )
{
	select_new_decomposition_ = select_new_decomposition;
}

} // end namespace ReducedSpaceSQPPack

#endif // MOOCHO_NO_BASIS_PERM_DIRECT_SOLVERS
