// ////////////////////////////////////////////////////////////////////
// exampleNLPDiagSetup.hpp

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "exampleNLPDiagSetup.hpp"
#include "ExampleVectorLib/src/MPIDenseVector.hpp"
#include "AbstractLinAlgPack/src/serial/interfaces/VectorSpaceSerial.hpp"
#include "AbstractLinAlgPack/src/abstract/tsfcore/VectorSpaceTSFCore.hpp"
#include "TSFCore/src/Core/TSFCoreSerialVectorSpaceDecl.hpp"
#include "OptionsFromStream.hpp"
#include "CommandLineProcessor.hpp"

///
int AbstractLinAlgPack::exampleNLPDiagSetup(
	int argc, char* argv[], MPI_Comm comm
	,MemMngPack::ref_count_ptr<const VectorSpace>   *vec_space
	,size_type *n, value_type *xo, bool *has_bounds, bool *dep_bounded
	)
{

	using std::endl;
	using std::setw;
	namespace mmp = MemMngPack;
	using mmp::ref_count_ptr;
	typedef AbstractLinAlgPack::size_type size_type;
	typedef AbstractLinAlgPack::value_type value_type;

	using AbstractLinAlgPack::VectorSpace;
	using AbstractLinAlgPack::Vector;
	using AbstractLinAlgPack::VectorMutable;

	using AbstractLinAlgPack::VectorSpaceTSFCore;

	using CommandLineProcessorPack::CommandLineProcessor;

	// Get an idea of what processors we have.
	int num_proc, proc_rank;
	MPI_Comm_size( comm, &num_proc );
	MPI_Comm_rank( comm, &proc_rank );

	// Get the size of the problem to solve
	*n = 4;
	// Get the starting point
	*xo = 0.1;
	// Determine if the NLP has bounds or not.
	*has_bounds = false;
	// Make the dependent or independent variables bounded.
	*dep_bounded = true;
	// Serial or parallel?
	bool in_parallel = false;
	// Use TSF?
	bool use_tsf = false;

	CommandLineProcessor  command_line_processor;
	
	command_line_processor.set_option( "n",  n,   "Global number of dependent (and independent) variables" );
	command_line_processor.set_option( "xo", xo,  "Initial guess of the solution" );
	command_line_processor.set_option(
		"has-bounds", "no-has-bounds", has_bounds
		,"Determine if the NLP has bounds or not" );
	command_line_processor.set_option(
		"dep-bounded", "indep-bounded", dep_bounded
		,"Determine if the dependent or independent variables are bounded" );
	command_line_processor.set_option(
		"in-parallel", "in-serial", &in_parallel
		,"Determine if computations are performed in parallel or not" );
	command_line_processor.set_option(
		"use-tsf", "no-use-tsf", &use_tsf
		,"Determine whether to use TSF vectors or not" );
	
	CommandLineProcessor::EParseCommandLineReturn
		parse_return = command_line_processor.parse_command_line(argc,argv,&std::cerr);
	
	if( parse_return != CommandLineProcessor::PARSE_SUCCESSFULL )
		return parse_return;

	// Create the vector space object to use.

	if(in_parallel) {
		//
		// Use parallel vectors!
		//
		// Determine the mapping of elements to processors for MPIDenseVectorSpace
		RTOp_index_type local_dim = (*n)/num_proc; // assume n > num_proc
		RTOp_index_type *ind_map  = new RTOp_index_type[num_proc];
		RTOp_index_type i_u = local_dim;
		for( int p = 0; p < num_proc; ++p, i_u += local_dim )
			ind_map[p] = i_u;
		ind_map[num_proc-1] = *n; // Make sure we don't go past n
		local_dim = ( proc_rank > 0
					  ? ind_map[proc_rank]-ind_map[proc_rank-1]
					  : ind_map[0] );
		*vec_space = mmp::rcp(new MPIDenseVectorSpace(MPI_COMM_WORLD,ind_map,true,1,*n));
	}
	else {
		//
		// Use serial vectors
		//
		if( use_tsf ) {
			*vec_space = mmp::rcp(new VectorSpaceTSFCore(mmp::rcp(new TSFCore::SerialVectorSpace<value_type>(*n))));
		}
		else {
			*vec_space = mmp::rcp(new AbstractLinAlgPack::VectorSpaceSerial(*n));
		}
	}
	
	return 0;
}