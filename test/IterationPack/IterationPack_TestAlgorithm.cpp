// /////////////////////////////////////////////////////////////
// TestAlgorithm.cpp
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

#pragma warning(disable : 4786)	// too long class name for debugger warning

#include <ostream>
#include <iomanip>

#include "../test/TestGeneralIterationPack.h"
#include "../include/Algorithm.h"
#include "../include/AlgorithmState.h"
#include "../test/AlgorithmStepTesting.h"
#include "../test/AlgorithmTrackTesting.h"
#include "debug.h"

namespace GeneralIterationPack {

// Implement minor loop 1
class MinorLoop1Step : public AlgorithmStepTesting {
public:

	MinorLoop1Step() : called_(false)
	{}

	bool do_step(Algorithm& algo, poss_type step_poss, EDoStepType type
		, poss_type assoc_step_poss)
	{
		if(called_) {
			algo.terminate(false);
			return false;
		}
		AlgorithmStepTesting::do_step(algo,step_poss,type,assoc_step_poss);
		algo.track().journal_out() << "\n* Jump to \"Step_2\"\n";
		algo.do_step_next("Step_2");
		called_ = true;
		return false;
	}

private:
	bool called_;
};

// Implement controled loop 1
class ControledLoop1Step : public AlgorithmStepTesting {
public:

	ControledLoop1Step() : called_(false)
	{}

	bool do_step(Algorithm& algo, poss_type step_poss, EDoStepType type
		, poss_type assoc_step_poss)
	{
		if(called_) {
			algo.terminate(false);
			return false;
		}
		AlgorithmStepTesting::do_step(algo,step_poss,type,assoc_step_poss);
		algo.track().journal_out() << "\n* Do step \"Step_2_p1\" then Jump to \"Step_1\"\n";
		algo.get_assoc_step(2,POST_STEP,1)->do_step(algo,2,DO_POST_STEP,1);
		algo.do_step_next("Step_1");
		called_ = true;
		return false;
	}

private:
	bool called_;
};

// Runtime Config change
class RuntimeConfigChangeStep : public AlgorithmStepTesting {
public:

	RuntimeConfigChangeStep() : called_(false)
	{}

	bool do_step(Algorithm& algo, poss_type step_poss, EDoStepType type
		, poss_type assoc_step_poss)
	{
		AlgorithmStepTesting::do_step(algo,step_poss,type,assoc_step_poss);
		if(!called_) {
			algo.track().journal_out() << "\n* Remove step \"Step_1\" then put it back then print steps\n";
			algo.begin_config_update();
			Algorithm::step_ptr_t step = algo.get_step(1);
			algo.remove_step(1);
			algo.insert_step(1,"Step_1",step);
			algo.end_config_update();
			called_ = true;
		}
		else {
			algo.track().journal_out() << "\n* Remove step 1, step_name = "
				<< algo.get_step_name(1) << " but don't put it back then print steps\n";
			algo.begin_config_update();
			algo.remove_step(1);
			algo.end_config_update();
			called_ = true;
		}
		algo.print_steps(algo.track().journal_out());
		return true;
	}

private:
	bool called_;
};

}	// end namespace GeneralIterationPack

bool GeneralIterationPack::TestingPack::TestAlgorithm(std::ostream* out) {

	using std::endl;
	using std::setw;
	
	std::ostream& _out = *out;
	// ToDo: RAB: 7/1/99: Modify for optional output when out == 0;

	try {

	bool success = true;
//	const int w = 15;
	_out << std::boolalpha;
		
	_out	<< "\n\n*************************\n"
			<< "*** Testing Algorithm ***\n"
			<< "*************************\n";

	Algorithm algo;
	
	Algorithm::state_ptr_t			state		= new AlgorithmState;
	Algorithm::track_ptr_t			track		= new AlgorithmTrackTesting(_out);

	algo.set_state( state );
	algo.set_track( track );

	Algorithm::step_ptr_t			step		= new AlgorithmStepTesting;
	Algorithm::step_ptr_t			assoc_step	= new AlgorithmStepTesting;

	algo.insert_step( 1, "Step_1", step );

	algo.insert_step( 2, "Step_2", step );

	algo.insert_assoc_step( 2, PRE_STEP, 1, "Step_2_m1", assoc_step );

	algo.insert_assoc_step( 2, POST_STEP, 1, "Step_2_p1", assoc_step );

	algo.insert_step( 3, "Step_3", step );

	algo.insert_step( 4, "Step_4", step );

	algo.insert_assoc_step( 4, PRE_STEP, 1, "Step_4_m2", assoc_step );

	algo.insert_assoc_step( 4, PRE_STEP, 2, "Step_4_m1", assoc_step );

	_out	<< "\n\n*** algo.print_algorithm(_out)\n\n";
	algo.print_algorithm(_out);

	// Test the major loop

	_out	<< "\n\n*** Test the major loop for two iterations ***\n";

	_out	<< "\nalgo.max_iter(2); algo.do_algorithm();\n\n";
	algo.max_iter(2);
	algo.do_algorithm();

	// Test Minor Loop 1

	_out	<< "\n\n*** Test Minor Loop 1 ***\n";

	_out	<< "\nalgo.remove_assoc_step( 4, PRE_STEP, 2 );\n";
	algo.remove_assoc_step( 4, PRE_STEP, 2 );

	_out	<< "\nalgo.insert_assoc_step( 4, PRE_STEP, 2, \"Step_4_m1\", new MinorLoop1Step );\n";
	algo.insert_assoc_step( 4, PRE_STEP, 2, "Step_4_m1", new MinorLoop1Step );

	_out	<< "\nalgo.state().k(0);\n";
	algo.state().k(0);

	_out	<< "\nalgo.max_iter(1);\n";
	algo.max_iter(1);

	_out	<< "\n\nalgo.print_steps(_out)\n\n";
	algo.print_steps(_out);

	_out	<< "\nalgo.do_algorithm();\n";
	algo.do_algorithm();

	// Test Controlled Loop 1

	_out	<< "\n\n*** Test Controlled Loop 1 ***\n";

	_out	<< "\nalgo.remove_assoc_step( 4, PRE_STEP, 2 );\n";
	algo.remove_assoc_step( 4, PRE_STEP, 2 );

	_out	<< "\nalgo.insert_assoc_step( 4, PRE_STEP, 2, \"Step_4_m1\", assoc_step );\n";
	algo.insert_assoc_step( 4, PRE_STEP, 2, "Step_4_m1", assoc_step );

	_out	<< "\nalgo.remove_assoc_step( 4, PRE_STEP, 1 );\n";
	algo.remove_assoc_step( 4, PRE_STEP, 1 );

	_out	<< "\nalgo.insert_assoc_step( 4, PRE_STEP, 1 , \"Step_4_m2\", new ControledLoop1Step );\n";
	algo.insert_assoc_step( 4, PRE_STEP, 1 , "Step_4_m2", new ControledLoop1Step );

	_out	<< "\n\nalgo.print_steps(_out)\n\n";
	algo.print_steps(_out);

	_out	<< "\nalgo.state.k(0);\n";
	algo.state().k(0);

	_out	<< "\nalgo.do_algorithm();\n";
	algo.do_algorithm();

	// Test runtime configuration change.

	_out	<< "\n\n*** Test runtime configuration change ***\n";

	_out	<< "\nalgo.remove_assoc_step( 4, PRE_STEP, 1 );\n";
	algo.remove_assoc_step( 4, PRE_STEP, 1 );

	_out	<< "\nalgo.insert_assoc_step( 4, PRE_STEP, 1, \"Step_4_m2\", new  RuntimeConfigChangeStep );\n";
	algo.insert_assoc_step( 4, PRE_STEP, 1, "Step_4_m2", new  RuntimeConfigChangeStep );

	_out	<< "\n\nalgo.print_steps(_out)\n\n";
	algo.print_steps(_out);

	_out	<< "\nalgo.state.k(0);\n";
	algo.state().k(0);

	_out	<< "\nalgo.max_iter(5);\n";
	algo.max_iter(5);

	try {
		_out	<< "* In the 5th iteration a Algorithm::InvalidConfigChange exception should be thrown.";
		_out	<< "\nalgo.do_algorithm();\n";
		algo.do_algorithm();
		success = false;
	}
	catch(Algorithm::InvalidConfigChange& excpt) {
		_out << "\nCaught a Algorithm::InvalidConfigChange&: " << excpt.what() << endl;
	}

	_out << "\n*** Congradulations, If you read this the tests for Algorithm"
			" seem to have been successful\n";

	return success;

	} // end try
	catch(const std::exception& excpt) {
		_out << "\nCaught a std::exception: " << excpt.what() << endl;
	}
	catch(...) {
		_out << "\nCaught an unknown exception\n";
	}

	return false;
}
