// ///////////////////////////////////////////////////////////////////
// TestGenMatrixClass.cpp
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

#include <iomanip>
#include <ostream>
#include <vector>
#include <typeinfo>

#include "TestLinAlgPack.hpp"
#include "LinAlgPack/src/GenMatrixClass.hpp"
#include "LinAlgPack/src/VectorOut.hpp"
#include "LinAlgPack/src/GenMatrixOut.hpp"
#include "LinAlgPack/src/MatVecCompare.hpp"

namespace {

using LinAlgPack::sqrt_eps;

// Check consistency of row(), col(), diag() and operator()().
template<class M_t>
void check_access( M_t& M, typename M_t::size_type row_offset, typename M_t::size_type col_offset
	, std::ostream* out, bool* success )
{
	if(out)
		*out	<< "Checking M(i,j) == M.row(i)(j) == M.col(j)(i) == "
				<< "M.diag(...)(...) == "
				<< "(i + "<<row_offset<<") + 0.1*(j+"<<col_offset<<") : ";

	bool result = true;

	for( typename M_t::size_type i = 1; i <= M.rows(); ++i ) {
		for( typename M_t::size_type j = 1; j <= M.rows(); ++j ) {
			const typename M_t::value_type
				Mij = M(i,j);
			typename M_t::value_type
				val = (i+row_offset)+0.1*(j+col_offset);
			if( ::fabs(Mij-val) > sqrt_eps ) {
				result = false;
				if(out) *out << "(M("<<i<<","<<j<<") -> "<<Mij<<") != "<<val<<std::endl;
			}
			if( Mij != (val = M.row(i)(j)) ) {
				result = false;
				if(out) *out << "M("<<i<<","<<j<<") != (M.row("<<i<<")("<<j<<") -> "<<val<<")\n";
			}
			if( Mij != (val = M.col(j)(i)) ) {
				result = false;
				if(out) *out << "M("<<i<<","<<j<<") != (M.col("<<j<<")("<<i<<") -> "<<val<<")\n";
			}
			const int k = ( i > j ? -i + j : j - i );
			const typename M_t::size_type k_i = ( i > j ? j : i );
			if( Mij != (val = M.diag(k)(k_i) ) ) {
				result = false;
				if(out) *out << "M("<<i<<","<<j<<") != (M.diag("<<k<<")("<<k_i<<") -> "<<val<<")\n";
			}
		}
	}
	if(out) *out << result << std::endl;
	if(!result) *success = false;
}

// Print out a string for overlap
const char* overlap_str( LinAlgPack::EOverLap overlap ) {
	switch(overlap) {
		case LinAlgPack::NO_OVERLAP:
			return "NO_OVERLAP";
		case LinAlgPack::SOME_OVERLAP:
			return "SOME_OVERLAP";
		case LinAlgPack::SAME_MEM:
			return "SAME_MEM";
	}
	return "Invalid value for EOverLap";
}

}	// end namespace

bool LinAlgPack::TestingPack::TestGenMatrixClass(std::ostream* out)
{

	using LinAlgPack::comp;
	using LinAlgPack::sqrt_eps;

	bool success = true;
	bool result;

	if(out)
		*out	<< "\n****************************************************"
				<< "\n*** Testing GenMatrix and GenMatrixSlice classes ***"
				<< "\n****************************************************\n"
				<< std::boolalpha;

	try {

	const size_type
		m = 6,
		n = 8;

	const value_type
		ptr[m*n] =
			 {	1.1,	2.1,	3.1,	4.1,	5.1,	6.1,
				1.2,	2.2,	3.2,	4.2,	5.2,	6.2,
				1.3,	2.3,	3.3,	4.3,	5.3,	6.3,
				1.4,	2.4,	3.4,	4.4,	5.4,	6.4,
				1.5,	2.5,	3.5,	4.5,	5.5,	6.5,
				1.6,	2.6,	3.6,	4.6,	5.6,	6.6,
				1.7,	2.7,	3.7,	4.7,	5.7,	6.7,
				1.8,	2.8,	3.8,	4.8,	5.8,	6.8	};

	// /////////////////////////////
	// Test Constructors

	if(out)
		*out	<< "\n***\n*** Testing constructors\n***\n";

	// GenMatrixSlice
	if(out) *out << "\nGenMatrixSlice gms1;\n";
	GenMatrixSlice gms1;
	if(out) *out << "gms1 =\n" << gms1;
	update_success( result = (gms1.rows() == 0 && gms1.cols() == 0 ), &success );
	if(out)
		*out	<< "((gms1.rows() -> "<<gms1.rows()
				<< ") == 0 && (gms1.cols() -> "<<gms1.cols()<<") == 0 ) : "
				<< result << std::endl;


	if(out) *out << "\nGenMatrixSlice gms2( const_cast<value_type*>(ptr), m*n, m, m, n );\n";
	const GenMatrixSlice gms2( const_cast<value_type*>(ptr), m*n, m, m, n );
	if(out) *out << "gms2 =\n" << gms2;

	if(out) *out << "\nGenMatrixSlice gms3( const_cast<GenMatrixSlice&>(gms2), Range1D(1,m), Range1D(1,n) );\n";
	const GenMatrixSlice gms3( const_cast<GenMatrixSlice&>(gms2), Range1D(1,m), Range1D(1,n) );
	if(out) *out << "gms3 =\n" << gms3;

	// GenMatrix

	if(out) *out << "\nGenMatrix gm1;\n";
	GenMatrix gm1;	
	if(out) *out << "gm1 =\n" << gm1;
	update_success( result = (gm1.rows() == 0 && gm1.cols() == 0 ), &success );
	if(out)
		*out	<< "((gm1.rows() -> "<<gm1.rows()
				<< ") == 0 && (gm1.cols() -> "<<gm1.cols()<<") == 0 ) : "
				<< result << std::endl;
	
	if(out) *out << "\nGenMatrix gm2(m,n);\n";
	GenMatrix gm2(m,n);
	if(out) *out << "gm2 =\n" << gm2;

	if(out) *out << "\nGenMatrix gm3(1.0,m,n);\n";
	GenMatrix gm3(1.0,m,n);
	if(out) *out << "gm3 =\n" << gm3;
	update_success( result = comp( gm3(), 1.0 ), &success );
	if(out) *out << "gm3 == 1.0 : " << result << std::endl;

	if(out) *out << "\nGenMatrix gm4(ptr,m,n);\n";
	GenMatrix gm4(ptr,m,n);
	if(out) *out << "gm4 =\n" << gm4;

	if(out) *out << "\nGenMatrix gm5(gms2);\n";
	GenMatrix gm5(gms2);
	if(out) *out << "gm5 =\n" << gm5;

	// ////////////////////////////
	// Test GenMatrixSlice binding

	if(out)
		*out	<< "\n***\n*** Testing GenMatrixSlice binding\n***\n";

	if(out) *out << "\ngms1.bind(gm4());\n";
	gms1.bind(gm4());
	if(out) *out << "gms1 =\n" << gms1;

	// ////////////////////////////
	// Test GenMatrix resizing

	if(out)
		*out	<< "\n***\n*** Testing GenMatrix resizing\n***\n";

	if(out) *out << "\ngm1.resize(m,n,1.0);\n";
	gm1.resize(m,n,1.0);
	if(out) *out << "gm1 =\n" << gm1;
	update_success( result = comp( gm1(), 1.0 ), &success );
	if(out) *out << "gm1 == 1.0 : " << result << std::endl;

	// ///////////////////////////////////////////////
	// Test row, col, diag access and element access

	// GenMatrixSlice

	if(out)
		*out	<< "\n***\n*** Testing row, col, diag access and element access\n***\n";

	if(out) *out << "\nLet M = gms1\n";
	check_access( gms1, 0, 0, out, &success );

	if(out) *out << "\nLet M = const_cast<const GenMatrixSlice&>(gms1)\n";
	check_access( const_cast<const GenMatrixSlice&>(gms1), 0, 0, out, &success );

	// GenMatrix

	if(out) *out << "\nLet M = gm4\n";
	check_access( gm4, 0, 0, out, &success );

	if(out) *out << "\nLet M = const_cast<const GenMatrix&>(gm4)\n";
	check_access( const_cast<const GenMatrix&>(gm4), 0, 0, out, &success );

	// ////////////////////////////
	// Test submatrix access

	if(out)
		*out	<< "\n***\n*** Testing submatrix access\n***\n";

	if(out) *out << "\nRange1D r_rng(2,m-1), c_rng(2,n-1);\n";
	Range1D r_rng(2,m-1), c_rng(2,n-1);

	// GenMatrixSlice

	if(out) *out << "\nLet M = const_cast<GenMatrixSlice&>(gms2)(r_rng,c_rng)\n";
	gms1.bind( const_cast<GenMatrixSlice&>(gms2)(r_rng,c_rng) );
	if(out) *out << "M =\n" << gms1;
	check_access( gms1, 1, 1, out, &success );

	if(out) *out << "\nLet M = const_cast<GenMatrixSlice&>(gms2(r_rng,c_rng))\n";
	gms1.bind( const_cast<GenMatrixSlice&>(gms2)(r_rng,c_rng) );
	if(out) *out << "M =\n" << gms1;
	check_access( gms1, 1, 1, out, &success );

	if(out) *out << "\nLet M = const_cast<GenMatrixSlice&>(gms2)(2,m-1,2,n-1)\n";
	gms1.bind(const_cast<GenMatrixSlice&>(gms2)(2,m-1,2,n-1) );
	if(out) *out << "M =\n" << gms1;
	check_access( gms1, 1, 1, out, &success );

	if(out) *out << "\nLet M = const_cast<GenMatrixSlice&>(gms2(2,m-1,2,n-1))\n";
	gms1.bind( const_cast<GenMatrixSlice&>(gms2)(2,m-1,2,n-1) );
	if(out) *out << "M =\n" << gms1;
	check_access( gms1, 1, 1, out, &success );

	// GenMatrix

	if(out) *out << "\nLet M = gm4(r_rng,c_rng)\n";
	gms1.bind( gm4(r_rng,c_rng) );
	if(out) *out << "M =\n" << gms1;
	check_access( gms1, 1, 1, out, &success );

	if(out) *out << "\nLet M = const_cast<const GenMatrixSlice&>(gm4)(r_rng,c_rng)\n";
	gms1.bind( const_cast<const GenMatrix&>(gm4)(r_rng,c_rng) );
	if(out) *out << "M =\n" << gms1;
	check_access( gms1, 1, 1, out, &success );

	if(out) *out << "\nLet M = gm4(2,m-1,2,n-1)\n";
	gms1.bind( gm4(2,m-1,2,n-1) );
	if(out) *out << "M =\n" << gms1;
	check_access( gms1, 1, 1, out, &success );

	if(out) *out << "\nLet M = const_cast<const GenMatrixSlice&>(gm4)(2,m-1,2,n-1)\n";
	gms1.bind( const_cast<const GenMatrix&>(gm4)(2,m-1,2,n-1) );
	if(out) *out << "M =\n" << gms1;
	check_access( gms1, 1, 1, out, &success );

	// ////////////////////
	// Test matrix overlap

	if(out)
		*out	<< "\n***\n*** matrix overlap\n***\n";

	EOverLap ovlap;

	// GenMatrixSlice

	if(out) *out << "(gms2.overlap(gms2) -> ";
	ovlap = gms2.overlap(gms2);
	result = update_success( ovlap == SAME_MEM, &success );
	if(out)	*out	<< overlap_str(ovlap) << ") == SAME_MEM : " << result << std::endl;

	if(out) *out << "(gms2.overlap(gms2(r_rng,c_rng)) -> ";
	ovlap = gms2.overlap(gms2(r_rng,c_rng));
	result = update_success( ovlap == SOME_OVERLAP, &success );
	if(out)	*out	<< overlap_str(ovlap) << ") == SOME_OVERLAP : " << result << std::endl;

	if(out) *out << "(gms2(1,m/2,1,n/2).overlap(gms2(m/2,m,n/2,n)) -> ";
	ovlap = gms2(1,m/2,1,n/2).overlap(gms2(m/2,m,n/2,n));
	result = update_success( ovlap == SOME_OVERLAP, &success );
	if(out)	*out	<< overlap_str(ovlap) << ") == SOME_OVERLAP : " << result << std::endl;

	if(out) *out << "(gms2(1,m/2,1,n/2).overlap(gms2(m/2+1,m,n/2+1,n)) -> ";
	ovlap = gms2(1,m/2,1,n/2).overlap(gms2(m/2+1,m,n/2+1,n));
	result = update_success( ovlap == NO_OVERLAP, &success );
	if(out)	*out	<< overlap_str(ovlap) << ") == NO_OVERLAP : " << result << std::endl;

	// GenMatrix

	if(out) *out << "(gm4.overlap(gm4) -> ";
	ovlap = gm4.overlap(gm4);
	result = update_success( ovlap == SAME_MEM, &success );
	if(out)	*out	<< overlap_str(ovlap) << ") == SAME_MEM : " << result << std::endl;

	if(out) *out << "(gm4.overlap(gm4(r_rng,c_rng)) -> ";
	ovlap = gm4.overlap(gm4(r_rng,c_rng));
	result = update_success( ovlap == SOME_OVERLAP, &success );
	if(out)	*out	<< overlap_str(ovlap) << ") == SOME_OVERLAP : " << result << std::endl;

	// //////////////////////////////////////////////////////
	// Test vector overlap (continuation of vector testing)
	//
	// ToDo: Finish this someday once you get it figured out.

	// ///////////////////////////
	// Test assignment operators

	if(out)
		*out	<< "\n***\n*** assignment operators\n***\n";

	// GenMatrixSlice

	if(out) *out << "\ngms1.bind(gm1());\n";
	gms1.bind(gm1());

	if(out) *out << "\ngms1 = 2.0;\n";
	gms1 = 2.0;
	if(out) *out << "gms1 =\n" << gms1;
	update_success( result = comp(gms1,2.0), &success );
	if(out) *out << "gms1 == 2.0 : " << result << std::endl; 

	if(out) *out << "\ngms1 = gms2;\n";
	gms1 = gms2;
	if(out) *out << "gms1 =\n" << gms1;
	update_success( result = comp(gms1,gms2), &success );
	if(out) *out << "gms1 == gms2 : " << result << std::endl; 

	// GenMatrix

	if(out) *out << "\ngm1 = 3.0;\n";
	gm1 = 3.0;
	if(out) *out << "gm1 =\n" << gm1;
	update_success( result = comp(gm1,3.0), &success );
	if(out) *out << "gm1 == 3.0 : " << result << std::endl; 

	if(out) *out << "\ngm1.resize(0,0); gm1 = gms2;\n";
	gm1.resize(0,0);
	gm1 = gms2;
	if(out) *out << "gm1 =\n" << gm1;
	update_success( result = comp(gm1,gms2), &success );
	if(out) *out << "gm1 == gms2 : " << result << std::endl; 

	if(out) *out << "\ngm1.resize(0,0); gm1 = gm4;\n";
	gm1.resize(0,0);
	gm1 = gm4;
	if(out) *out << "gm1 =\n" << gm1;
	update_success( result = comp(gm1,gm4), &success );
	if(out) *out << "gm1 == gm4 : " << result << std::endl; 

	} // end try
	catch( const std::exception& excpt ) {
		success = false;
		if(out)
			(*out)	<< "\nError, a standard exception was thrown: "
					<< typeid(excpt).name() << ": "
					<< excpt.what() << std::endl; 
	}
	catch(...) {
		success = false;
		if(out)
			(*out)	<< "\nError, an unknown exception was thrown\n";
	}

	if(out) {
		if(success)
			(*out)
				<< "\n*** Congradulations, GenMatrix and GenMatrixSlice seem to check out. ***\n";
		else
			(*out)
				<< "\n*** Oops, all of the tests for GenMatrix and GenMatrixSlice "
					"where not successful. ***\n";
	}

	return success;
}

