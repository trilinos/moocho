// ///////////////////////////////////////////////////////////////////////
// VectorSpaceTester.h
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

#ifndef VECTOR_SPACE_TESTER_H
#define VECTOR_SPACE_TESTER_H

#include <iosfwd>

#include "AbstractLinAlgPack/include/AbstractLinAlgPackTypes.h"
#include "StandardMemberCompositionMacros.h"

namespace AbstractLinAlgPack {

///
/** Testing class for \c VectorSpace, \c VectorWithOp and \c VectorWithOpMutable.
 *
 * The purpose of this class is to test a \c VectorSpace object and the
 * \c VectorWithOpMutable objects that it creates.  The testing function
 * \c check_vector_space() calls all of the methods defined in the interfaces
 * \c %VectorSpace, \c %VectorWithOp and \c %VectorWithOpMutable and checks
 * many of the post conditions but not all.  It would be very difficult to 
 * completely verify every postcondition in every situation. 
 *
 * The behavior of the testing function check_vector_space() is strongly influenced
 * by a set of options (see \c VectorSpaceTester()).
 *
 * When writting new vector implementations, a developer is likely to spend a lot
 * of time debuggin while in this testing function.
 */
class VectorSpaceTester {
public:

	/// Members for option \c print_all_tests() (see StandardMemberCompositionMacros.h).
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, print_all_tests )
#ifdef DOXYGEN_COMPILE
		;
#endif		
	/// Members for option \c print_vectors() (see StandardMemberCompositionMacros.h).
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, print_vectors )
#ifdef DOXYGEN_COMPILE
		;
#endif		
	/// Members for option \c throw_exception() (see StandardMemberCompositionMacros.h).
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, throw_exception )
#ifdef DOXYGEN_COMPILE
		;
#endif		
	/// Members for option \c num_random_tests() (see StandardMemberCompositionMacros.h).
	STANDARD_MEMBER_COMPOSITION_MEMBERS( size_type, num_random_tests )
#ifdef DOXYGEN_COMPILE
		;
#endif		
	/// Members for option \c () warning_tol(see StandardMemberCompositionMacros.h).
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, warning_tol )
#ifdef DOXYGEN_COMPILE
		;
#endif		
	/// Members for option \c error_tol() (see StandardMemberCompositionMacros.h).
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, error_tol )
#ifdef DOXYGEN_COMPILE
		;
#endif		

	///
	/** Constructor (set default options).
	 *
	 * These default options are appropriate for even the largest vector spaces.
	 */
	VectorSpaceTester(
		bool         print_all_tests  = false
		,bool        print_vectors    = false
		,bool        throw_exception  = true
		,size_type   num_random_tests = 4
		,value_type  warning_tol      = 1e-14
		,value_type  error_tol        = 1e-10
		);

	///
	virtual ~VectorSpaceTester() {}

	///
	/** Run a vector space and the vectors it creates through a set of comprehensive tets.
	 *
	 * @param  space  [in] The vector space object to test.
	 * @param  out    [in/out] If <tt>out != NULL</tt> then output will be sent to <tt>*out</tt>.
	 *
	 * The behavior of this function greatly depends on a number of options (see \c VectorSpaceTester()
	 * for the default values for these options).  Access functions to set these options are provided
	 * by the prototypes of the macro <tt>STANDARD_MEMBER_COMPOSITION_MEMBERS()</tt>.
	 * <ul>
	 * <li> <b><tt>print_all_tests(bool)</tt></b>:  If <tt>print_all_tests() == true</tt>, then some output will be sent to
	 *      <tt>*out</tt> for every test performed.  This is useful to see all of tests that are performed and
	 *      in debugging.
	 * <li> <b><tt>print_vectors(bool)</tt></b>:  If <tt>print_vectors() == true</tt>, then all of the vectors will be printed
	 *      that are created durring the tests.  This option is really only needed durring initial debugging
	 *      and should only be used with small vector spaces since it will produce a lot of <tt>O(space.dim())</tt>
	 *      output.
	 * <li> <b><tt>throw_exception(bool)</tt></b>:  If <tt>throw_exception() == true</tt>, then any object that throws
	 *      an unexpected exception will cause that exception to be thrown clear of of this function.  If
	 *      <tt>out != NULL</tt> then the <tt>what()</tt> string will be printed to <tt>*out</tt> before the exception
	 *      is rethrown.  If <tt>throw_exception() == false</tt>, then all exceptions will be caught, printed to 
	 *      <tt>*out</tt> and then <tt>false</tt> is returned from the function.
	 * <li> <b><tt>num_random_tests(int)</tt></b>:  This is the number of random tests to perform per category of test.
	 *      A higher number will result is better validation but will consume more CPU time.
	 * <li> <b><tt>warning_tol(value_type)</tt></b>:  Any test with a relative error greater than <tt>warning_tol()</tt> will
	 *      result in a warning message printed to <tt>*out</tt>.
	 * <li> <b><tt>error_tol(value_type)</tt></b>:  Any test with a relative error greater than <tt>erfor_tol()</tt> will
	 *      result in an error message printed to <tt>*out</tt> and the function will immediatly return <tt>false</tt>.
	 * </ul>
	 */
	virtual bool check_vector_space(
		const VectorSpace &space
		,std::ostream     *out
		) const;

private:

	///
	void check_test(value_type err, std::ostream* out, bool* success) const;

}; // end class VectorSpaceTester

} // end namespace AbstractLinAlgPack

#endif // VECTOR_SPACE_TESTER_H