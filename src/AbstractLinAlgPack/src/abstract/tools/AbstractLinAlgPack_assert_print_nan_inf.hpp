// ///////////////////////////////////////////////////////////
// assert_print_nan_inf.hpp
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

#ifndef ASSERT_PRINT_NAN_INF_H
#define ASSERT_PRINT_NAN_INF_H

#include <stdexcept>

#include "AbstractLinAlgPackTypes.hpp"

namespace AbstractLinAlgPack {

///
class NaNInfException : public std::runtime_error
{public: NaNInfException(const std::string& what_arg) : std::runtime_error(what_arg) {}};

///
/** This function asserts if a value_type scalare is a NaN or Inf and optionally
 * prints out these entires.
 * 
 * @param  val             [in] Value the check
 * @param  name            [in] Name of the scale variable for output purposes
 * @param  throw_excpt     [in] If true and is found to be a NaN or Inf
 *                         then a NaNInfException excetion is thrown after
 *                         any output.
 * @param  out             [in/out] If out==NULL then not output is produced.
 *                         If out!=NULL and val is not
 *                         NaN or Inf, then no output is produced.
 *                         If out!=NULL and val is
 *                         NaN or Inf then this will be printed before any
 *                         execption is thrown.
 *									
 * @return Returns true if val is not NaN or Inf.  If val
 * is NaN or Inf then false will be returned unless a
 * excetion NaNInfException was thrown (throw_except==true).
 */
bool assert_print_nan_inf( const value_type& val, char name[]
	, bool throw_excpt, std::ostream* out );

///
/** This function asserts if a vector has any NaN or inf entries and optionally
 * prints out these entires.
 * 
 * @param  v              [in]	DVector slice to check
 * @param  name           [in]	Name of the vector for output purposes
 * @param  throw_excpt    [in]	If true and an entry is found to be a NaN or Inf
 *                        then a NaNInfException excetion is thrown after
 *                        any output.
 * @param  out            [in/out]	If out==NULL then not output is produced.
 *                        If out!=NULL and none of the entries is
 *                        NaN or Inf, then no output is produced.
 *                        If out!=NULL then any entries that are
 *                        NaN or Inf will be printed before any
 *                        execption is thrown.
 *
 * @return Returns true none of the entries are NaN or Inf.  If one of the
 * entries is NaN or Inf then false will be returned unless an
 * excetion was thrown (throw_except==true).
 */
bool assert_print_nan_inf( const Vector& v, char name[]
	, bool throw_excpt, std::ostream* out );

}	// end namespace AbstractLinAlgPack

#endif // ASSERT_PRINT_NAN_INF_H
