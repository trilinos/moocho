// ////////////////////////////////////////////////////////////////////////////////////////////////
// LinAlgPackOptions.hpp
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
//
// Options for LinAlgPack compilation
//

#ifndef LINALGPACK_OPTIONS_H
#define LINALGPACK_OPTIONS_H

#include "LinAlgPackDebugAcronyms.hpp"
#include "extended_value_type.hpp"
#include "fortran_types.hpp"

#if !defined(LINALGPACK_NO_CHECKS)

/** @name {\bf LinAlgPack Options}.
  *
  * The header file LinAlgPackOptions.hpp contains the defines for several macros that
  * determine how the library is built.  The user should comment out any
  * macros that her or she does not want to be defined.  The definition of
  * these macros cause the library code to assert the preconditions documented
  * for each of the member and non-member functions and throw the listed exceptions
  * if they are not satisfied.  Precondtions are supposed to be the 
  * responcibility of the client code so the user may only want to define
  * these macros during debugging for better program verification.
  * If the user checks all of the preconditions listed in this documentation for the calls
  * to all functions then the checks performed by the library are redundant.
  */
//@{

///
/** If defined the library code checks to see if subscripts are in bounds for element access
  * an subregion indexing.  If the preconditions for the subscripting operations are
  * not satisfied then the listed exceptions will be thrown.
  */
#ifndef LINALGPACK_CHECK_RANGE
#define LINALGPACK_CHECK_RANGE 1
#endif

///
/** If defined the library code checks to see if the sizes of rhs arguments in expressions are compatible.
  * The exception std::length_error will be thrown if rhs sizes are not compatible.
  */
#ifndef LINALGPACK_CHECK_RHS_SIZES
#define LINALGPACK_CHECK_RHS_SIZES 1
#endif

///
/** If defined the library code checks to see if VectorSlice and GenMatrixSlice objects have valid constructions.
  * If they do not have valid constructions then an exception will be thrown.  The operation of these
  * checks may depend on the definition of the macro \Ref{LINALGPACK_CHECK_RANGE}.
  */
#ifndef LINALGPACK_CHECK_SLICE_SETUP
#define LINALGPACK_CHECK_SLICE_SETUP 1
#endif

#endif

namespace LinAlgPack{
/// Typedef for the value type of elements that is used for the library.
typedef FortranTypes::f_dbl_prec		value_type;
/// Typedef for the index type of elements that are used by the library
typedef FortranTypes::f_int				index_type;
/// Typedef for the size type of elements that are used by the library
typedef	FortranTypes::f_int				size_type;

}

//@}

#endif // LINALGPACK_OPTIONS_H
