// /////////////////////////////////////////////////////////////////////////////
// PermIn.hpp
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

#ifndef PERM_IN_H
#define PERM_IN_H

#include <istream>

#include "LinAlgPackTypes.hpp"

namespace LinAlgPack {

///
/** Input stream operator for IVector used as a permutation array.
  *
  * The input format is:\\
  * #size#\\
  * #i1 i2 i3 ... isize#\\
  */
std::istream& operator>>(std::istream& istrm, IVector& perm);

}	// end namespace LinAlgPack

#endif // PERM_IN_H
