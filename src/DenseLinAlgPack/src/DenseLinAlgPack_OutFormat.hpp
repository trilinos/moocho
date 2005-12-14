// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef LINALGPACK_OUT_FORMAT_H
#define LINALGPACK_OUT_FORMAT_H

#include "DenseLinAlgPack_IOFormat.hpp"

namespace DenseLinAlgPack {

///
/* * Output stream operator function for const_bound_format objects.
  *
  * This template function performs the following tasks.
  * \begin{enumeration}
  * <li> Saves the formating state of #os#
  * <li> Sets the formating state of #os# to that stored in #bf.f()#
  * <li> Calls #output(os, bf.obj(), bf.extra_flags().flags())#
  * <li> Resets the streams formating state to its original.
  * \end{enumeration}
  *
  * The original formating state of #os# is preserved even if an exception
  * is thrown. 
  */
template<class T>
std::ostream& operator<<(std::ostream& os, const LinAlgPackIO::const_bound_format<T>& bf)
{
	using LinAlgPackIO::ios_format_memento;
	ios_format_memento old_format = ios_format_memento::save_format(os);
	try {
		bf.f().set_format(os);
		output( os, bf.obj(), bf.f().extra_flags().flags() );
	}
	catch(...) {
		old_format.set_format(os);
		throw;
	}
	old_format.set_format(os);
	return os;
}

/// Force a type conversion from #bound_format<T># to #const_bound_format<T># to call #operator<<#().
template<class T>
inline std::ostream& operator<<(std::ostream& os, const LinAlgPackIO::bound_format<T>& bf) {
	return operator<<( os, LinAlgPackIO::const_bound_format<T>( bf.f(), bf.obj() ) );
}

}	// end namespace DenseLinAlgPack

#endif // LINALGPACK_OUT_FORMAT_H
