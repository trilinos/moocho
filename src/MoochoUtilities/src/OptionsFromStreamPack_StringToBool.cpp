// //////////////////////////////////////////////////////////////////////
// StringToBool.cpp
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

#include <string.h>
#include <sstream>

#include "StringToBool.hpp"
#include "OptionsFromStreamExceptions.hpp"

bool OptionsFromStreamPack ::StringToBool( const char* opt_name, const char* str ) {

	if( !::strcmp( str, "true" ) )
		return true; 
	if( !::strcmp( str, "false" ) )
		return false;
	std::ostringstream omsg;
	omsg	<< "StringToBool(...): "
			<< "Error, the string \"" << str << "\" must be \"true\" or \"false\" for \""
			<< opt_name << "\"";
	throw InputException( omsg.str() );

	return false;	// will never execute.
}
