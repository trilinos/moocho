// ///////////////////////////////////////////////////////////
// random_vector.cpp
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

#include <stdlib.h>

#include <stdexcept>

#include "LinAlgPack/include/random_vector.h"
#include "LinAlgPack/include/VectorClass.h"

void LinAlgPack::seed_random_vector_generator( unsigned int s )
{
	srand(s);
}

void LinAlgPack::random_vector( value_type l, value_type u, VectorSlice* v )
{
	if(!v)
		throw std::invalid_argument( "random_vector(...) : Error, "
			"v can not be NULL" );
	if( l > u )
		throw std::invalid_argument( "random_vector(...) : Error, "
			"l can not be greater than u" );
	for( VectorSlice::iterator itr = v->begin(); itr != v->end(); )
		*itr++ = l + (double(rand())/RAND_MAX) * (u -l);
}
