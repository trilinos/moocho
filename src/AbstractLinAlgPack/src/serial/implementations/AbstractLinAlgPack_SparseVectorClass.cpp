// /////////////////////////////////
// SparseVectorClass.cpp
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

#include "SparseLinAlgPack/include/SparseVectorClassDecl.h"

void SparseLinAlgPack::SparseVectorUtilityPack::assert_is_sorted(bool is_sorted)
{
	if(!is_sorted)
		throw NotSortedException("SparseVector***<> : The sparse vector is not sorted.");

}
