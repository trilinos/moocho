// //////////////////////////////////////////////////////////////////////////////////
// MatrixBase.cpp
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

#include "AbstractLinAlgPack/include/MatrixBase.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"

namespace AbstractLinAlgPack {

size_type MatrixBase::rows() const
{
	return this->space_cols().dim();
}

size_type MatrixBase::cols() const
{
	return this->space_rows().dim();
}

size_type MatrixBase::nz() const
{
	return rows() * cols();
}

}	// end namespace AbstractLinAlgPack