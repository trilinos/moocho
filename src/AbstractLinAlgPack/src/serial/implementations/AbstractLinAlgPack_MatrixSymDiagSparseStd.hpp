// ///////////////////////////////////////////////////////////////////
// MatrixSymDiagonalSparseStd.hpp
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

#ifndef SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_STD_H
#define SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_STD_H

#include "MatrixSymDiagonalSparse.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/SpVectorClass.hpp"

namespace SparseLinAlgPack {

///
/** Concrete subclass for a serial symmetric diagonal matrix with many zeros on the diagonal.
  *
  * The underlying diagonal vector is sorted and determines the dimensions of the
  * matrix.
  *
  * The default constructor, copy constructor are allowed.
  */
class MatrixSymDiagonalSparseStd: virtual public MatrixSymDiagonalSparse {
public:

	/** @name Constructors/initializes */
	//@{

	/// Construct uninitialized
	MatrixSymDiagonalSparseStd()
	{}

	///
	/** Construct the diagonal.
	  */
	MatrixSymDiagonalSparseStd( const SpVectorSlice& diag );

	///
	/** Reinitialize the diagonal.
	  */
	void initialize( const SpVectorSlice& diag );

	//@}

	/** @name Overridden from MatrixOp */
	//@{

	///
	MatrixOp& operator=(const MatrixOp& m);

	//@}

	/** Overridden from MatrixDiagonalSparse */
	//@{

	///
	const SpVectorSlice diag() const;

	//@}

private:
	
	SpVector	diag_;

};	// end class MatrixDiagonalSparse

}	// end namespace SparseLinAlgPack

#endif	// SPARSE_LINALG_PACK_MATRIX_DIAGONAL_SPARSE_STD_H
