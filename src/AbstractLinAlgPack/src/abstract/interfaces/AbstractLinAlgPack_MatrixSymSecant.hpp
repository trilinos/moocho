// //////////////////////////////////////////////////////////////////////////////////
// MatrixSymSecantUpdateable.hpp
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

#ifndef MATRIX_SYM_SECANT_UPDATEABLE_H
#define MATRIX_SYM_SECANT_UPDATEABLE_H

#include <stdexcept>

#include "MatrixSymInitDiagonal.hpp"

namespace AbstractLinAlgPack {

///
/** Mix-in interface for all polymorphic symmetric matrices that support secant updating.
 *
 * This interface is ment to be incorrporated in with a concrete <tt>AbstractLinAlgPack::MatrixWithOp</tt>
 * object that can implement some secant updating method.  Note that this is purely abstract interface\
 * and can be used in any application.
 *
 * Note that the methods <tt>AbstractLinAlgPack::MatrixSymInitDiagonal::init_identity()</tt> and
 * <tt>AbstractLinAlgPack::MatrixSymInitDiagonal::init_diagonal()</tt> do not state any postconditions
 * on the state of \c this after they are performed.  Subclasses of this interface also do not have
 * adhere to the obvious strick postconditions that these methods suggest but they should do the
 * "right thing" for the application.  If the client needs the obvious strict postconditions
 * for correct behavior, then the client is wise to test to see if \c this really is the
 * identity matrix or is a diagonal matrix (these can be cheap tests).
 */
class MatrixSymSecantUpdateable
	: virtual public AbstractLinAlgPack::MatrixSymInitDiagonal // doxygen needs the full name
{
public:

	///
	class UpdateSkippedException : public std::runtime_error
	{public: UpdateSkippedException(const std::string& what_arg) : std::runtime_error(what_arg) {}};

	///
	/** Perform a secant update of the matrix.
	 *
	 * The update is a secant update:
	 *
	 * <tt>B_hat * s = y</tt>
	 *
	 * It is assumed that <tt>||B_hat - B||</tt> (here \c B is the matrix before the
	 * update and \c B_hat is the matrix after the update) is not too large.
	 *
	 * The update vectors \c s and \c y may be used for workspace and are therefore
	 * not gaurented to be preserved.
	 *
	 * The vector \c Bs may be set by the client to <tt>B*s</tt>.  This may help the
	 * implementing subclass from having to compute it also.  Again, \c Bs may
	 * be used as workspace by subclass so it may change.
	 *
	 * If the update is not performed then an "UpdateSkippedException" will be
	 * thrown with an informative error message enclosed.
	 *
	 * Subclasses may also throw other unspecified exceptions but they should all
	 * be derived from <tt>std::exception</tt>.
	 * 
	 * @param  s   [in/work] On input must be set to the \c s vector (see above).
	 *             May be used as workspace.
	 * @param  y   [in/work] On input must be set to the \c y vector (see above)
	 *             May be used as workspace.
	 * @param  Bs  [in/work] On input (if not \c NULL), \c Bs may be set to the
	 *             <tt>B*s</tt> vector (see above).  If <tt>Bs == NULL</tt> on
	 *             input, then the subclass implementation will do without.
	 *             May be used as workspace.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>s != NULL</tt> (throw <tt>???</tt>)
	 * <li> \c s must be compatible with \c space_rows() of the underlying matrix.(throw <tt>???</tt>)
	 * <li> <tt>y != NULL</tt> (throw <tt>???</tt>)
	 * <li> \c y must be compatible with \c space_cols() of the underlying matrix.(throw <tt>???</tt>)
	 * <li> [<tt>Bs != NULL</tt>] \c Bs must be compatible with \c space_cols() of the underlying matrix.(throw <tt>???</tt>)
	 * </ul>
	 *
	 * Postconidtons:<ul>
	 * <li> <tt>(*this) * x \approx y</tt>.  In other words, we expect the secant condition will hold.
	 *      In almost every application, this is required for correct behavior so clients should
	 *      meet this condition.
	 * </ul>
	 */
	virtual void secant_update(
		VectorWithOpMutable     *s
		,VectorWithOpMutable    *y
		,VectorWithOpMutable    *Bs = NULL
		) = 0;
	
};	// end class MatrixSymSecantUpdateable 

}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_SYM_SECANT_UPDATEABLE_H
