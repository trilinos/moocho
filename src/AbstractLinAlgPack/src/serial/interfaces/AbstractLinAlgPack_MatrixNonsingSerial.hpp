// ///////////////////////////////////////////////////////////////////
// MatrixNonsingularSerial.hpp
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

#ifndef SLAP_MATRIX_NONSINGULAR_SERIAL_H
#define SLAP_MATRIX_NONSINGULAR_SERIAL_H

#include "SparseLinAlgPackTypes.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/MatrixNonsing.hpp"

namespace SparseLinAlgPack {

///
/** Abstract base class for all <tt>AbstractLinAlgPack::MatrixNonsing</tt> objects
 * implemented in shared memory space.
 *
 * This base class does a mapping from fully abstract linear algebra to shared memory
 * linear algebra.
 *
 * These methods should not be called directly but instead should be called through
 * the line \ref MatrixNonsingularSerial_funcs "non-member functions" that are provided.
 */
class MatrixNonsingularSerial
	: virtual public AbstractLinAlgPack::MatrixNonsing // doxygen needs full name
{
public:

	/** @name Level-2 BLAS */
	//@{

	/// v_lhs	= inv(op(M_rhs1)) * vs_rhs2
	virtual void V_InvMtV(
		DVector* v_lhs, BLAS_Cpp::Transp trans_rhs1
		,const DVectorSlice& vs_rhs2) const;
	/// vs_lhs	= inv(op(M_rhs1)) * vs_rhs2
	virtual void V_InvMtV(
		DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
		,const DVectorSlice& vs_rhs2) const = 0;
	/// v_lhs	= inv(op(M_rhs1)) * sv_rhs2
	virtual void V_InvMtV(
		DVector* v_lhs, BLAS_Cpp::Transp trans_rhs1
		,const SpVectorSlice& sv_rhs2) const;
	/// vs_lhs	= inv(op(M_rhs1)) * sv_rhs2
	virtual void V_InvMtV(
		DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
		,const SpVectorSlice& sv_rhs2) const;
	/// result	= vs_rhs1' * inv(op(M_rhs2)) * vs_rhs3
	virtual value_type transVtInvMtV(
		const DVectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2, const DVectorSlice& vs_rhs3) const;
	/// result	= sv_rhs1' * inv(op(M_rhs2)) * sv_rhs3
	virtual value_type transVtInvMtV(
		const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const;

	//		end Level-2 BLAS
	//@}

	/** @name Level-3 BLAS */
	//@{

	/// gm_lhs	= alpha * inv(op(M_rhs1)) * op(gms_rhs2) (right)
	virtual void M_StInvMtM(
		DMatrix* gm_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const DMatrixSlice& gms_rhs2, BLAS_Cpp::Transp trans_rhs2 ) const;
	/// gms_lhs	= alpha * inv(op(M_rhs1)) * op(gms_rhs2) (right)
	virtual void M_StInvMtM(
		DMatrixSlice* gms_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const DMatrixSlice& gms_rhs2, BLAS_Cpp::Transp trans_rhs2 ) const;
	/// gm_lhs	= alpha * op(gms_rhs1) * inv(op(M_rhs2)) (left)
	virtual void M_StMtInvM(
		DMatrix* gm_lhs, value_type alpha
		,const DMatrixSlice& gms_rhs1, BLAS_Cpp::Transp trans_rhs1
		,BLAS_Cpp::Transp trans_rhs2 ) const;
	/// gms_lhs	= alpha * op(gms_rhs1) * inv(op(M_rhs2)) (left)
	virtual void M_StMtInvM(
		DMatrixSlice* gms_lhs, value_type alpha
		,const DMatrixSlice& gms_rhs1, BLAS_Cpp::Transp trans_rhs1
		,BLAS_Cpp::Transp trans_rhs2 ) const;
	/// gm_lhs	= alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right)
	virtual void M_StInvMtM(
		DMatrix* gm_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixWithOpSerial& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2 ) const;
	/// gms_lhs	= alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right)
	virtual void M_StInvMtM(
		DMatrixSlice* gms_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixWithOpSerial& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2 ) const;
	/// gm_lhs	= alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left)
	virtual void M_StMtInvM(
		DMatrix* gm_lhs, value_type alpha
		,const MatrixWithOpSerial& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		,BLAS_Cpp::Transp trans_rhs2 ) const;
	/// gms_lhs	= alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left)
	virtual void M_StMtInvM(
		DMatrixSlice* gms_lhs, value_type alpha
		,const MatrixWithOpSerial& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		,BLAS_Cpp::Transp trans_rhs2 ) const;

	//		end Level-3 BLAS
	//@}

	/** Overridden from MatrixNonsing */
	//@{

	/// v_lhs	= inv(op(M_rhs1)) * vs_rhs2
	void V_InvMtV(
		VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		,const Vector& v_rhs2) const;
	/// v_lhs	= inv(op(M_rhs1)) * sv_rhs2
	void V_InvMtV(
		VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		,const SpVectorSlice& sv_rhs2) const;
	/// result	= vs_rhs1' * inv(op(M_rhs2)) * vs_rhs3
	value_type transVtInvMtV(
		const Vector& v_rhs1
		,BLAS_Cpp::Transp trans_rhs2
		,const Vector& v_rhs3) const;
	/// m_lhs = alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right).
	void M_StInvMtM(
		MatrixOp* m_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
		) const;
	/// m_lhs = alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left).
	void M_StMtInvM(
		MatrixOp* m_lhs, value_type alpha
		,const MatrixOp& mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
		,BLAS_Cpp::Transp trans_rhs2
		) const;

	//@}

};	// end class MatrixNonsingularSerial

/** \defgroup MatrixNonsingularSerial_funcs MatrixNonsingularSerial nonmember inline functions.
 *
 * These nonmember functions allow operations to be called on \c MatrixNonsingularSerial objects
 * in similar manner to those in \c DenseLinAlgPack.
 */
//@{

/** @name Level-2 BLAS */
//@{

/// v_lhs	= inv(op(M_rhs1)) * vs_rhs2
inline void V_InvMtV(DVector* v_lhs, const MatrixNonsingularSerial& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2)
{
	M_rhs1.V_InvMtV(v_lhs,trans_rhs1,vs_rhs2);
}

/// vs_lhs	= inv(op(M_rhs1)) * vs_rhs2
inline void V_InvMtV(DVectorSlice* vs_lhs, const MatrixNonsingularSerial& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const DVectorSlice& vs_rhs2)
{
	M_rhs1.V_InvMtV(vs_lhs,trans_rhs1,vs_rhs2);
}

/// v_lhs	= inv(op(M_rhs1)) * sv_rhs2
inline void V_InvMtV(DVector* v_lhs, const MatrixNonsingularSerial& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2)
{
	M_rhs1.V_InvMtV(v_lhs,trans_rhs1,sv_rhs2);
}

/// vs_lhs	= inv(op(M_rhs1)) * sv_rhs2
inline void V_InvMtV(DVectorSlice* vs_lhs, const MatrixNonsingularSerial& M_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const SpVectorSlice& sv_rhs2)
{
	M_rhs1.V_InvMtV(vs_lhs,trans_rhs1,sv_rhs2);
}

/// result	= vs_rhs1' * inv(op(M_rhs2)) * vs_rhs3
inline value_type transVtInvMtV(const DVectorSlice& vs_rhs1, const MatrixNonsingularSerial& M_rhs2
	, BLAS_Cpp::Transp trans_rhs2, const DVectorSlice& sv_rhs3)
{
	return M_rhs2.transVtInvMtV(vs_rhs1,trans_rhs2,sv_rhs3);
}

/// result	= sv_rhs1' * inv(op(M_rhs2)) * sv_rhs3
inline value_type transVtInvMtV(const SpVectorSlice& sv_rhs1, const MatrixNonsingularSerial& M_rhs2
	, BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3)
{
	return M_rhs2.transVtInvMtV(sv_rhs1,trans_rhs2,sv_rhs3);
}

//		end Level-2 BLAS
//@}

/** @name Level-3 BLAS */
//@{

/// gm_lhs	= alpha * inv(op(M_rhs1)) * op(gms_rhs2) (right)
inline void M_StInvMtM(
	DMatrix* gm_lhs, value_type alpha
	,const MatrixNonsingularSerial& M_rhs1,   BLAS_Cpp::Transp trans_rhs1
	,const DMatrixSlice&          gms_rhs2, BLAS_Cpp::Transp trans_rhs2
	)
{
	M_rhs1.M_StInvMtM(gm_lhs,alpha,trans_rhs1,gms_rhs2,trans_rhs2);
}

/// gms_lhs	= alpha * inv(op(M_rhs1)) * op(gms_rhs2) (right)
inline void M_StInvMtM(
	DMatrixSlice* gms_lhs, value_type alpha
	,const MatrixNonsingularSerial& M_rhs1,   BLAS_Cpp::Transp trans_rhs1
	,const DMatrixSlice&          gms_rhs2, BLAS_Cpp::Transp trans_rhs2
	)
{
	M_rhs1.M_StInvMtM(gms_lhs,alpha,trans_rhs1,gms_rhs2,trans_rhs2);
}

/// gm_lhs	= alpha * op(gms_rhs1) * inv(op(M_rhs2)) (left)
inline void M_StMtInvM(
	DMatrix* gm_lhs, value_type alpha
	,const DMatrixSlice&          gms_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixNonsingularSerial& M_rhs2,   BLAS_Cpp::Transp trans_rhs2
	)
{
	M_rhs2.M_StMtInvM(gm_lhs,alpha,gms_rhs1,trans_rhs1,trans_rhs2);
}

/// gms_lhs	= alpha * op(gms_rhs1) * inv(op(M_rhs2)) (left)
inline void M_StMtInvM(
	DMatrixSlice* gms_lhs, value_type alpha
	,const DMatrixSlice&          gms_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixNonsingularSerial& M_rhs2,   BLAS_Cpp::Transp trans_rhs2
	)
{
	M_rhs2.M_StMtInvM(gms_lhs,alpha,gms_rhs1,trans_rhs1,trans_rhs2);
}

/// gm_lhs	= alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right)
inline void M_StInvMtM(
	DMatrix* gm_lhs, value_type alpha
	,const MatrixNonsingularSerial& M_rhs1,   BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOpSerial&      mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	)
{
	M_rhs1.M_StInvMtM(gm_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2);
}

/// gms_lhs	= alpha * inv(op(M_rhs1)) * op(mwo_rhs2) (right)
inline void M_StInvMtM(
	DMatrixSlice* gms_lhs, value_type alpha
	,const MatrixNonsingularSerial& M_rhs1,   BLAS_Cpp::Transp trans_rhs1
	,const MatrixWithOpSerial&      mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
	)
{
	M_rhs1.M_StInvMtM(gms_lhs,alpha,trans_rhs1,mwo_rhs2,trans_rhs2);
}

/// gm_lhs	= alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left)
inline void M_StMtInvM(
	DMatrix* gm_lhs, value_type alpha
	,const MatrixWithOpSerial&      mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixNonsingularSerial& M_rhs2,   BLAS_Cpp::Transp trans_rhs2
	)
{
	M_rhs2.M_StMtInvM(gm_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2);
}

/// gms_lhs	= alpha * op(mwo_rhs1) * inv(op(M_rhs2)) (left)
inline void M_StMtInvM(
	DMatrixSlice* gms_lhs, value_type alpha
	,const MatrixWithOpSerial&      mwo_rhs1, BLAS_Cpp::Transp trans_rhs1
	,const MatrixNonsingularSerial& M_rhs2,   BLAS_Cpp::Transp trans_rhs2
	)
{
	M_rhs2.M_StMtInvM(gms_lhs,alpha,mwo_rhs1,trans_rhs1,trans_rhs2);
}

//		end Level-3 BLAS
//@}

//		end Inline non-member operation functions
//@}

}	// end namespace SparseLinAlgPack

#endif	// SLAP_MATRIX_NONSINGULAR_SERIAL_H
