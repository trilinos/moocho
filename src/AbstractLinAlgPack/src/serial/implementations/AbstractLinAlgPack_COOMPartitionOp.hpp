// //////////////////////////////////////////////////////////////////////////////////
// COOMPartitionOp.h

#ifndef COOM_PARTITION_OP_H
#define COOM_PARTITION_OP_H

#include "COOMatrixTmplOp.h"

namespace SparseLinAlgPack {

// ///////////////////////////////////////////////////////////////////////////
/** @name {\bf Linear algebra operations for Partition<> (Level 2,3 BLAS)}.
  *
  * These are inline functions that call the linear algebra functions templated
  * for the COOMatrixTemplateInterface.  Here the #COOM# is replaced with
  * the more standard #M# for a more common usage as with LinAlgPack.
  */
//@{

// ////////////////////////////////////////////////////////////////////
/** @name {\bf Level-1 BLAS matrix-matrix element-wise operations}.
  */
//@{

/// gms_lhs += alpha * op(coom_rhs) (time = O(coom_rhs.nz()), space = O(1)).
template <class T_Indice, class T_Value>
inline void Mp_StM( GenMatrixSlice* gms_lhs, value_type alpha
	, const COOMatrixPartitionedViewUtilityPack::Partition<T_Indice,T_Value>& coom_rhs
	, BLAS_Cpp::Transp trans_rhs )
{
	Mp_StCOOM(gms_lhs, alpha, coom_rhs, trans_rhs);
}

//@}

// ///////////////////////////////////////////////////////////////////////
/** @name {\bf Level-2 BLAS (vector-matrtix) Liner Algebra Operations}.
  *
  * These are functions that implement linear algebra operations that have Vector or
  * VectorSlice objects as lhs operands and a mixture of matrix and VectorSlice 
  * rhs operands.
  */
//@{

/// vs_lhs += alpha * op(coom_rhs1) * vs_rhs2 (BLAS xGEMV) (time = O(coom_rhs.nz()), space = O(1)).
template <class T_Indice, class T_Value>
inline void Vp_StMtV( VectorSlice* vs_lhs, value_type alpha
	, const COOMatrixPartitionedViewUtilityPack::Partition<T_Indice,T_Value>& coom_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2 )
{
	Vp_StCOOMtV(vs_lhs, alpha, coom_rhs1, trans_rhs1, vs_rhs2);
}

/// vs_lhs += alpha * op(coom_rhs1) * vs_rhs2 (BLAS xGEMV) (time = O(coom_rhs.nz()), space = O(1)).
template <class T_Indice, class T_Value>
inline void Vp_StMtV( VectorSlice* vs_lhs, value_type alpha
	, const COOMatrixPartitionedViewUtilityPack::Partition<T_Indice,T_Value>& coom_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2, value_type beta )
{
	if( beta == 0.0 ) {
		*vs_lhs = 0.0;	// We must specifically assign this in case uninitialized memory
						// was used and you might get NaN by accident (this happened to me!).
	}
	else {
		LinAlgPack::Vt_S(vs_lhs,beta);
	}
	Vp_StCOOMtV(vs_lhs, alpha, coom_rhs1, trans_rhs1, vs_rhs2);
}

//@}

// /////////////////////////////////////////////////////////////////////////////////
/** @name {\bf Level-3 BLAS (matrix-matrix) Linear Algebra Operations}.
  */
//@{

/// gms_lhs += alpha * op(coom_rhs1) * op(gms_rhs2) (right) (BLAS xGEMM).
template <class T_Indice, class T_Value>
inline void Mp_StMtM( GenMatrixSlice* gms_lhs, value_type alpha
	, const COOMatrixPartitionedViewUtilityPack::Partition<T_Indice,T_Value>& coom_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2, BLAS_Cpp::Transp trans_rhs2 )
{
	Mp_StCOOMtM(gms_lhs, alpha, coom_rhs1, trans_rhs1, gms_rhs2, trans_rhs2);
}

/// gms_lhs += alpha * op(gms_rhs1) * op(coom_rhs2) (left) (BLAS xGEMM).
template <class T_Indice, class T_Value>
inline void Mp_StMtM( GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1
	, const COOMatrixPartitionedViewUtilityPack::Partition<T_Indice,T_Value>& coom_rhs2
	, BLAS_Cpp::Transp trans_rhs2 )
{
	Mp_StMtCOOM(gms_lhs, alpha, gms_rhs1, trans_rhs1, coom_rhs2, trans_rhs2);
}

///
/** gms_lhs = alpha * op(coom_rhs1) * op(M2_rhs2) + beta * gms_lhs (right).
  *
  * Was part of LinAlgOpPack but had to be removed becasue of unexplained bug.
  */
template <class T_Indice, class T_Value, class M2>
inline
void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha
	, const COOMatrixPartitionedViewUtilityPack::Partition<T_Indice,T_Value>& coom_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2
	, value_type beta)
{
	Mp_MtM_assert_sizes(	  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
							, coom_rhs1.rows(), coom_rhs1.cols(), trans_rhs1
							, M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2 );
	if( beta == 0.0 ) {
		*gms_lhs = 0.0;	// We must specifically assign this in case uninitialized memory
						// was used and you might get NaN by accident (this happened to me!).
	}
	else {
		LinAlgPack::Mt_S(gms_lhs,beta);
	}
	Mp_StMtM(gms_lhs,alpha,coom_rhs1,trans_rhs1,M2_rhs2,trans_rhs2);
}

///
/** gms_lhs = alpha * op(M1_rhs1) * op(coom_rhs2) + beta * gms_lhs (right).
  *
  * Was part of LinAlgOpPack but had to be removed becasue of unexplained bug.
  */
template <class M1, class T_Indice, class T_Value>
inline
void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const M1& M1_rhs1
	, BLAS_Cpp::Transp trans_rhs1
	, const COOMatrixPartitionedViewUtilityPack::Partition<T_Indice,T_Value>& coom_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta)
{
	Mp_MtM_assert_sizes(	  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
							, M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1
							, coom_rhs2.rows(), coom_rhs2.cols(), trans_rhs2 );
	if( beta == 0.0 ) {
		*gms_lhs = 0.0;	// We must specifically assign this in case uninitialized memory
						// was used and you might get NaN by accident (this happened to me!).
	}
	else {
		LinAlgPack::Mt_S(gms_lhs,beta);
	}
	Mp_StMtM(gms_lhs,alpha,M1_rhs1,trans_rhs1,coom_rhs2,trans_rhs2);
}

//@}

//@}

// ///////////////////////////////////////////////////////////////////////////
/** @name {\bf Linear algebra operations for TransposedPartition<> (Level 2,3 BLAS)}.
  *
  * These are inline functions that call the linear algebra functions templated
  * for the COOMatrixTemplateInterface.  Here the #COOM# is replaced with
  * the more standard #M# for a more common usage as with LinAlgPack.
  */
//@{

// ////////////////////////////////////////////////////////////////////
/** @name {\bf Level-1 BLAS matrix-matrix element-wise operations}.
  */
//@{

/// gms_lhs += alpha * op(coom_rhs) (time = O(coom_rhs.nz()), space = O(1)).
template <class T_Indice, class T_Value>
inline void Mp_StM( GenMatrixSlice* gms_lhs, value_type alpha
	, const COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>& coom_rhs
	, BLAS_Cpp::Transp trans_rhs )
{
	Mp_StCOOM(gms_lhs, alpha, coom_rhs, trans_rhs);
}

//@}

// ///////////////////////////////////////////////////////////////////////
/** @name {\bf Level-2 BLAS (vector-matrtix) Liner Algebra Operations}.
  *
  * These are functions that implement linear algebra operations that have Vector or
  * VectorSlice objects as lhs operands and a mixture of matrix and VectorSlice 
  * rhs operands.
  */
//@{

/// vs_lhs += alpha * op(coom_rhs1) * vs_rhs2 (BLAS xGEMV) (time = O(coom_rhs.nz()), space = O(1)).
template <class T_Indice, class T_Value>
inline void Vp_StMtV( VectorSlice* vs_lhs, value_type alpha
	, const COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>& coom_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2 )
{
	Vp_StCOOMtV(vs_lhs, alpha, coom_rhs1, trans_rhs1, vs_rhs2);
}

/// vs_lhs += alpha * op(coom_rhs1) * vs_rhs2 (BLAS xGEMV) (time = O(coom_rhs.nz()), space = O(1)).
template <class T_Indice, class T_Value>
inline void Vp_StMtV( VectorSlice* vs_lhs, value_type alpha
	, const COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>& coom_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const VectorSlice& vs_rhs2, value_type beta )
{
	if( beta == 0.0 ) {
		*vs_lhs = 0.0;	// We must specifically assign this in case uninitialized memory
						// was used and you might get NaN by accident (this happened to me!).
	}
	else {
		LinAlgPack::Vt_S(vs_lhs,beta);
	}
	Vp_StCOOMtV(vs_lhs, alpha, coom_rhs1, trans_rhs1, vs_rhs2);
}

//@}

// /////////////////////////////////////////////////////////////////////////////////
/** @name {\bf Level-3 BLAS (matrix-matrix) Linear Algebra Operations}.
  */
//@{

/// gms_lhs += alpha * op(coom_rhs1) * op(gms_rhs2) (right) (BLAS xGEMM).
template <class T_Indice, class T_Value>
inline void Mp_StMtM( GenMatrixSlice* gms_lhs, value_type alpha
	, const COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>& coom_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2, BLAS_Cpp::Transp trans_rhs2 )
{
	Mp_StCOOMtM(gms_lhs, alpha, coom_rhs1, trans_rhs1, gms_rhs2, trans_rhs2);
}

/// gms_lhs += alpha * op(gms_rhs1) * op(coom_rhs2) (left) (BLAS xGEMM).
template <class T_Indice, class T_Value>
inline void Mp_StMtM( GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
	, BLAS_Cpp::Transp trans_rhs1
	, const COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>& coom_rhs2
	, BLAS_Cpp::Transp trans_rhs2 )
{
	Mp_StMtCOOM(gms_lhs, alpha, gms_rhs1, trans_rhs1, coom_rhs2, trans_rhs2);
}

///
/** gms_lhs = alpha * op(coom_rhs1) * op(M2_rhs2) + beta * gms_lhs (right).
  *
  * Was part of LinAlgOpPack but had to be removed becasue of unexplained bug.
  */
template <class T_Indice, class T_Value, class M2>
inline
void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha
	, const COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>& coom_rhs1
	, BLAS_Cpp::Transp trans_rhs1, const M2& M2_rhs2, BLAS_Cpp::Transp trans_rhs2
	, value_type beta)
{
	Mp_MtM_assert_sizes(	  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
							, coom_rhs1.rows(), coom_rhs1.cols(), trans_rhs1
							, M2_rhs2.rows(), M2_rhs2.cols(), trans_rhs2 );
	if( beta == 0.0 ) {
		*gms_lhs = 0.0;	// We must specifically assign this in case uninitialized memory
						// was used and you might get NaN by accident (this happened to me!).
	}
	else {
		LinAlgPack::Mt_S(gms_lhs,beta);
	}
	Mp_StMtM(gms_lhs,alpha,coom_rhs1,trans_rhs1,M2_rhs2,trans_rhs2);
}

///
/** gms_lhs = alpha * op(M1_rhs1) * op(coom_rhs2) + beta * gms_lhs (right).
  *
  * Was part of LinAlgOpPack but had to be removed becasue of unexplained bug.
  */
template <class M1, class T_Indice, class T_Value>
inline
void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const M1& M1_rhs1
	, BLAS_Cpp::Transp trans_rhs1
	, const COOMatrixPartitionedViewUtilityPack::TransposedPartition<T_Indice,T_Value>& coom_rhs2
	, BLAS_Cpp::Transp trans_rhs2, value_type beta)
{
	Mp_MtM_assert_sizes(	  gms_lhs->rows(), gms_lhs->cols(), BLAS_Cpp::no_trans
							, M1_rhs1.rows(), M1_rhs1.cols(), trans_rhs1
							, coom_rhs2.rows(), coom_rhs2.cols(), trans_rhs2 );
	if( beta == 0.0 ) {
		*gms_lhs = 0.0;	// We must specifically assign this in case uninitialized memory
						// was used and you might get NaN by accident (this happened to me!).
	}
	else {
		LinAlgPack::Mt_S(gms_lhs,beta);
	}
	Mp_StMtM(gms_lhs,alpha,M1_rhs1,trans_rhs1,coom_rhs2,trans_rhs2);
}

//@}

//@}

} // end namespace SparseLinAlgPack

#endif	// COOM_PARTITION_OP_H
