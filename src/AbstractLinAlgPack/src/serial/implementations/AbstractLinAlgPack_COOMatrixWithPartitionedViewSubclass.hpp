// //////////////////////////////////////////////////////////////////////////////////
// COOMatrixWithPartitionedViewSubclass.hpp
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

#ifndef COO_MATRIX_WITH_PARTITIONED_VIEW_SUBCLASS_H
#define COO_MATRIX_WITH_PARTITIONED_VIEW_SUBCLASS_H

#include "MatrixWithOpConcreteEncap.hpp"
#include "COOMatrixWithPartitionedView.hpp"
#include "ConvertToSparseCompressedColumn.hpp"
#include "MatrixConvertToSparseFortranCompatible.hpp"
#include "COOMPartitionOut.hpp"

namespace SparseLinAlgPack {

///
/** Implementation of MatrixWithOp abstract interface for COOMatrixWithPartitionedView.
  *
  * Warning:  The return values of rows() and cols() may change if when the partitioned
  * view is setup it does not include the entire sparse matrix but this is unlikely.
  */
class COOMatrixWithPartitionedViewSubclass
	: public MatrixWithOpConcreteEncap<COOMatrixWithPartitionedView>
		, public ConvertToSparseCompressedColumn
		, public MatrixConvertToSparseFortranCompatible
{
public:

	///
	COOMatrixWithPartitionedViewSubclass()
	{}

	///
	COOMatrixWithPartitionedViewSubclass(const COOMatrixWithPartitionedView& m)
		: MatrixWithOpConcreteEncap<COOMatrixWithPartitionedView>(m)
	{}

	// /////////////////////////////////////////////////////
	// Overridden from Matrix
	size_type nz() const;

	// /////////////////////////////////////////////////////
	// Overridden from MatrixWithOp

	///
	std::ostream& output(std::ostream& out) const;

	// /////////////////////////////////////////////////////
	/** @name Level-1 BLAS */
	//@{

	/// (1) gms_lhs += alpha * op(M_rhs) (BLAS xAXPY)
	void Mp_StM(GenMatrixSlice* gms_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs) const;

	//		end Level-1 BLAS
	//@}

	// ////////////////////////////////////////////////////
	/** @name Level-2 BLAS */
	//@{

	/// (2) vs_lhs = alpha * op(M_rhs1) * vs_rhs2 + beta * vs_lhs (BLAS xGEMV)
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;

	/// (3) vs_lhs = alpha * op(M_rhs1) * sv_rhs2 + beta * vs_lhs (BLAS xGEMV)
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const SpVectorSlice& sv_rhs2, value_type beta) const;

	/// (4) result = vs_rhs1' * op(M_rhs2) * vs_rhs3
	value_type transVtMtV(const VectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const VectorSlice& vs_rhs3) const;

	/// (5) result = sv_rhs1' * op(M_rhs2) * sv_rhs3
	value_type transVtMtV(const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
		, const SpVectorSlice& sv_rhs3) const;

	//		end Level-2 BLAS
	//@}

	// ////////////////////////////////////////////////////
	/** @name Level-3 BLAS */
	//@{

	/// (6) gms_lhs = alpha * op(M_rhs1) * op(gms_rhs2) + beta * gms_lhs (right) (xGEMM)
	void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha
		, BLAS_Cpp::Transp trans_rhs1, const GenMatrixSlice& gms_rhs2
		, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

	/// (7) gms_lhs = alpha * op(gms_rhs1) * op(M_rhs2) + beta * gms_lhs (left) (xGEMM)
	void Mp_StMtM(GenMatrixSlice* gms_lhs, value_type alpha, const GenMatrixSlice& gms_rhs1
		, BLAS_Cpp::Transp trans_rhs1, BLAS_Cpp::Transp trans_rhs2, value_type beta) const;

	//		end Level-3 BLAS
	//@}

	// /////////////////////////////////////////////////////
	// Overridden from ConvertToSparseCompressedColumn

	///
	size_type num_in_column(
		  BLAS_Cpp::Transp					trans
		, size_type							col_offset
		, const IVector::value_type*		col_perm
		, size_type*						num_in_col	) const;
		
	///
	void insert_nonzeros(
		  BLAS_Cpp::Transp					trans
		, value_type						alpha
		, size_type							row_offset
		, size_type							col_offset
		, const IVector::value_type*		row_perm
		, const IVector::value_type*		col_perm
		, size_type*						next_nz_in_col
		, FortranTypes::f_dbl_prec*			D_val
		, FortranTypes::f_int*				D_row_i			) const;

	///
	value_type insert_scaled_nonzeros(
		  BLAS_Cpp::Transp					trans
		, value_type						scaled_max_ele
		, size_type							row_offset
		, size_type							col_offset
		, const IVector::value_type*		row_perm
		, const IVector::value_type*		col_perm
		, size_type*						next_nz_in_col
		, FortranTypes::f_dbl_prec*			D_val
		, FortranTypes::f_int*				D_row_i			) const;

	// ////////////////////////////////////////////////////////////
	// Overridden from MatrixConvertToSparseFortranCompatible

	///
	FortranTypes::f_int num_nonzeros( EExtractRegion extract_region ) const;

	///
	void coor_extract_nonzeros(
		  EExtractRegion extract_region
		, const FortranTypes::f_int len_Aval
			, FortranTypes::f_dbl_prec Aval[]
		, const FortranTypes::f_int len_Aij
			, FortranTypes::f_int Arow[]
			, FortranTypes::f_int Acol[]
			, const FortranTypes::f_int row_offset
			, const FortranTypes::f_int col_offset
		 ) const;

};	// end class COOMatrixWithPartitionedViewSubclass

}	// end namespace SparseLinAlgPack 

#endif	// COO_MATRIX_WITH_PARTITIONED_VIEW_SUBCLASS_H
