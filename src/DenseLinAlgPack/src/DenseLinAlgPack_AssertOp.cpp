// //////////////////////////////////////////////////////////////////////////////////
// LinAlgPackAssertOp.cpp

#include <stdexcept>
#include <string>

#include "../include/LinAlgPackAssertOp.h"

#ifdef LINALGPACK_CHECK_RHS_SIZES

void LinAlgPack::Vp_V_assert_sizes(size_type v_lhs_size, size_type v_rhs_size)
{
	if(v_lhs_size != v_rhs_size)
		throw std::length_error("Vp_V_assert_sizes(...) : The sizes of v_lhs and v_rhs "
			"in the operation v_lhs += op v_rhs do not match");
}

void LinAlgPack::VopV_assert_sizes(size_type v_rhs1_size, size_type v_rhs2_size)
{
	if(v_rhs1_size != v_rhs2_size)
		throw std::length_error("VopV_assert_sizes(...) : The sizes of v_rhs1 and v_rhs2 "
			"in the operation v_rhs1 op v_rhs2 do not match");
}

void LinAlgPack::Mp_M_assert_sizes(size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
	, size_type m_rhs_rows, size_type m_rhs_cols, BLAS_Cpp::Transp trans_rhs)
{
	if(		rows(m_lhs_rows,m_lhs_cols,trans_lhs) != rows(m_rhs_rows,m_rhs_cols,trans_rhs)
		||	cols(m_lhs_rows,m_lhs_cols,trans_lhs) != cols(m_rhs_rows,m_rhs_cols,trans_rhs) )
	{
		throw std::length_error("Mp_M_assert_sizes(...) : The sizes of m_lhs and m_rhs "
			"in the operation op(m_lhs) += op op(m_rhs) do not match");
	}
}

void LinAlgPack::MopM_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
	, size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{
	if(		rows(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != rows(m_rhs2_rows,m_rhs2_cols,trans_rhs2)
		||	cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != cols(m_rhs2_rows,m_rhs2_cols,trans_rhs2) )
	{
		throw std::length_error("Mp_M_assert_sizes(...) : The sizes of m_rhs1 and m_rhs2 "
			"in the operation op(m_rhs1) op op(m_rhs2) do not match");
	}
}

void LinAlgPack::MtV_assert_sizes(size_type m_rhs1_rows, size_type m_rhs1_cols
	, BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size)
{
	if(cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != v_rhs2_size)
		throw std::length_error("MtV_assert_sizes(...) : The number of columns in "
			"m_rhs1 and the size of v_rhs2 in the operation v_lhs += op(m_rhs1) * v_rhs2 "
			"do not match");
}

void LinAlgPack::Vp_MtV_assert_sizes(size_type v_lhs_size, size_type m_rhs1_rows
	, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1, size_type v_rhs2_size)
{
	if(cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != v_rhs2_size)
		throw std::length_error("Vp_MtV_assert_sizes(...) : The number of columns in"
			" m_rhs1 and the size of v_rhs2 in the operation v_lhs += op(m_rhs1) * v_rhs2"
			" do not match");
	if(rows(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != v_lhs_size)
		throw std::length_error("Vp_MtV_assert_sizes(...) : The number of rows in"
			" m_rhs1 and the size of v_lhs in the operation v_lhs += op(m_rhs1) * v_rhs2"
			" do not match");
}

void LinAlgPack::MtM_assert_sizes(
	  size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
	, size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{
	if(cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != rows(m_rhs2_rows,m_rhs2_cols,trans_rhs2))
		throw std::length_error("MtM_assert_sizes(...) : The number of columns in"
			" m_rhs1 and the number of rows in m_rhs2 in the operation"
			" op(m_lhs) += op(m_rhs1) * op(m_rhs2) do not match");
}

void LinAlgPack::Mp_MtM_assert_sizes(
	  size_type m_lhs_rows, size_type m_lhs_cols, BLAS_Cpp::Transp trans_lhs
	, size_type m_rhs1_rows, size_type m_rhs1_cols, BLAS_Cpp::Transp trans_rhs1
	, size_type m_rhs2_rows, size_type m_rhs2_cols, BLAS_Cpp::Transp trans_rhs2)
{
	if(cols(m_rhs1_rows,m_rhs1_cols,trans_rhs1) != rows(m_rhs2_rows,m_rhs2_cols,trans_rhs2))
		throw std::length_error("Mp_MtM_assert_sizes(...) : The number of columns in"
			" m_rhs1 and the number of rows in m_rhs2 in the operation"
			" op(m_lhs) += op(m_rhs1) * op(m_rhs2) do not match");
	if(rows(m_lhs_rows,m_lhs_cols,trans_lhs) != rows(m_rhs1_rows,m_rhs1_cols,trans_rhs1))
		throw std::length_error("Mp_MtM_assert_sizes(...) : The number of rows in"
			" m_lhs and the number of rows in m_rhs1 in the operation"
			" op(m_lhs) += op(m_rhs1) * op(m_rhs2) do not match");
	if(cols(m_lhs_rows,m_lhs_cols,trans_lhs) != cols(m_rhs2_rows,m_rhs2_cols,trans_rhs2))
		throw std::length_error("Mp_MtM_assert_sizes(...) : The number of columns in"
			" m_lhs and the number of columns in m_rhs1 in the operation"
			" op(m_lhs) += op(m_rhs1) * op(m_rhs2) do not match");
}

#endif
