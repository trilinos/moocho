// //////////////////////////////////////////////////////////////////////////////////
// MatrixSymPosDefLBFGS.h
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

#ifndef MATRIX_SYM_POS_DEF_LBFGS_H
#define MATRIX_SYM_POS_DEF_LBFGS_H

#include <vector>

#include "MatrixSymSecantUpdateable.h"
#include "MatrixSymAddDelUpdateable.h"
#include "SparseLinAlgPack/include/MatrixSymWithOpFactorized.h"
#include "LinAlgPack/include/GenMatrixAsTriSym.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ConstrainedOptimizationPack {

///
/** Implementation of limited Memory BFGS matrix.
 *
 * The function set_num_updates_stored(l) must be called first to set the maximum number of
 * the most recent updates that can  be stored.  The storage requirements for this class are
 * O( n_max*l + l*l ) which is O(n_max*l) when n_max >> l which is expected.
 *
 * This implementation is based on:
 *
 * Byrd, Nocedal, and Schnabel, "Representations of quasi-Newton matrices
 * and their use in limited memory methods", Mathematical Programming, 63 (1994)
 * 
 * Consider BFGS updates of the form:
 \begin{verbatim}

 ( B^{k-1}, s^{k-1}, y^{k-1} ) -> B^{k}
 
 where:
 
 B^{k} = B^{k-1} - ( (B*s)*(B*s)' / (s'*B*s) )^{k-1} + ( (y*y') / (s'*y) )^{k-1}
 
 B <: R^(n x n)
 s <: R^(n)
 y <: R^(n)
 
 \end{verbatim}
 * Given that we start from the same initial matrix #Bo#, the updated matrix #B^{k}#
 * will be the same independent of the order the #(s^{i},y^{i})# updates are added.
 *
 * Now let us consider limited memory updating.  For this implementation we set:
 \begin{verbatim}

 Bo = ( 1 / gamma_k ) * I
 
 where:
               / (s^{k-1}'*y^{k-1})/(y^{k-1}'*y^{k-1})           :  if auto_rescaling() == true
     gamma_k = |
	           \  alpha from last call to init_identity(n,alpha) :  otherwise
 
 \end{verbatim}
 * Now let us define the matrices #S# and #Y# that store the update vectors
 * #s^{i}# and #y^{i}# for #i = 1 ... m_bar#:
 \begin{verbatim}

 S = [ s^{1}, s^{2},...,s^{m_bar} ] <: R^(n x m)
 Y = [ y^{1}, y^{2},...,y^{m_bar} ] <: R^(n x m)
 
 \end{verbatim}
 * Here we are only storing the #m_bar <= m# most recent update vectors
 * and their ordering is arbitrary.  The columns #S(:,k_bar)# and #Y(:,k_bar)#
 * contain the most recent update vectors.  The next most recent vectors
 * are to the left (i.e. #p = k_bar-1#) and so forth until #p = 1#.  Then
 * the next most recent update vectors start at #m_bar# and move to the
 * left until you reach the oldest update vector stored at column #k_bar+1#.
 * This is all client need to know in order to reconstruct the updates
 * themselves.
 *
 * This class allows matrix vector products x = B*y and the inverse
 * matrix vector products x = inv(B)*y to be performed at a cost of
 * about O(n*m^2).
 *
 * In addition, the class supports the MatixSymAddDelUpdateable interface
 * with a few major restrictions.  This allows the client to add and
 * remove rows and columns from the matrix.
 */
class MatrixSymPosDefLBFGS
	: public MatrixSymWithOpFactorized
	, public MatrixSymSecantUpdateable
	, public MatrixSymAddDelUpdateable
{
public:

	// //////////////////////////////////////////////
	// Constructors and initializers

	/// Calls initial_setup(,,,)
	MatrixSymPosDefLBFGS(
		size_type   max_size           = 0
	    ,size_type  m                  = 10
		,bool       maintain_original  = true
		,bool       maintain_inverse   = true
		,bool       auto_rescaling     = false
		);

	///
	/** Set whether automatic rescaling is used or not.
	  *
	  * This function must be called before a BFGS update is performed
	  * in order for it to take effect for that update.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, auto_rescaling )

	///
	/** Initial setup for the matrix.
	  *
	  * This function must be called before init_identity(n)
	  * is called.  When this function is called all current
	  * updates are lost and the matrix becomes uninitialized.
	  *
	  * @param  max_size
	  *            [in] If max_size > 0 then this is the max size
	  *            the matrix is allowed to become.  If max_size == 0
	  *            then this size will be determined by one of the
	  *            initialization methods.
	  * @param  m  [in] Max number of recent update vectors stored
	  * @param  maintain_original
	  *            [in] If true then quantities needed to compute
	  *            x = Bk*y will be maintained, otherwise they
	  *            will not be unless needed.  This is to save
	  *            computational costs in case matrix vector
	  *            products will never be needed.  However,
	  *            if a matrix vector product is needed then
	  *            these quantities will be computed on the fly
	  *            in order to satisfy the request.
	  * @param  maintain_inverse
	  *            [in] If true then quantities needed to compute
	  *            x = inv(Bk)*y = x = Hk*y will be maintained
	  *            , otherwise they will not be unless needed.
	  *            This is to save computational costs in case
	  *            inverse matrix vector products will never be needed.
	  *            However, if the inverse product is ever needed
	  *            then the needed quantities will be computed
	  *            on the fly in order to satisfiy the request.
	  *            Because it takes so little extra work to maintain
	  *            the quantities needed for Hk it is recommended
	  *            to always set this to true.
	  * @param  auto_rescaling
	  *            [in] See intro.
	  */
	 void initial_setup(
		 size_type   max_size           = 0
		 ,size_type  m                  = 10
		 ,bool       maintain_original  = true
		 ,bool       maintain_inverse   = true
		 ,bool       auto_rescaling     = false
		 );

	// //////////////////////////////////
	// Representation access

	///
	size_type m() const;
	///
	size_type m_bar() const;
	///
	size_type k_bar() const;
	///
	value_type gamma_k() const;
	///
	const GenMatrixSlice S() const;
	///
	const GenMatrixSlice Y() const;
	///
	bool maintain_original() const;
	///
	bool maintain_inverse() const;
	/// Returns the total number of successful secant updates performed.
	size_type num_secant_updates() const;

	// /////////////////////////////////////////////////////
	// Overridden from Matrix

	///
	size_type rows() const;

	// /////////////////////////////////////////////////////////
	/** @name Overridden from MatrixWithOp */
	//@{

	///
	std::ostream& output(std::ostream& out) const;
	///
	MatrixWithOp& operator=(const MatrixWithOp& m);
	///
	void Vp_StMtV(VectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2, value_type beta) const;

	//@}

	// ////////////////////////////////////////////////////////////
	/** @name Overridden from MatrixWithOpFactorized */
	//@{

	///
	void V_InvMtV( VectorSlice* v_lhs, BLAS_Cpp::Transp trans_rhs1
		, const VectorSlice& vs_rhs2) const;

	//@}

	// ///////////////////////////////////////////////////////////
	/** @name Overridden from MatrixSymSecantUpdateable */
	//@{

	///
	void init_identity( size_type n, value_type alpha );
	///
	/** Actually this calls init_identity( (&diag)->size(), norm_inf(diag) ).
	  *
	  * This initialization is not convienent for this implementation.
	  * Besides, when we are using automatric rescaling (auto_rescaling == true)
	  * then this will really not matter much anyway.
	  */
	void init_diagonal( const VectorSlice& diag );
	///
	void secant_update(VectorSlice* s, VectorSlice* y, VectorSlice* Bs);

	//		end Overridden from MatrixSymSecantUpdateable
	//@}

	// ////////////////////////////////////////////////////////
	/** @name Overridden from MatrixSymAddDelUpdateble */
	//@{

	/// This is fine as long as alpha > 0.0.
	void initialize(
		value_type         alpha
		,size_type         max_size
		);
	/// Sorry, this will throw an exception!
	void initialize(
		const sym_gms      &A
		,size_type         max_size
		,bool              force_factorization
		,Inertia           inertia
		,PivotTolerances   pivot_tols
		);
	///
	size_type max_size() const;
	/// Returns (0,0,rows())
	Inertia inertia() const;
	/// Will set rows() == 0
	void set_uninitialized();
	///
	/** Augment the matrix to add a row and column.
	 *
	 * This function is very limited in what it will do.
	 * It will throw exceptions if alpha <= 0.0 or t != NULL
	 * or add_eigen_val == EIGEN_VAL_NEG or this->rows() == this->max_size().
	 * The obvious postconditions for this function will only technically
	 * be satisfied if alpha == this->gamma_k().
	 */
	void augment_update(
		const VectorSlice  *t
		,value_type        alpha
		,bool              force_refactorization
		,EEigenValType     add_eigen_val
		,PivotTolerances   pivot_tols
		);
	/// Should always succeed unless user gives a wrong value for drop_eigen_val.
	void delete_update(
		size_type          jd
		,bool              force_refactorization
		,EEigenValType     drop_eigen_val
		,PivotTolerances   pivot_tols
		);
	
	//@}
	
private:

	// //////////////////////////////////
	// Private types

	// //////////////////////////////////
	// Private data members

	bool        maintain_original_;  // If true, qualities needed for Bk will be maintained
	bool        original_is_updated_;// If true, qualities needed for Bk are already updated
	bool        maintain_inverse_;   // If true, quantities needed for Hk will be maintained
	bool        inverse_is_updated_; // If true, quantities needed for Hk are already updated

	size_type	n_max_,	// The maximum size the matrix is allowed to become.
				n_,		// Size of the matrix.  If 0 then is uninitialized
				m_,		// Maximum number of update vectors that can be stored.
				m_bar_,	// Current number of update vectors being stored.
						// 0 <= m_bar <= m
				k_bar_,	// Position of the most recently stored update vector in S & Y
						// 1 <= k_bar <= m_bar
		        num_secant_updates_; // Records the number of secant updates performed
	value_type	gamma_k_;// Scaling factor for Bo = (1/gamma_k) * I.

	GenMatrix	S_,		// (n_max x m) Matrix of stored update vectors = [ s1, ..., sm ]
						// S(:,k_bar) is the most recently stored s update vector
				Y_,		// (n_max x m) Matrix of stored update vectors = [ y1, ..., ym ]
						// Y(:,k_bar) is the most recently stored y update vector
				STY_,	// (m x m) The matrix S'Y
				STSYTY_;// ((m+1) x (m+1)) The strictly upper triangular part stores the
						// upper triangular part Y'Y and the strictly lower triangular
						// part stores the lower triangular part of S'S.  The diagonal
						// can be used for workspace.

	mutable bool		Q_updated_;	// True if Q has been updated for the most current update.
	mutable GenMatrix	QJ_;		// Used to store factorization of the schur complement of Q.

	mutable Vector		work_;	// workspace for performing operations.

	// //////////////////////////////////
	// Private member functions

	// Access to important matrices.

	///
	const tri_gms R() const;
	/// Strictly lower triangular part of L
	const tri_gms Lb() const;
	///
	GenMatrixSlice STY();
	///
	const GenMatrixSlice STY() const;
	///
	sym_gms STS();
	///
	const sym_gms STS() const;
	///
	sym_gms YTY();
	///
	const sym_gms YTY() const;
	/// y = inv(Q) * x
	void V_invQtV( VectorSlice* y, const VectorSlice& x ) const;
	/// y += D * x
	void Vp_DtV( VectorSlice* y, const VectorSlice& x ) const;

	// Updates

	/// Update Q
	void update_Q() const;

	///
	void assert_initialized() const;

};	// end class MatrixSymPosDefLBFGS

// //////////////////////////////////////////////
// Inline member functions

inline
size_type MatrixSymPosDefLBFGS::m() const
{
	return m_;
}

inline
size_type MatrixSymPosDefLBFGS::m_bar() const
{
	return m_bar_;
}

inline
size_type MatrixSymPosDefLBFGS::k_bar() const
{
	return k_bar_;
}

inline
value_type MatrixSymPosDefLBFGS::gamma_k() const
{
	return gamma_k_;
}

inline
const GenMatrixSlice MatrixSymPosDefLBFGS::S() const
{
	return S_(1,n_,1,m_bar_);
}

inline
const GenMatrixSlice MatrixSymPosDefLBFGS::Y() const
{
	return Y_(1,n_,1,m_bar_);
}

inline
bool MatrixSymPosDefLBFGS::maintain_original() const
{
	return maintain_original_;
}

inline
bool MatrixSymPosDefLBFGS::maintain_inverse() const
{
	return maintain_inverse_;
}

inline
size_type MatrixSymPosDefLBFGS::num_secant_updates() const
{
	return num_secant_updates_;
}


}	// end namespace ConstrainedOptimizationPack 

#endif	// MATRIX_SYM_POS_DEF_LBFGS_H
