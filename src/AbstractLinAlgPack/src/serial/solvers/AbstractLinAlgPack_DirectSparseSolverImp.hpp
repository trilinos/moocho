// /////////////////////////////////////////////////////////////
// DirectSparseSolverImp.hpp
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

#ifndef	DIRECT_SPARSE_FORTRAN_COMPATIBLE_SOLVER_IMP_H
#define DIRECT_SPARSE_FORTRAN_COMPATIBLE_SOLVER_IMP_H

#include "DirectSparseSolver.hpp"
#include "SparseLinAlgPack/src/VectorSpaceSerial.hpp"

namespace SparseSolverPack {

class DirectSparseSolverImp;

///
/** Implementation node class for \c DirectSparseSolver that takes
 * care of the memory management details.
 *
 * ToDo: Finish documentation!
 *
 * Subclasses must override the following methods which are pure
 * virtual and are not implemented here: \c basis_matrix_factory(),
 * \c estimated_fillin_ratio(), \c create_fact_struc(),
 * \c create_fact_nonzeros(), \c imp_analyze_and_factor() and
 * \c imp_factor().
 */
class DirectSparseSolverImp : public DirectSparseSolver {
public:

	/** @name Protected types */
	//@{

	///
	/** Abstract class for objects that represent the factorization
	  * nonzeros of a particular matrix.
	  *
	  * This storage for the nonzeros can be reused over and over again
	  * for factored matrices with the same structure.
	  */
	class FactorizationNonzeros {
	public:
		///
		virtual ~FactorizationNonzeros() {}
	};

	///
	/** Implementation node subclass that combines factorization structure and
	 * factorization nonzeros into a single basis matrix object.
	 *
	 * The only methods that subclasses must override for a complete
	 * basis matrix object are:
	 * 
	 * <tt>AbstractLinAlgPack::MatrixNonsingular::V_InvMtV()</tt>,
	 * and  <tt>create_matrix()</tt>.
	 *
	 * Note that is if very important that subclasses do not maintain their
	 * own copies of smart reference counting pointer objects to the
	 * factorization structure or factorization nonzeros.  Instead, the 
	 * methods \c get_fact_struc() and \c get_fact_nonzeros() should
	 * be called inside of the implementation of the the \c V_InvMtV()
	 * method and then \c dynamic_cast<> used to get at the concreate
	 * objects.
	 */
	class BasisMatrixImp : public BasisMatrix {
	public:

		/** @name Public types */
		//@{

		///
		typedef MemMngPack::ref_count_ptr<FactorizationNonzeros> fact_nonzeros_ptr_t;

		//@}

		/** @name Access */
		//@{

		///
		/** Return a reference to a smart pointer to the object that represents
		 * the factorization nonzeros.
		 *
		 * Returning a reference to a \c ref_count_ptr<> object verses returning
		 * a \c ref_count_ptr<> object itself is critical so that we can rely on
		 * \c ref_count_ptr<>::count() to tell us how many clients have a reference
		 * to this object.
		 */
		virtual const fact_nonzeros_ptr_t&  get_fact_nonzeros() const;

		//@}

		/** @name Overridden from MatrixBase */
		//@{

		///
		const VectorSpace& space_cols() const;
		///
		const VectorSpace& space_rows() const;
		///
		size_type rows() const;
		///
		size_type cols() const;

		//@}

		/** @name Overridden from MatrixNonsinguar */
		//@{

		///
		mat_mns_mut_ptr_t clone_mns();

		//@}

		/** @name Overridden from BasisMatrix */
		//@{

		///
		virtual const fact_struc_ptr_t&  get_fact_struc() const;

		//@}

	protected:

		/** @name Constructors/initailizers */
		//@{

		///
		/** Default initializers to uninitialized.
		 */
		BasisMatrixImp();

		///
		/** Calls <tt>this->initialize()</tt>
		 */
		BasisMatrixImp(
			size_type                      dim
			,const fact_struc_ptr_t        &fact_struc
			,const fact_nonzeros_ptr_t     &fact_nonzeros
			);

		///
		/** Initialize given initialized factorization structure and factorization nonzeros objects.
		 */
		virtual void initialize(
			size_type                      dim
			,const fact_struc_ptr_t        &fact_struc
			,const fact_nonzeros_ptr_t     &fact_nonzeros
			);

		///
		/** Make uninitialized.
		 *
		 * Postconditions:<ul>
		 * <li> <tt>this->dim() == 0</tt>
		 * <li> <tt>this->get_fact_struc().get() == NULL</tt>
		 * <li> <tt>this->get_fact_nonzeros().get() == NULL</tt>
		 * </ul>
		 */
		void set_uninitialized();

		//@}

		/** @name Factory methods to be overridden by subclasses */
		//@{

		/// Called by \c this->clone-mns(). 
		virtual MemMngPack::ref_count_ptr<BasisMatrixImp> create_matrix() const = 0;

		//@}

	private:

		///
		/** Allow only DirectSparseSolverImp to initialize objects of this type.
		 * Important !!!!! Even though DirectSparseSolverImp has access to these
		 * private members it is strictly not to access them directly !!!!
		 */
		friend class DirectSparseSolverImp;

#ifdef DOXYGEN_COMPILE
		FactorizationStructure    *fact_struc;
		FactorizationNonzeros     *fact_nonzeros;
#else
		size_type                 dim_;
		fact_struc_ptr_t          fact_struc_;
		fact_nonzeros_ptr_t       fact_nonzeros_;
		VectorSpaceSerial         vec_space_;
#endif

	}; // end class BasisMatrixImp

	//@}

	/** @name Overridden from DirectSparseSolver */
	//@{

	///
	void analyze_and_factor(
		const SparseLinAlgPack::MatrixConvertToSparse   &A
		,DenseLinAlgPack::IVector                            *row_perm
		,DenseLinAlgPack::IVector                            *col_perm
		,size_type                                      *rank
		,BasisMatrix                                    *basis_matrix
		,std::ostream                                   *out
		);
	///
	void factor(
		const SparseLinAlgPack::MatrixConvertToSparse   &A
		,BasisMatrix                                    *basis_matrix
		,const BasisMatrix::fact_struc_ptr_t            &fact_struc
		,std::ostream                                   *out
		);
	///
	const BasisMatrix::fact_struc_ptr_t& get_fact_struc() const;
	///
	void set_uninitialized();

	//@}

protected:

	/** @name Protected pure virtual methods to be overridden by concrete direct solver subclasses */
	//@{

	///
	/** Create a new, uninitialized \c FactorizationStructure object.
	 */
	virtual const MemMngPack::ref_count_ptr<FactorizationStructure> create_fact_struc() const = 0;

	///
	/** Create a new, uninitialized \c FactorizationNonzeros object.
	 */
	virtual const MemMngPack::ref_count_ptr<FactorizationNonzeros> create_fact_nonzeros() const = 0;

	///
	/** Called to implement the \c analyze_and_factor() without having to worry about
	 * memory mangagment details.
	 *
	 * ToDo: Finish documentation!
	 */	
	virtual void imp_analyze_and_factor(
		const SparseLinAlgPack::MatrixConvertToSparse   &A
		,FactorizationStructure                         *fact_struc
		,FactorizationNonzeros                          *fact_nonzeros
		,DenseLinAlgPack::IVector                            *row_perm
		,DenseLinAlgPack::IVector                            *col_perm
		,size_type                                      *rank
		,std::ostream                                   *out            = NULL
		) = 0;

	///
	/** Called to implement the \c analyze_and_factor() without having to worry about
	 * memory mangagment details.
	 *
	 * ToDo: Finish documentation!
	 */	
	virtual void imp_factor(
		const SparseLinAlgPack::MatrixConvertToSparse   &A
		,const FactorizationStructure                   &fact_struc
		,FactorizationNonzeros                          *fact_nonzeros
		,std::ostream                                   *out            = NULL
		) = 0;

	//@}

private:

#ifdef DOXYGEN_COMPILE
		FactorizationStructure           *fact_struc;
#else
		BasisMatrix::fact_struc_ptr_t    fact_struc_;
		size_type                        rank_;
#endif

};	// end class DirectSparseSolverImp 

}	// end namespace SparseSolverPack 

#endif	// DIRECT_SPARSE_FORTRAN_COMPATIBLE_SOLVER_IMP_H
