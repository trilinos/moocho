// ///////////////////////////////////////////////////////////////////////
// NLPSerialPreprocessExplJac.h
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

#ifndef NLP_SERIAL_PREPROCESS_EXPL_JAC_H
#define NLP_SERIAL_PREPROCESS_EXPL_JAC_H

#include <valarray>

#include "NLPSerialPreprocess.h"
#include "NLPFirstOrderInfo.h"
#include "LinAlgPack/include/VectorClass.h"
#include "AbstractFactory.h"

namespace NLPInterfacePack {

///
/** NLP node subclass complementing \c NLPSerialPreprocess for explicit Jacobians.
 *
 * This subclass does a lot of work.  It has to consider several different
 * types of variability.  The matrices \c Gc and \c Gh that are computed must
 * take into consideration whether or not inequalities are converted
 * to equalities (<tt>convert_inequ_to_equ</tt>) and the permutation
 * of the entries according to the current basis selection.
 *
 *  When <tt>convert_inequ_to_equ == false</tt> then:
 \verbatim
    
    Gc = P_var * Gc_orig * P_equ'   (only if m_orig  > 0)

    Gh = P_var * Gh_orig            (only if mI_orig > 0)
 \endverbatim
 * However, then <tt>convert_inequ_to_equ == true</tt> then:
 \verbatim

    Gc = P_var * [  Gc_orig    Gh_orig   ] * P_equ'
                 [     0          -I     ]

    Gf_full = Empty
 \endverbatim
 *
 * ToDo: Finish documentation!
 *
 * <b>Subclass developers</b>
 *
 * Subclass developer's don't have to worry about slack variables or basis
 * permutations.  A concreate subclass just has to override the functions
 * that defined the original %NLP (see the tutorial example %NLP ???).
 *
 * In addition to the methods that must be overridden in \c NLPSerialPreprocess
 * (<A HREF="classNLPInterfacePack_1_1NLPSerialPreprocess.html#must_override">see</A>)
 * the following methods must be overridden as well: \c imp_Gc_nz_orig(), \c imp_Gh_nz_orig(),
 * \c imp_calc_Gc_orig(), \c imp_calc_Gh_orig().
 */
class NLPSerialPreprocessExplJac
	: virtual public NLPSerialPreprocess
	, virtual public NLPFirstOrderInfo
{
public:

	/** @name Public types */
	//@{
	
	///
	typedef MemMngPack::ref_count_ptr<
		const MemMngPack::AbstractFactory<MatrixWithOp> >    factory_mat_ptr_t;

	//@}

	/** @name Constructors / initializers */
	//@{

	///
	/** Calls <tt>this->set_mat_factories()</tt>.
	 */
	NLPSerialPreprocessExplJac(
		const factory_mat_ptr_t     &factory_Gc_full = MemMngPack::null
		,const factory_mat_ptr_t    &factory_Gh_full = MemMngPack::null
		);

	///
	/** Initialize with matrix factories for original matrices \c Gc and \c Gh.
	 *
	 * These matrix types will be used for \c AbstractLinAlgPack::MatrixPermAggr::mat_orig()
	 * returned by the initialized \c Gc and \c Gh objects respectively.
	 *
	 * @param  factory_Gc_full
	 *                [in] Smart pointer to matrix factory for \c Gc_full.  If
	 *                <tt>factory_Gc_full.get() == NULL</tt> then the concrete matrix
	 *                type ??? will be used as the default.
	 * @param  factory_Gh_full
	 *                [in] Smart pointer to matrix factory for \c Gh_full.  If
	 *                <tt>factory_Gh_full.get() == NULL</tt> then the concrete matrix
	 *                type ??? will be used as the default. */
	void set_mat_factories(
		const factory_mat_ptr_t     &factory_Gc_full
		,const factory_mat_ptr_t    &factory_Gh_full
		);

	//@}

	/** @name Overridden public members from NLP */
	//@{

	///
	void initialize(bool test_setup);	
	///
	bool is_initialized() const;

	//@}

	/** @name Overridden public members from NLPFirstOrderInfo */
	//@{
	
	///
	const mat_fcty_ptr_t factory_Gc() const;
	///
	const mat_fcty_ptr_t factory_Gh() const;
	/// Validates the type of Gc is correct
	void set_Gc(MatrixWithOp* Gc);
	/// Validates the type of Gh is correct
	void set_Gh(MatrixWithOp* Gh);

	//@}

	/** @name Overridden public members from NLPVarReductPerm */
	//@{

	///
	bool get_next_basis(
		Permutation*  P_var,   Range1D* var_dep
		,Permutation* P_equ,   Range1D* equ_decomp
		,Permutation* P_inequ, Range1D* inequ_decomp
		);
	///
	void set_basis(
		const Permutation   &P_var,   const Range1D  &var_dep
		,const Permutation  *P_equ,   const Range1D  *equ_decomp
		,const Permutation  *P_inequ, const Range1D  *inequ_decomp
		);

	//@}


protected:

	/** @name Overridden protected members from NLPFirstOrderInfo */
	//@{

	///
	void imp_calc_Gc(
		const VectorWithOp& x, bool newx
		,const FirstOrderInfo& first_order_info
		) const;
	///
	void imp_calc_Gh(
		const VectorWithOp& x, bool newx
		,const FirstOrderInfo& first_order_info
		) const;
	
	//@}

	/** @name Protected types */
	//@{

	///
	/** Struct for zero and explicit first order quantities that subclass must fill in.
	 *
	 * When computing Gc and/or Gh, the subclass can be instructed to set the row and columns
	 * index arrays by setting \c Gc_ivect==NULL and \c Gh_ivect==NULL or not respecitively.
	 */
	struct FirstOrderExplInfo {
		///
		typedef std::valarray<value_type>    val_t;
		///
		typedef std::valarray<index_type>    ivect_t;
		//
		typedef std::valarray<index_type>    jvect_t;
		///
		FirstOrderExplInfo()
			:Gc_val(NULL), Gc_ivect(NULL), Gc_jvect(NULL)
			,Gh_val(NULL), Gh_ivect(NULL), Gh_jvect(NULL)
			,f(NULL)
		{}
		///
		FirstOrderExplInfo(
			index_type* Gc_nz_in, val_t* Gc_val_in, ivect_t* Gc_ivect_in, jvect_t* Gc_jvect_in
			,index_type* Gh_nz_in, val_t* Gh_val_in, ivect_t* Gh_ivect_in, jvect_t* Gh_jvect_in
			,const ObjGradInfoSerial& obj_grad
			)
			:Gc_nz(Gc_nz_in), Gc_val(Gc_val_in), Gc_ivect(Gc_ivect_in), Gc_jvect(Gc_jvect_in)
			,Gh_nz(Gh_nz_in), Gh_val(Gh_val_in), Gh_ivect(Gh_ivect_in), Gh_jvect(Gh_jvect_in)
			,Gf(obj_grad.Gf), f(obj_grad.f), c(obj_grad.c), h(obj_grad.h)
		{}
		///
		size_type*    Gc_nz;
		///
		val_t*        Gc_val;
		///
		ivect_t*      Gc_ivect;
		///
		jvect_t*      Gc_jvect;
		///
		size_type*    Gh_nz;
		///
		val_t*        Gh_val;
		///
		ivect_t*      Gh_ivect;
		///
		jvect_t*      Gh_jvect;
		///
		Vector*       Gf;
		///
		value_type*   f;
		///
		Vector*       c;
		///
		Vector*       h;
	}; // end struct FirstOrderExplInfo

	//@}

	/** @name Pure virtual template methods to be defined by subclasses */
	//@{

	///
	/** Return the number of nonzero elements in \c Gc before elements are removed for fixed variables.
	  *
	  * The value returned from this method before the first time \c imp_calc_Gc() is called
	  * is an upper estimate of the number of nonzeros.  To get the actual number
	  * of nonzeros, call this function again after \c imp_calc_Gc() has been called.
	  */
	virtual size_type imp_Gc_nz_orig() const = 0;

	///
	/** Return the number of nonzero elements in \c Gh before elements are removed for fixed variables.
	  *
	  * The value returned from this method before the first time \c imp_calc_Gh() is called
	  * is an upper estimate of the number of nonzeros.  To get the actual number
	  * of nonzeros, call this function again after \c imp_calc_Gh() has been called.
	  */
	virtual size_type imp_Gh_nz_orig() const = 0;

	///
	/** Calculate the COOR matrix for the gradient for all of the \a c(x) constaints in the orig %NLP.
	 *
	 * @param x       [in]  Unknown vector (size n_full).
	 * @param newx    [in]  True if is a new point.
	 * @param first_order_expl_info
	 *                [out] Pointers to zero and first order quantities .
	 *                On output, <tt>*first_order_expl_info.Gc_nz</tt> must be set to the actual
	 *                number of nonzero elements in \c Gc and the array of nonzero entry
	 *                values <tt>*first_order_expl_info.Gc_val</tt> must also be set.  In addition,
	 *                if <tt>this->multi_calc() == true</tt> then
	 *                any of the other quantities pointed to in \c first_order_expl_info may be set on
	 *                output, but are not guaranteed to be.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>first_order_expl_info.Gc_nz != NULL</tt>
	 * <li> <tt>first_order_expl_info.Gc_val != NULL</tt>
	 * <li> <tt>(first_order_expl_info.Gc_ivect != NULL) == (first_order_expl_info.Gc_jvect != NULL)</tt> 
	 * </ul>
	 *
	 * Preconditions:<ul>
	 * <li> <tt>*first_order_expl_info.Gc_nz</tt> is updated to number of nonzero elements set in
	 *      <tt>*first_order_expl_info.Gc_val</tt>.
	 * <li> <tt>(*first_order_expl_info.Gc_val)[k]</tt>, for <tt>k = 1...*first_order_expl_info.Gc_nz</tt>
	 *      is set to the nonzero entry values in \c Gc.
	 * <li> [<tt>first_order_expl_info.Gc_ivect != NULL</tt>]
	 *      <tt>(*first_order_expl_info.Gc_ivect)[k]</tt>, for <tt>k = 1...*first_order_expl_info.Gc_nz</tt>
	 *      is set to the row indexes for the nonzero entires in \c Gc.
	 * <li> [<tt>first_order_expl_info.Gc_jvect != NULL</tt>]
	 *      <tt>(*first_order_expl_info.Gc_jvect)[k]</tt>, for <tt>k = 1...*first_order_expl_info.Gc_nz</tt>
	 *      is set to the column indexes for the nonzero entires in \c Gc.
	 * </ul>
	 *
	 * Note that duplicate entires with the same row and column indexes are allowed.  In this case, the
	 * matrix entries are considered to be summed.
	 */
	virtual void imp_calc_Gc_orig(
		const VectorSlice& x_full, bool newx
		, const FirstOrderExplInfo& first_order_expl_info
		) const = 0;

	///
	/** Calculate the COOR matrix for the gradient for all of the \c h(x) constaints in the orig %NLP.
	 *
	 * @param x       [in]  Unknown vector (size n_full).
	 * @param newx    [in]  True if is a new point.
	 * @param first_order_expl_info
	 *                [out] Pointers to zero and first order quantities .
	 *                On output, <tt>*first_order_expl_info.Gh_nz</tt> must be set to the actual
	 *                number of nonzero elements in \c Gh and the array of nonzero entry
	 *                values <tt>*first_order_expl_info.Gh_val</tt> must also be set.  In addition,
	 *                if <tt>this->multi_calc() == true</tt> then
	 *                any of the other quantities pointed to in \c first_order_expl_info may be set on
	 *                output, but are not guaranteed to be.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>first_order_expl_info.Gh_nz != NULL</tt>
	 * <li> <tt>first_order_expl_info.Gh_val != NULL</tt>
	 * <li> <tt>(first_order_expl_info.Gh_ivect != NULL) == (first_order_expl_info.Gh_jvect != NULL)</tt> 
	 * </ul>
	 *
	 * Preconditions:<ul>
	 * <li> <tt>*first_order_expl_info.Gh_nz</tt> is updated to number of nonzero elements set in
	 *      <tt>*first_order_expl_info.Gh_val</tt>.
	 * <li> <tt>(*first_order_expl_info.Gh_val)[k]</tt>, for <tt>k = 1...*first_order_expl_info.Gh_nz</tt>
	 *      is set to the nonzero entry values in \c Gh.
	 * <li> [<tt>first_order_expl_info.Gh_ivect != NULL</tt>]
	 *      <tt>(*first_order_expl_info.Gh_ivect)[k]</tt>, for <tt>k = 1...*first_order_expl_info.Gh_nz</tt>
	 *      is set to the row indexes for the nonzero entires in \c Gh.
	 * <li> [<tt>first_order_expl_info.Gh_jvect != NULL</tt>]
	 *      <tt>(*first_order_expl_info.Gh_jvect)[k]</tt>, for <tt>k = 1...*first_order_expl_info.Gh_nz</tt>
	 *      is set to the column indexes for the nonzero entires in \c Gh.
	 * </ul>
	 *
	 * Note that duplicate entires with the same row and column indexes are allowed.  In this case, the
	 * matrix entries are considered to be summed.
	 */
	virtual void imp_calc_Gh_orig(
		const VectorSlice& x_full, bool newx
		, const FirstOrderExplInfo& first_order_expl_info
		) const = 0;

	//@}

	/** @name Protected member functions for subclasses to use */
	//@{

	/// Assert if we have been initizlized (throws UnInitialized)
	void assert_initialized() const;

	///
	const FirstOrderExplInfo first_order_expl_info() const;

	//@}

private:

	// ////////////////////////////////////////
	// Private data members
	
	bool initialized_;              // Flag for if the NLP has has been properly initialized
	bool test_setup_;               // Flag for if to test the setup of things or not

	factory_mat_ptr_t   factory_Gc_full_;
	factory_mat_ptr_t   factory_Gh_full_;
	mat_fcty_ptr_t      factory_Gc_;
	mat_fcty_ptr_t      factory_Gh_;

	mutable size_type   Gc_nz_orig_;    // Number of nonzeros in the original NLP Gc
	mutable size_type   Gh_nz_orig_;    // Number of nonzeros in the original NLP Gh
	mutable size_type   Gc_nz_full_;    // Number of nonzeros in the full NLP Gc
	mutable size_type   Gh_nz_full_;    // Number of nonzeros in the full NLP Gh
	mutable FirstOrderExplInfo::val_t    Gc_val_orig_;   // Storage for explicit nonzeros of full Gc
	mutable FirstOrderExplInfo::ivect_t  Gc_ivect_orig_;
	mutable FirstOrderExplInfo::jvect_t  Gc_jvect_orig_;
	mutable FirstOrderExplInfo::val_t    Gh_val_orig_;   // Storage for explicit nonzeros of orig Gh
	mutable FirstOrderExplInfo::ivect_t  Gh_ivect_orig_;
	mutable FirstOrderExplInfo::jvect_t  Gh_jvect_orig_;

	mutable bool                         Gc_perm_new_basis_updated_;  // Flag for if a new basis was set!
	mutable bool                         Gh_perm_new_basis_updated_;  // Flag for if a new basis was set!

	// ////////////////////////////
	// Private member functions

	//
	void imp_calc_Gc_or_Gh(
		bool calc_Gc
		,const VectorWithOp& x, bool newx
		,const FirstOrderInfo& first_order_info
		) const;

	//
	void imp_fill_jacobian_entries(
		size_type           n             // [in]
		,size_type          n_full        // [in]
		,bool               load_struct   // [in] If true, then the structure is loaded also
		,const index_type   col_offset    // [in] Offset for filled column indexes
		,const value_type   *val_full     // [in] Values (!=NULL)
		,const value_type   *val_full_end // [in] Values end (!=NULL)
		,const index_type   *ivect_full   // [in] Row indexes (!=NULL)
		,const index_type   *jvect_full   // [in] Column indexes (!=NULL)
		,index_type         *nz           // [in/out] Number of nonzeros added (!=NULL)            
		,value_type         *val_itr      // [out] Values to fill (!=NULL)
		,index_type         *ivect_itr    // [out] Row indexes (can be NULL if load_struct == false)
		,index_type         *jvect_itr    // [out] Column indexes  (can be NULL if load_struct == false)
		) const;

};	// end class NLPSerialPreprocessExplJac

// ///////////////////////////
// inline members

inline
const NLPSerialPreprocessExplJac::FirstOrderExplInfo
NLPSerialPreprocessExplJac::first_order_expl_info() const
{
	return FirstOrderExplInfo(
		&Gc_nz_orig_,&Gc_val_orig_,&Gc_ivect_orig_,&Gc_jvect_orig_
		,&Gh_nz_orig_,&Gh_val_orig_,&Gh_ivect_orig_,&Gh_jvect_orig_
		,obj_grad_orig_info()
		);
}

}	// end namespace NLPInterfacePack 

#endif // NLP_SERIAL_PREPROCESS_EXPL_JAC_H
