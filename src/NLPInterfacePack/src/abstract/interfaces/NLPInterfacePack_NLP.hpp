// ///////////////////////////////////////////////////////////////////////
// NLP.h
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

#ifndef NLP_H
#define NLP_H

#include <stdexcept>
#include <string>

#include "NLPInterfacePackTypes.h"
#include "StandardCompositionRelationshipsPack.h"
#include "ref_count_ptr.h"

namespace NLPInterfacePack {

///
/** NLP interface class {abstract}.
 *
 * This class represents an abstract interface to a general nonlinear
 * programming problem of the form: \\
 * \\
 \begin{verbatim}
     min     f(x)
     s.t.    c(x) = 0
            hl <= h(x) <= hu
	        xl <= x <= xu
	where:
	        x    <: R^n
	        c(x) <: R^n -> R^m 
	        h(x) <: R^n -> R^mI 
 \end{verbatim}
 * In the above form, none of the variables are fixed between bounds (strictly
 * xl < xu).
 *
 * The Lagrangian for this problem is defined by:\\
 \begin{verbatim}
	L = f(x) + lambda' * c(x)
        + lambdaI_l' * ( hl - h(x) ) + lambdaI_u' * ( h(x) - hu )
        + nul' * ( xl - x ) + nuu' * ( x - xu )
 \end{verbatim}
 *
 * The optimality conditions are given by:
 \begin{verbatim}
	del(L,x)      = del(f,x) + del(c,x) * lambda + del(h,x) * lambdaI + nu = 0
	del(L,lambda) = c(x) = 0
	  where:
        lambdaI = lambdaI_u - lambdaI_l
		nu = nuu - nul
		lambdaI_l(j) * ( h(x)(j) - hl(j) ) = 0,    for j = 1...mI
		lambdaI_u(j) * ( h(x)(j) - hu(j) ) = 0,    for j = 1...mI
		nuu(i) * ( x(i) - xu(i) ) = 0,      for i = 1...n
		nuu(i) * ( x(i) - xu(i) ) = 0,      for i = 1...n
  \end{verbatim}
  *
  * What is unique about this interface is that the vector objects are hidden behind
  * abstact interfaces.  Clients can create vectors from the various vector spaces
  * using the \Ref{VectorSpace} objects returned from \ref{space_x}#()# (dim #n#),
  * \ref{space_c}#()# (dim #m#) and \ref{space_h}#()# (dim #mI#).
  */
class NLP {
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<const VectorSpace>  vec_space_ptr_t;

	/** @name exceptions */
	//@{

	/// Thrown if any member functions are called before initialize() has been called.
	class UnInitialized : public std::logic_error
	{public: UnInitialized(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown from #initialize()# if some logical error occured
	class InvalidInitialization : public std::logic_error
	{public: InvalidInitialization(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if an incompatible object is used
	class IncompatibleType : public std::logic_error
	{public: IncompatibleType(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if #xl()# or #xu()# are called and #num_bounds_x() == 0# 
	class NoBoundsOnVariables : public std::logic_error
	{public: NoBoundsOnVariables(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}

	/** @name Constructors */
	//@{

	/// Initialize to no reference set to calculation quanities
	NLP();

	/// Destructor that cleans all the memory it owns
	virtual ~NLP();

	//@}
	
	///
	/** Initialize the NLP for its first use.
	  *
	  * Postconditions:  If this operation returns without throwing an exception then all of
	  * the member functions can be accessed.
	  *
	  * Note that subclasses must call this function to reset what needs to be
	  * reset here.
	  */
	virtual void initialize();

	///
	/** Return if this is initialized.
	  */
	virtual bool is_initialized() const = 0;
	
	/** @name Dimensionality. */
	//@{

	/// Return the number of variables
	virtual size_type n() const = 0;
	/// Return the number of general equality constraints
	virtual size_type m() const = 0;
	/// Return the number of general inequality constraints
	virtual size_type mI() const = 0;

	//@}

	/** @name Vector Space objects. */
	//@{

	///
	/** Vector space object for unknown variables x (dimension n).
	 */
	virtual vec_space_ptr_t space_x() const = 0;

	///
	/** Vector space object for general equality constraints c(x) (dimension m).
	 *
	 * If #m() ==0# then this method should return #return.get() == NULL#.
	 */
	virtual vec_space_ptr_t space_c() const = 0;

	///
	/** Vector space object for general inequality constraints c(x) (dimension mI).
	 *
	 * If #mI() == 0# then this method should return #return.get() == NULL#.
	 */
	virtual vec_space_ptr_t space_h() const = 0;

	//@}

	/** @name Access Variable Bounds and the Initial Point */
	//@{

	///
	/** Returns the number of variables in #x(i)# for which #xl(i)> -infinite_bound()#
	 * or #xu(i) < +infinite_bound()#.
	 */
	virtual size_type num_bounded_x() const = 0;

	///
	/** Set if the initial point must be within the bounds.
	  */
	virtual void force_xinit_in_bounds(bool force_xinit_in_bounds) = 0;

	///
	/** Returns if the initial point must be within the bounds.
	  */
	virtual bool force_xinit_in_bounds() const = 0;

	///
	/** Returns a reference to the vector of the initial guess for the solution #x#.
	 *
	 * @return  #return.dim() == this->n()# and will be compatible with vectors from
	 * #this->space_x()->create_member()#.
	 */
	virtual const VectorWithOp& xinit() const = 0;

	///
	static value_type infinite_bound();

	/** @name Bounds on the variables #x#.
	 *
	 * If #this->num_bounded_x() == 0# then these methods should
	 * not be called or the excetpion #NoBoundsOnVariables# will be
	 * thrown.  These vectors will be compatible with the vector space
	 * #this->space_x()#.
	 */
	//@{

	///
	/** Returns the lower bounds on the variables #x#.
	 *
	 * Any bounds that are non-existant will return #xl().get_ele(i) == -infinite_bound()#.
	 */
	virtual const VectorWithOp& xl() const = 0;
	///
	/** Returns a reference to the vector of upper bounds on the variables #x#.
	 *
	 * Any bounds that are non-existant will return #xu().get_ele(i) == +infinite_bound()#.
	 */
	virtual const VectorWithOp& xu() const = 0;

	//@}

	/** @name Bounds on the general inequality constraints #h(x)#.
	 *
	 * If #this->mI() == 0# then these methods should
	 * not be called or the excetpion #NoBoundsOnVariables# will be
	 * thrown.  These vectors will be compatible with the vector space
	 * #this->space_h()#.
	 */
	//@{

	///
	/** Returns the lower bounds on the inequality constraints #h(x)#.
	 *
	 * Any bounds that are non-existant will return #hl().get_ele(j) == -infinite_bound()#.
	 */
	virtual const VectorWithOp& hl() const = 0;
	///
	/** Returns upper bounds on the inequality constraints #h(x)#.
	 *
	 * Any bounds that are non-existant will return #hu().get_ele(j) == +infinite_bound()#.
	 */
	virtual const VectorWithOp& hu() const = 0;

	//@}

	///
	/** Get the initial value of the Lagrange multipliers lambda.
	  *
	  * By default this function just sets them to zero.
	  *
	  * Postconditions: \begin{itemize}
	  * \item [lambda != NULL] lambda->size() == this->n()
	  * \item [nu != NULL] nu->size() == this->n()
	  * \end{itemize}
	  *
	  * @param lambda  [out] Pointer to lagrange multipliers for equalities.
	  *                lambda == NULL is allowed in which case it will not
	  *                be set.  Must have been created by #this->space_c()->create_member()#.
	  *                Must be NULL if m() == 0.
	  * @param lambdaI [out] Pointer to lagrange multipliers for inequalities.
	  *                lambdaI == NULL is allowed in which case it will not
	  *                be set.  Must have been created by #this->space_h()->create_member()#.
	  *                Must be NULL if mI() == 0;
	  * @param nu      [out] Pointer to lagrange multipliers for bounds.
	  *                nu == NULL is allowed in which case it will not
	  *                be set.  Must have been created by #this->space_x()->create_member()#.
	  *                Must be NULL if num_bounded_x() == 0.
	  */
	virtual void get_init_lagrange_mult(
		VectorWithOpMutable*   lambda
		,VectorWithOpMutable*  lambdaI
		,VectorWithOpMutable*  nu
		) const;

	//@}

	/** @name <<std aggr>> stereotype member functions
	  */
	//@{

	/** @name <<std aggr>> members for the objective function value, f.
	 *
	 * The is a reference to a scalar storage location that is updated when
	 * #calc_f(...)# is called or when #f# is updated as a side effect
	 * when another "calc" member funciton is called.
	 */
	//@{

	///
	virtual void set_f(value_type* f);
	///
	virtual value_type* get_f();
	///
	virtual value_type& f();
	///
	virtual const value_type& f() const;

	//@}

	/** @name <<std aggr>> members for the equality constaints vector #c#.
	 *
	 * All of these methods throw exceptions if #m() == 0#.
	 *
	 * The is a reference to a vector storage location that is updated when
	 * #calc_c(...)# is called or when #c# is updated as a side effect
	 * when another "calc" member funciton is called.  The vector #c# must
	 * have been created by #this->space_c()->create_member()#.
	 */
	//@{

	///
	virtual void set_c(VectorWithOpMutable* c);
	///
	virtual VectorWithOpMutable* get_c();
	///
	virtual VectorWithOpMutable& c();
	///
	virtual const VectorWithOp& c() const;

	//@}
 
	/** @name <<std aggr>> members for the general inequality constaints vector #h#.
	 *
	 * All of these methods throw exceptions if #mI() == 0#.
	 *
	 * The is a reference to a vector storage location that is updated when
	 * #calc_h(...)# is called or when #h# is updated as a side effect
	 * when another "calc" member funciton is called.  The vector #h# must
	 * have been created by #this->space_h()->create_member()#.
	 */
	//@{

	///
	virtual void set_h(VectorWithOpMutable* h);
	///
	virtual VectorWithOpMutable* get_h();
	///
	virtual VectorWithOpMutable& h();
	///
	virtual const VectorWithOp& h() const;

	//@}

	//@}

	/** @name Calculation Members.
	  *
	  * These members coordinate the calculation of quanaties at various points #x#.
	  * The #set_info(info)#, where #info# = #f#,  #c# or #h#, are used to set references to
	  * storage to quanities that are updated when #calc_info(x,newx)# memebers are called.
	  * See the <<std comp>> stereotype definition.
	  *
	  * Preconditions: \begin{itemize}
	  * \item #this->get_info() != NULL# (throw #NoRefSet#)
	  * \end{itemize}
	  *
	  * Postconditions: \begin{itemize}
	  * \item #this->info()# returns a const reference to the value to the updated
	  *		quantity at #x#
	  * \end{itemize}
	  *
	  * @param	x		[in] Point at which to calculate the object function #f#.  This
	  *						must have been created by #this->space_x()->create_member()#.
	  * @param	newx	[in] true: The values of #x# are the same as the last call to a 
	  *						#this->calc_info(x,newx)# member.
	  *						false: The values of #x# where not the same as the last call to a
	  *						#this->calc_info(x,newx)# member.
	  *						Default value = true
	  */
	//@{

	///
	/** Set whether the subclass can update multiple quanities to save recalculations.
	  *
	  * The default implementation does nothing.
	  *
	  * @param	set_mult_calc	[in] true: The subclass is allowed to update multiple quantities if it is
	  *							more efficient to do so.  For example, calling calc_f(...) to update 'f'
	  *							may result in an update of 'c' if it is more efficent to do so.
	  *							false: The subclass is not allowed to update multiple quantities.
	  */
	virtual void set_mult_calc(bool mult_calc) const;
	///
	/** Query whether the NLP is allowed to perform multiple updates.
	 *
	 * The default implementation just returns false.
	 */
	virtual bool mult_calc() const;
	///
	/** Set the scaling of the objective function.
	  */
	virtual void scale_f( value_type scale_f ) = 0;
	///
	/** Get the scaling being used for the objective function.
	  */
	virtual value_type scale_f() const = 0;
	///
	/** Update the value for #f# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then storage reference for #c# and/or #h# may also be changed
	  * but is not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_f(const VectorWithOp& x, bool newx = true) const;
	///
	/** Update the vector for #c# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then storage reference for #f# and/or #h# may also be changed
	  * but is not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_c(const VectorWithOp& x, bool newx = true) const;
	///
	/** Update the vector for #h# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then storage reference for #f# and/or #c# may also be changed
	  * but is not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_h(const VectorWithOp& x, bool newx = true) const;

	//@}

	///
	/** Used by the solver to report the final solution and multipliers.
	  *
	  * Call this function to report the final solution of the
	  * unknows x and the Lagrange multipliers for the
	  * equality constriants #lambda# and the varaible bounds
	  * #nu#.  If either of the multipliers
	  * are not known then you can pass #NULL# in for them.
	  *
	  * The default behavior is to just ignore this.
	  */
	virtual void report_final_solution(
		const VectorWithOp&    x
		,const VectorWithOp*   lambda
		,const VectorWithOp*   lambdaI
		,const VectorWithOp*   nu
		,bool                  optimal
		) const;

	/** @name Objective and constraint function counts.
	  *
	  * These functions can be called to find out how many evaluations
	  * the client requested since initialize() was called.
	  *
	  * These do not include any function evaluations that may have
	  * been used for finite difference evaluations or anything
	  * like that.
	  *
	  * Also, if the client calls calc_info( x, false ) several times
	  * then the count will be incremented for each call even though
	  * the actual quantity may not actually be
	  * recalculated each time.  This information is not known here
	  * in this base class but the subclasses can overide this behavior
	  * if they want to.
	  */
	//@{

	///
	/** Gives the number of object function f(x) evaluations called by the solver
	  * since initialize() was called.
	  */
	virtual size_type num_f_evals() const;
	///
	/** Gives the number of constraint function c(x) evaluations called by the solver
	  * since initialize() was called.  Throws exception if #this->m() == 0#.
	  */
	virtual size_type num_c_evals() const;
	///
	/** Gives the number of constraint function h(x) evaluations called by the solver
	  * since initialize() was called.  Throws exception if #this->mI() == 0#.
	  */
	virtual size_type num_h_evals() const;

	//@}

protected:

	///
	/** Struct for objective and constriants (pointer)
	 */
	struct ZeroOrderInfo {
	public:
		///
		ZeroOrderInfo() : f(NULL), c(NULL)
		{}
		///
		ZeroOrderInfo( value_type* f_in, VectorWithOpMutable* c_in, VectorWithOpMutable* h_in )
			: f(f_in), c(c_in), h(h_in)
		{}
		/// Pointer to objective function #f# (may be NULL if not set)
		value_type*           f;
		/// Pointer to constraints residule #c# (may be NULL if not set)
		VectorWithOpMutable*  c;
		/// Pointer to constraints residule #h# (may be NULL if not set)
		VectorWithOpMutable*  h;
	}; // end struct ZeroOrderInfo

	/// Return pointer to set quantities
	const ZeroOrderInfo zero_order_info() const;

	/** @name Protected methods to be overridden by subclasses
	 */
	//@{

	///
	/** Overridden to compute f(x) and perhaps c(x) and/or h(x) (if multiple calculaiton = true).
	 *
	 * Preconditions:\begin{itemize}
	 * \item #x.space().is_compatible(*this->space_x())# (throw #IncompatibleType#)
	 * \item #zero_order_info.f != NULL# (throw #std::invalid_argument#)
	 * \end{itemize}
	 *
	 * Postconditions:\begin{itemize}
	 * \item #*zero_order_info.f# is updated to f(x).
	 * \end{itemize}
	 *
	 * @param x       [in]  Unknown vector (size n).
	 * @param newx    [in]  True if is a new point.
	 * @param zero_order_info
	 *                [out] Pointers to f, c and h.  If #zero_order_info.c != NULL# and #this->mult_calc() == true#
	 *                      then it is allowed that #*zero_order_info.c# will be set on output.
	 *                      If #zero_order_info.h != NULL# and #this->mult_calc() == true#
	 *                      then it is allowed that #*zero_order_info.h# will be set on output.
	 */
	virtual void imp_calc_f(const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const = 0;
	///
	/** Overridden to compute c(x) and perhaps f(x) and/or h(x) (if multiple calculaiton = true).
	 *
	 * Preconditions:\begin{itemize}
	 * \item #x.space().is_compatible(*this->space_x())# (throw #IncompatibleType#)
	 * \item #zero_order_info.c != NULL# (throw #std::invalid_argument#)
	 * \end{itemize}
	 *
	 * Postconditions:\begin{itemize}
	 * \item #*zero_order_info.c# is updated to c(x).
	 * \end{itemize}
	 *
	 * @param x       [in]  Unknown vector (size n).
	 * @param newx    [in]  True if is a new point.
	 * @param zero_order_info
	 *                [out] Pointers to f, c and h.  If #zero_order_info.f != NULL# and #this->mult_calc() == true#
	 *                      then it is allowed that #*zero_order_info.f# will be set on output.
	 *                      If #zero_order_info.h != NULL# and #this->mult_calc() == true#
	 *                      then it is allowed that #*zero_order_info.h# will be set on output.
	 */
	virtual void imp_calc_c(const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const = 0;
	///
	/** Overridden to compute h(x) and perhaps f(x) and/or c(x) (if multiple calculaiton = true).
	 *
	 * Preconditions:\begin{itemize}
	 * \item #x.space().is_compatible(*this->space_x())# (throw #IncompatibleType#)
	 * \item #zero_order_info.h != NULL# (throw #std::invalid_argument#)
	 * \end{itemize}
	 *
	 * Postconditions:\begin{itemize}
	 * \item #*zero_order_info.h# is updated to h(x).
	 * \end{itemize}
	 *
	 * @param x       [in]  Unknown vector (size n).
	 * @param newx    [in]  True if is a new point.
	 * @param zero_order_info
	 *                [out] Pointers to f, c and h.  If #zero_order_info.f != NULL# and #this->mult_calc() == true#
	 *                      then it is allowed that #*zero_order_info.f# will be set on output.
	 *                      If #zero_order_info.c != NULL# and #this->mult_calc() == true#
	 *                      then it is allowed that #*zero_order_info.c# will be set on output.
	 */
	virtual void imp_calc_h(const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const = 0;

	//@}

	/// Assert referece has been set for a quanity
	template<class T>
	void assert_ref_set(T* p, std::string info) const {
		StandardCompositionRelationshipsPack::assert_role_name_set(p, false, info);
	}

private:
	mutable ZeroOrderInfo           first_order_info_;
	mutable size_type				num_f_evals_;
	mutable size_type				num_c_evals_;
	mutable size_type				num_h_evals_;
	
};	// end class NLP

// /////////////////
// Inline members

inline
const NLP::ZeroOrderInfo NLP::zero_order_info() const
{
	return first_order_info_;
}

}	// end namespace NLPInterfacePack 

#endif // NLP_H
