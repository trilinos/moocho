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
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "StandardCompositionRelationshipsPack.h"
#include "ref_count_ptr.h"

namespace NLPInterfacePack {

///
/** %NLP interface class {abstract}.
 *
 * <b>Overview:</b>
 *
 * This class represents an abstract interface to a general nonlinear
 * programming problem of the form:
 \verbatim

     min    f(x)
     s.t.   c(x) = 0
            hl <= h(x) <= hu
	        xl <= x <= xu
	where:
	        x    <: R^n
	        c(x) <: R^n -> R^m 
	        h(x) <: R^n -> R^mI 
 \endverbatim
 * In the above form, none of the variables are fixed between bounds (strictly
 * xl < xu).  It is allowed however for <tt>m == 0</tt> and/or <tt>mI == 0</tt>
 * for the elimination of more general constriants.  It is also allowed for
 * <tt>n == m</tt> in which case <tt>this</tt> represents a fully determined
 * system of nonlinear equaltions.  In any case, an objective function is always
 * included in the formutation and will impact solution algorithms.
 *
 * Special types of NLPs are identified as: <ol>
 * <li> Fully general %NLP :
 *      <ul><li><tt>( xl != -Inf || xu != +inf ) && ( m  > 0 && mI  > 0 )</tt></ul>
 * <li> General equality only constrained %NLP :
 *      <ul><li><tt>( xl == -Inf && xu == +Inf ) && ( m  > 0 && mI == 0 )</tt></ul>
 * <li> General inequality constrained %NLP :
 *      <ul><li><tt>( xl ==  ??? && xu ==  ??? ) && ( m == 0 && mI  > 0 )</tt></ul>
 * <li> Bound constrained %NLP :
 *      <ul><li><tt>( xl != -Inf || xu != +Inf ) && ( m == 0 && mI == 0 )</tt></ul>
 * <li> Unconstrained %NLP :
 *      <ul><li><tt>( xl == -Inf && xu == +Inf ) && ( m == 0 && mI == 0 )</tt></ul>
 * <li> Nonlinear Equations (NLE) :
 *      <ul><tt>n == m</tt><li></ul>
 * </ol>
 *
 * Note that in (6) above, that this allows for bounds on the variables and general
 * inequality constriants.  If some of the equations in <tt>c(x)</tt> are linearly
 * dependent (but consistent) then it is possible that some of these extra inequalities
 * my be active at the solution but in general they can not be (unless they are degenerate).
 * An optimization algorithm may refuse to solve some of the above problems but this
 * interface allows all of these possibilities.
 *
 * The Lagrangian for this problem is defined by:
 \verbatim

	L = f(x) + lambda' * c(x)
        + lambdaI_l' * ( hl - h(x) ) + lambdaI_u' * ( h(x) - hu )
        + nul' * ( xl - x ) + nuu' * ( x - xu )
 \endverbatim
 * The optimality conditions are given by:
 \verbatim

	del(L,x)      = del(f,x) + del(c,x) * lambda + del(h,x) * lambdaI + nu = 0
	del(L,lambda) = c(x) = 0
	  where:
        lambdaI = lambdaI_u - lambdaI_l
		nu = nuu - nul
		lambdaI_l(j) * ( h(x)(j) - hl(j) ) = 0,    for j = 1...mI
		lambdaI_u(j) * ( h(x)(j) - hu(j) ) = 0,    for j = 1...mI
		nuu(i) * ( x(i) - xu(i) ) = 0,      for i = 1...n
		nuu(i) * ( x(i) - xu(i) ) = 0,      for i = 1...n
 \endverbatim
 * What is unique about this interface is that the vector objects are hidden behind
 * abstact interfaces.  Clients can create vectors from the various vector spaces
 * using the <tt>\ref AbstractLinAlgPack::VectorSpace "VectorSpace"</tt> objects returned from
 * <tt>this->space_x()</tt> (dim \c n), <tt>this->space_c()</tt> (dim \c m) and
 * <tt>this->space_h()</tt> (dim \c mI).
 *
 * <b>Client Usage:</b>
 *
 * Before an %NLP object can be used, the <tt>initialize()</tt> method must be called to
 * make sure that all of the initializations needed for the NLP have been performed.
 * This method also resets counters an other information.  Before calling <tt>initialize()</tt>
 * the client can specify whether the initial point for \c x must be in bounds by calling
 * <tt>force_xinit_in_bounds(bool)</tt>.
 *
 * Smart reference counted pointers to the three vector spaces for \a x, \a c(x) and \a h(x)
 * are returned by the methods <tt>space_x()</tt>, <tt>space_c()</tt> and <tt>space_h()</tt>
 * respectively.  The vector space objects returned by these methods are ment to be more
 * than transient.  In fact, it is expected that these vector space objects should remain
 * valid for the entire run of an NLP algorithm.  Only if the underlying NLP is changed in
 * a fundamental way (i.e. \c n, \c m or \c mI changes) should the vector space objects returned
 * from these function become invalid.  In this case the client must call these methods again
 * to get updated vector space objects.
 *
 * The dimensionality of the NLP is returned by the methods <tt>n()</tt>, <tt>m()</tt>
 * and <tt>>mI()</tt> but they have default implementations based on <tt>space_x()</tt>,
 * <tt>space_c()</tt> and <tt>space_h()</tt> respectively.
 *
 * The number of variables \c x with finite bounds is returned by the method <tt>num_bounded_x()</tt>.
 * If <tt>num_bounded_x() > 0</tt> then the methods <tt>xl()</tt> and <tt>xu()</tt> return references
 * to vector objects reprresenting these bounds.  A lower bound is considered infinite if
 * <tt>xl().get_ele(i) == -infinite_bound()</tt> and an upper bound is considered infinite if
 * <tt>xu().get_ele(i) == +infinite_bound()</tt>.
 *
 * If <tt>mI() > 0</tt>, the methods \c hl() and \c hu() return references to the upper and lower
 * bounds to the general inequality constraints \a h(x).  While it is expected that 
 * <tt>hl().get_ele(j) != -infinite_bound() || hu().get_ele(j) != +infinite_bound</tt> for <tt>j = 1...mI()</tt>
 * this is not required by this interface.  On the other hand it seems silly to define general inequality constriants
 * that are not bounded but there may be some reason to include these that makes things easier for the
 * implementor of the NLP subclass.
 *
 * The initial guess for the unknowns \a x (primal variables) is returned by the method <tt>xinit()</tt>.
 * Vectors containing the initial guesses for the Lagrange multipliers can be obtained by calling the
 * method <tt>get_init_lagrange_mult()</tt>.  
 * 
 * The bread and butter of an %NLP interface is the calculation of the functions that define the objective
 * and constraints and various points \a x and the methods \c calc_f(), \c calc_c() and \c calc_h()
 * accomplish this.  The quantities that these functions update must be set prior by calling the methods
 * \c set_f(), \c set_c() and \c set_h() respectively.  It may seem strange not to pass these quantities
 * directly to the calculation functions but there is a good reason why it is done this way.
 * The reason is that this interface supports the efficient update of mutiple quantities at a given
 * point \a x.  For example, in many %NLPs some of the same terms are shared between the constriants
 * functions and objective function.  Therefore it is more efficient to compute \a f(x), \a c(x) and \a h(x)
 * simultaneously rather than computing them separately.  In order to allow for this possibility,
 * the client can set the desired quantities (i.e. \c set_f(), \c set_c() and \c set_h()) and then
 * calling the method \c set_mult_calc(true) prior to calling \c calc_f(), \c calc_c() and \c calc_h().
 * Another reason for structuring this %NLP interface this way is when automatic differentiation is used
 * to compute derivatives (see \c NLPFirstOrderInfo) the function values are computed for free.
 *
 * Once an optimization algorithm has the solution (or gives up with a suboptimal point), it should
 * report this solution to the %NLP object using the method \c report_final_solution().
 *
 * Finally, the client can get the counts for the number of function evaluations since \c initialize()
 * was called using the methods \c num_f_evals(), \c num_c_evals() and \c num_h_evals().
 * These counts do not include any function evaluations that may have been used internally for finite
 * difference evaluations or anything of that nature.  Also, if the client calls <tt>calc_info(x,false)</tt>
 * (where \c info = \c f, \c c or \c h) several times then the default implementation will increment the count
 * for each call even though the actual quantity may not actually be recalculated each time (i.e. if <tt>newx==false</tt>).
 * This information is not known here in this base class but the subclasses can overide this behavior if desired.
 *
 * <b>Subclass developer's notes:</b>
 *
 * The calcuation methods \c calc_f(), \c calc_c() and \c calc_h() have default implementations in this base
 * class that should meet the needs of all the subclasses.  These method implementations call the protected
 * pure virtual methods \c imp_calc_f(), \c imp_calc_c() and \c imp_calc_h() to compute the actual quantities.
 * Subclasses must override these methods (in addition to several methods from the public interface).
 * Pointers to the quantities to be updated are passed to these methods to the subclasses in the form of
 * an aggregate \c ZeroOrderInfo object that is returned from \c the protected method zero_order_info().
 * This ensures that the only interaction between an NLP base object and its subclass objects is through
 * member functions and never through pulic or protected data members.
 *
 * <A NAME="must_override"></A>
 * The following methods must be overridden by a subclass in order to create a concrete NLP object:
 * \c force_xinit_in_bounds(bool), \c force_xinit_in_bounds(), \c is_initialized(),
 * \c space_x(), \c space_c(), \c space_h(), \c num_bounded_x(), \c xl(), xu(), \c hl(), \c hu(),
 * \c xinit(), \c scale_f(value_type), \c scale_f().
 *
 * <A NAME="should_override"></A>
 * The following methods should be overridden by most subclasses but do not have to be: \c initialize(),
 * \c n(), \c m(), \c mI(), \c get_init_lagrange_mult(), \c set_mult_calc(bool), \c mult_calc(),
 * \c report_final_solution().
 *
 * The following methods should never have to be overridden by most subclasses except in some very
 * strange situations: \c set_f(), \c get_f(), \c f(), \c set_c(), \c get_c(), \c c(),
 * \c set_h(), \c get_h(), \c h(), \c calc_f(), \c calc_c(), \c calc_h(), \c num_f_evals(),
 * \c num_c_evals(), \c num_h_evals().
 *
 * <b>Additional notes:</b>
 *
 * This is where the bounds on the variables can play a very
 * critical part.  It is desirable for the functions \a f(x), \a c(x) and \a h(x) to be defined and
 * relatively well behaved (i.e. smooth, continous, differentiable etc.) in the region <tt>xl <= x <= xu</tt>.
 * While this is not always possible, an NLP can often be reformulated to have this properly.  For example, suppose
 * there are constraints has the form:
 \verbatim

 log(x(1) - x(5)) - 4 == 0
 x(1) - x(5) >= 1e-8
 \endverbatim
 * where \c x(1) and x(5) are unbounded.  It is clear that if <tt>x(1) < x(5)</tt> that this constraint will be
 * undefined and will return \c NaN on most computers (if you are lucky).  The constraint <tt>x(1) < x(5)</tt> is
 * very hard to enforce at every iteration in an NLP solver so that this will not happen   A better approach would be
 * to add an extra varaible (say \c x(51) for an %NLP with <tt>n == 50</tt>) and add an extra constraint:
 \verbatim

 log(x(51)) == 0
 x(51) - (x(1) - x(5)) == 0
 x(51) >= 1e-8
 \endverbatim
 * In the above expanded formulation the simple bound <tt>x(51) >= 1e-8</tt> is easy to inforce and these undefined
 * regions can be avoided.  While the property that \a f(x), \a c(x) and \a h(x) being bounded for all
 * <tt>x <: { x | xl <= x <= x}</tt> is a desireable properly, this is not required by this interface.  As a result
 * the client should be prepaired to deal with return values of \c NaN or \c Inf for \c f, \c c and \c h.
 */
class NLP {
public:

	typedef AbstractLinAlgPack::VectorWithOp         VectorWithOp;         // doxygen likes typedef?
	typedef AbstractLinAlgPack::VectorWithOpMutable  VectorWithOpMutable;  // doxygen likes typedef?
	
	///
	typedef ReferenceCountingPack::ref_count_ptr<const VectorSpace>  vec_space_ptr_t;

	/** @name exceptions */
	//@{

	/// Thrown if any member functions are called before initialize() has been called.
	class UnInitialized : public std::logic_error
	{public: UnInitialized(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown from <tt>initialize()</tt> if some logical error occured
	class InvalidInitialization : public std::logic_error
	{public: InvalidInitialization(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if an incompatible object is used
	class IncompatibleType : public std::logic_error
	{public: IncompatibleType(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown some bounds do not existe
	class NoBounds : public std::logic_error
	{public: NoBounds(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}

	/// Value for an infinite bound.
	static value_type infinite_bound();

	/** @name Constructors, Destructor */
	//@{

	/// Initialize to no reference set to calculation quanities
	NLP();
	/// Destructor that cleans all the memory it owns
	virtual ~NLP();

	//@}
	
	/** @name NLP initialization */
	//@{

	///
	/** Set if the initial point must be within the bounds.
	 *
	 * This method must be called before <tt>this->initialize()</tt> is called.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->is_initialized() == false</tt>
	 * </ul>
	 */
	virtual void force_xinit_in_bounds(bool force_xinit_in_bounds) = 0;
	///
	/** Returns if the initial point must be within the bounds.
	  */
	virtual bool force_xinit_in_bounds() const = 0;
	///
	/** Initialize the NLP for it is used.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt>
	 * <li> [<tt>this->force_xinit_in_bounds()==true && this->num_bounded_x() > 0</tt>]
	 *      <tt>this->xl() <= this->xinit() <= this->xu()</tt>
	 * <li> <tt>this->num_f_evals() == 0</tt>
	 * <li> [<tt>this->m() > 0</tt>] <tt>this->num_c_evals() == 0</tt>
	 * <li> [<tt>this->mI() > 0</tt>] <tt>this->num_h_evals() == 0</tt>
	 * </ul>
	 *
	 * Note that subclasses must call this function to reset what needs to be
	 * reset in this base object.
	 */
	virtual void initialize();
	///
	/** Return if <tt>this</tt> is initialized.
	  */
	virtual bool is_initialized() const = 0;

	//@}

	/** @name Dimensionality. */
	//@{

	///
	/** Return the number of variables.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Default implementation returns <tt>this-space_x()->dim()</tt>.
	 */
	virtual size_type n() const;
	///
	/** Return the number of general equality constraints.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Default implementation returns <tt>this->space_c().get() != NULL ? this-space_c()->dim() : 0</tt>.
	 */
	virtual size_type m() const;
	///
	/** Return the number of general inequality constraints.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Default implementation returns <tt>this->space_h().get() != NULL ? this-space_h()->dim() : 0</tt>.
	 */
	virtual size_type mI() const;

	//@}

	/** @name Vector Space objects. */
	//@{

	///
	/** Vector space object for unknown variables x (dimension n).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get() != NULL</tt>
	 * </ul>
	 */
	virtual vec_space_ptr_t space_x() const = 0;
	///
	/** Vector space object for general equality constraints c(x) (dimension m).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->m() > 0</tt>] <tt>return.get() != NULL</tt>
	 * <li> [<tt>this->m() == 0</tt>] <tt>return.get() == NULL</tt>
	 * </ul>
	 */
	virtual vec_space_ptr_t space_c() const = 0;

	///
	/** Vector space object for general inequality constraints c(x) (dimension mI).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->mI() > 0</tt>] <tt>return.get() != NULL</tt>
	 * <li> [<tt>this->mI() == 0</tt>] <tt>return.get() == NULL</tt>
	 * </ul>
	 */
	virtual vec_space_ptr_t space_h() const = 0;

	//@}

	/** @name Bounds on the unknown variables x. */
	//@{

	///
	/** Returns the number of variables in <tt>x(i)</tt> for which <tt>xl(i)> -infinite_bound()</tt>
	 * or <tt>xu(i) < +infinite_bound()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual size_type num_bounded_x() const = 0;
	///
	/** Returns the lower bounds on the variables <tt>x</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->num_bounded_x() > 0</tt> (throw <tt>NoBounds</tt>)
	 * </ul>
	 *
	 * Any bounds that are non-existant will return <tt>this->xl().get_ele(i) == -NLP::infinite_bound()</tt>.
	 */
	virtual const VectorWithOp& xl() const = 0;
	///
	/** Returns a reference to the vector of upper bounds on the variables <tt>x</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->num_bounded_x() > 0</tt> (throw <tt>NoBounds</tt>)
	 * </ul>
	 *
	 * Any bounds that are non-existant will return <tt>this->xu().get_ele(i) == +NLP::infinite_bound()</tt>.
	 */
	virtual const VectorWithOp& xu() const = 0;

	//@}

	/** @name Bounds on the general inequality constraints h(x). */
	//@{

	///
	/** Returns the lower bounds on the inequality constraints <tt>h(x)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->mI() > 0</tt> (throw <tt>NoBounds</tt>).
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.space().is_compatible(*this->space_h()) == true</tt>
	 * </ul>
	 *
	 * Any bounds that are non-existant will return <tt>this->hl().get_ele(j) == -NLP::infinite_bound()</tt>.
	 */
	virtual const VectorWithOp& hl() const = 0;
	///
	/** Returns upper bounds on the inequality constraints <tt>h(x)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->mI() > 0</tt> (throw <tt>NoBounds</tt>).
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.space().is_compatible(*this->space_h()) == true</tt>
	 * </ul>
	 *
	 * Any bounds that are non-existant will return <tt>this->hu().get_ele(j) == +NLP::infinite_bound()</tt>.
	 */
	virtual const VectorWithOp& hu() const = 0;

	//@}

	/** @name Initial guess of NLP solution */
	//@{

	///
	/** Returns a reference to the vector of the initial guess for the solution <tt>x</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.space().is_compatible(*this->space_x()) == true)</tt>
	 * </ul>
	 */
	virtual const VectorWithOp& xinit() const = 0;
	///
	/** Get the initial value of the Lagrange multipliers lambda.
	 *
	 * By default this function just sets them to zero.
	 *
	 *
	 * @param lambda  [out] Pointer to lagrange multipliers for equalities.
	 *                lambda == NULL is allowed in which case it will not
	 *                be set.  Must have been created by <tt>this->space_c()->create_member()</tt>.
	 *                Must be NULL if m() == 0.
	 * @param lambdaI [out] Pointer to lagrange multipliers for inequalities.
	 *                lambdaI == NULL is allowed in which case it will not
	 *                be set.  Must have been created by <tt>this->space_h()->create_member()</tt>.
	 *                Must be NULL if mI() == 0;
	 * @param nu      [out] Pointer to lagrange multipliers for bounds.
	 *                nu == NULL is allowed in which case it will not
	 *                be set.  Must have been created by <tt>this->space_x()->create_member()</tt>.
	 *                Must be NULL if num_bounded_x() == 0.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>lambda != NULL</tt>] <tt>lambda->size() == this->n()</tt>
	 * <li> [<tt>nu != NULL</tt>] <tt>nu->size() == this->n()</tt>
	 * </ul>
	 */
	virtual void get_init_lagrange_mult(
		VectorWithOpMutable*   lambda
		,VectorWithOpMutable*  lambdaI
		,VectorWithOpMutable*  nu
		) const;

	//@}

	/** @name <<std aggr>> members for the objective function value f(x). */
	//@{

	///
	/** Set a pointer to an value to be updated when <tt>this->calc_f()</tt> is called.
	 *
	 * @param  f  [in] Pointer to objective function value.  May be \c NULL.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get_f() == f</tt>
	 * </ul>
	 */
	virtual void set_f(value_type* f);
	///
	/** Return pointer passed to <tt>this->set_f()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual value_type* get_f();
	///
	/** Returns non-<tt>const</tt> <tt>*this->get_f()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_f() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual value_type& f();
	///
	/** Returns <tt>const</tt> <tt>*this->get_f()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_f() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
 	virtual const value_type& f() const;

	//@}

	/** @name <<std aggr>> members for the residual of the general equality constriants c(x). */
	//@{

	///
	/** Set a pointer to a vector to be updated when <tt>this->calc_c()</tt> is called.
	 *
	 * @param  c  [in] Pointer to constraint residual vector.  May be \c NULL.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> [<tt>c != NULL</tt>] <tt>c->space().is_compatible(*this->space_c()) == true</tt>
	 *      (throw <tt>VectorSpaceBase::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get_c() == c</tt>
	 * </ul>
	 */
	virtual void set_c(VectorWithOpMutable* c);
	///
	/** Return pointer passed to <tt>this->set_c()</tt>.
	 */
	virtual VectorWithOpMutable* get_c();
	///
	/** Returns non-<tt>const</tt> <tt>*this->get_c()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_c() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual VectorWithOpMutable& c();
	///
	/** Returns <tt>const</tt> <tt>*this->get_c()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_c() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual const VectorWithOp& c() const;

	//@}

	/** @name <<std aggr>> members for the residual of the general inequality constriants h(x). */
	//@{

	///
	/** Set a pointer to a vector to be updated when <tt>this->calc_h()</tt> is called.
	 *
	 * @param  h  [in] Pointer to constraint residual vector.  May be \c NULL.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> [<tt>h != NULL</tt>] <tt>h->space().is_compatible(*this->space_h()) == true</tt>
	 *      (throw <tt>VectorSpaceBase::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get_h() == h</tt>
	 * </ul>
	 */
	virtual void set_h(VectorWithOpMutable* h);
	///
	/** Return pointer passed to <tt>this->set_h()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual VectorWithOpMutable* get_h();
	///
	/** Returns non-<tt>const</tt> <tt>*this->get_h()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_h() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual VectorWithOpMutable& h();
	///
	/** Returns <tt>const</tt> <tt>*this->get_h()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_h() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual const VectorWithOp& h() const;

	//@}

	/** @name Calculation Members. */
	//@{

	///
	/** Set whether the subclass can update multiple quanities to save recalculations.
	 *
	 * The default implementation does nothing.
	 *
	 * @param  set_mult_calc
	 *                [in] If \c true the subclass is allowed to update multiple quantities if it is
	 *                more efficient to do so.  For example, calling calc_f(...) to update 'f'
	 *                may result in an update of 'c' if it is more efficent to do so.
	 *                If \c false,  the subclass is not allowed to update multiple quantities.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual void set_mult_calc(bool mult_calc) const;
	///
	/** Query whether the NLP is allowed to perform multiple updates.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The default implementation just returns false.
	 */
	virtual bool mult_calc() const;
	///
	/** Set the scaling of the objective function.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->scale_f() == true</tt>
	 * </ul>
	 */
	virtual void scale_f( value_type scale_f ) = 0;
	///
	/** Return the scaling being used for the objective function.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual value_type scale_f() const = 0;
	///
	/** Update the value for the objective <tt>f</tt> at the point <tt>x</tt> and put it in the stored reference.
	 *
	 * @param  x     [in] Point at which to calculate the object function <tt>f</tt>.
	 * @param  newx  [in] (default \c true) If \c true, the values in \c x are the same as
	 *               the last call to a <tt>this->calc_*(x,newx)</tt> member.
	 *               If \c false, the values in \c x are not the same as the last call to a
	 *               <tt>this->calc_*(x,newx)</tt> member.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpaceBase::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>this->get_f() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->f()</tt> is updated to \a f(x)
	 * </ul>
	 *
	 * If <tt>set_mult_calc(true)</tt> was called then storage reference for <tt>c</tt> and/or <tt>h</tt> may also be changed
	 * but is not guarentied to be.  But no other quanities from possible subclasses are allowed
	 * to be updated as a side effect.
	 */ 
	virtual void calc_f(const VectorWithOp& x, bool newx = true) const;
	///
	/** Update the constraint residual vector for <tt>c</tt> at the point <tt>x</tt> and put it in the stored reference.
	 *
	 * @param  x     [in] Point at which to calculate residual to the equality constraints <tt>c</tt>.
	 * @param  newx  [in] (default \c true) If \c true, the values in \c x are the same as
	 *               the last call to a <tt>this->calc_*(x,newx)</tt> member.
	 *               If \c false, the values in \c x are not the same as the last call to a
	 *               <tt>this->calc_*(x,newx)</tt> member.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpaceBase::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>this->get_c() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->c()</tt> is updated to \a c(x)
	 * </ul>
	 *
	 * If <tt>set_mult_calc(true)</tt> was called then storage reference for <tt>f</tt> and/or <tt>h</tt> may also be changed
	 * but is not guarentied to be.  But no other quanities from possible subclasses are allowed
	 * to be updated as a side effect.
	 */ 
	virtual void calc_c(const VectorWithOp& x, bool newx = true) const;
	///
	/** Update the vector for <tt>h</tt> at the point <tt>x</tt> and put it in the stored reference.
	 *
	 * @param  x     [in] Point at which to calculate residual to the inequality constraints <tt>h</tt>.
	 * @param  newx  [in] (default \c true) If \c true, the values in \c x are the same as
	 *               the last call to a <tt>this->calc_*(x,newx)</tt> member.
	 *               If \c false, the values in \c x are not the same as the last call to a
	 *               <tt>this->calc_*(x,newx)</tt> member.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpaceBase::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>this->get_h() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->h()</tt> is updated to \a h(x)
	 * </ul>
	 *
	 * If <tt>set_mult_calc(true)</tt> was called then storage reference for <tt>f</tt> and/or <tt>c</tt> may also be changed
	 * but is not guarentied to be.  But no other quanities from possible subclasses are allowed
	 * to be updated as a side effect.
	 */ 
	virtual void calc_h(const VectorWithOp& x, bool newx = true) const;

	//@}

	/** @name Report final solution */
	//@{

	///
	/** Used by the solver to report the final solution and multipliers.
	 *
	 * Call this function to report the final solution of the
	 * unknows x and the Lagrange multipliers for the
	 * equality constriants <tt>lambda</tt> and the varaible bounds
	 * <tt>nu</tt>.  If any of the Lagrange multipliers
	 * are not known then you can pass <tt>NULL</tt> in for them.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The default behavior is to just ignore this.
	 */
	virtual void report_final_solution(
		const VectorWithOp&    x
		,const VectorWithOp*   lambda
		,const VectorWithOp*   lambdaI
		,const VectorWithOp*   nu
		,bool                  is_optimal
		) const;

	//@}

	/** @name Objective and constraint function evaluation counts. */
	//@{

	///
	/** Gives the number of object function f(x) evaluations called by the solver
	 * since initialize() was called.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual size_type num_f_evals() const;
	///
	/** Gives the number of constraint function c(x) evaluations called by the solver
	 * since initialize() was called.  Throws exception if <tt>this->m() == 0</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual size_type num_c_evals() const;
	///
	/** Gives the number of constraint function h(x) evaluations called by the solver
	 * since initialize() was called.  Throws exception if <tt>this->mI() == 0</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual size_type num_h_evals() const;

	//@}

protected:

	///
	/** Struct for objective and constriants (pointer).
	 *
	 * Objects of this type are passed on to subclasses and contain pointers to
	 * quantities to be updated.
	 */
	struct ZeroOrderInfo {
	public:
		///
        ZeroOrderInfo() : f(NULL), c(NULL), h(NULL)
		{}
		///
		ZeroOrderInfo( value_type* f_in, VectorWithOpMutable* c_in, VectorWithOpMutable* h_in )
			: f(f_in), c(c_in), h(h_in)
		{}
		/// Pointer to objective function <tt>f</tt> (may be NULL if not set)
		value_type*           f;
		/// Pointer to constraints residule <tt>c</tt> (may be NULL if not set)
		VectorWithOpMutable*  c;
		/// Pointer to constraints residule <tt>h</tt> (may be NULL if not set)
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
	 * Preconditions:<ul>
	 * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
	 * <li> <tt>zero_order_info.f != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*zero_order_info.f</tt> is updated to \a f(x).
	 * </ul>
	 *
	 * @param x       [in]  Unknown vector (size n).
	 * @param newx    [in]  True if is a new point.
	 * @param zero_order_info
	 *                [out] Pointers to \c f, \c c and \c h.
	 *                On output, <tt>*zero_order_info.f</tt> is updated to \a f(x)
	 *                If <tt>this->mult_calc() == true</tt> then
	 *                any of the other quantities pointed to in \c zero_order_info may be set on
	 *                output, but are not guaranteed to be.
	 */
	virtual void imp_calc_f(const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const = 0;
	///
	/** Overridden to compute c(x) and perhaps f(x) and/or h(x) (if multiple calculaiton = true).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
	 * <li> <tt>zero_order_info.c != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*zero_order_info.c</tt> is updated to c(x).
	 * </ul>
	 *
	 * @param x       [in]  Unknown vector (size n).
	 * @param newx    [in]  True if is a new point.
	 * @param zero_order_info
	 *                [out] Pointers to \c f, \c c and \c h.
	 *                On output, <tt>*zero_order_info.c</tt> is updated to \a c(x)
	 *                If <tt>this->mult_calc() == true</tt> then
	 *                any of the other quantities pointed to in \c zero_order_info may be set on
	 *                output, but are not guaranteed to be.
	 */
	virtual void imp_calc_c(const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const = 0;
	///
	/** Overridden to compute h(x) and perhaps f(x) and/or c(x) (if multiple calculaiton = true).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
	 * <li> <tt>zero_order_info.h != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*zero_order_info.h</tt> is updated to h(x).
	 * </ul>
	 *
	 * @param x       [in]  Unknown vector (size n).
	 * @param newx    [in]  True if is a new point.
	 * @param zero_order_info
	 *                [out] Pointers to \c f, \c c and \c h.
	 *                On output, <tt>*zero_order_info.h</tt> is updated to \a h(x)
	 *                If <tt>this->mult_calc() == true</tt> then
	 *                any of the other quantities pointed to in \c zero_order_info may be set on
	 *                output, but are not guaranteed to be.
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
