// ///////////////////////////////////////////////////////////////////
// MA28Solver.h
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
//
// This is a class for encapsulsting the MA28 package as
// and object so that multiple sparse matrix systems can
// be solved at the same time.

#ifdef SPARSE_SOLVER_PACK_USE_MA28

#ifndef SSP_MA28_SOLVER_H
#define SSP_MA28_SOLVER_H

#include "MA28CommonBlockEncap.h"

namespace MA28_Cpp {

// Adobt some of the declarations from MA29_CppDecl
using MA28_CppDecl::f_int;
using MA28_CppDecl::f_logical;
using MA28_CppDecl::f_real;
using MA28_CppDecl::f_dbl_prec;

///
/** MA28 Basic Encapsulation Class.
  *
  * Each object of this class represents a specific MA28 package.
  * Each object encapsulates the common block data for MA28.
  * This class has the same interface functions as MA28 (#ma28ad#
  * , #ma28bd#, and #ma28cd#).  It also has functions for
  * getting and retrieving the common block data from the 
  * MA28 common blocks.
  *
  */
class MA28Solver {
public:

	/// Construct a solver object that is initialized with the default common block data variables.
	MA28Solver();

	/// Construct a solver object that is initialized with the common block data of another solver.
	MA28Solver(const MA28Solver& s);

	// MA28 interface functions

	///
	void ma28ad(const f_int& n, const f_int& nz, f_dbl_prec a[], const f_int& licn
		, f_int irn[], const f_int& lirn, f_int icn[], const f_dbl_prec& u
		, f_int ikeep[], f_int iw[], f_dbl_prec w[], f_int& iflag)
	{	
		set_common_block_data();
		MA28_CppDecl::ma28ad(n,nz,a,licn,irn,lirn,icn,u,ikeep,iw,w,iflag);
		get_common_block_data();
	}

	///
	void ma28bd(const f_int& n, const f_int& nz, f_dbl_prec a[], const f_int& licn
		, const f_int ivect[], const f_int jvect[], const f_int icn[]
		, const f_int ikeep[], f_int iw[], f_dbl_prec w[], f_int& iflag)
	{	
		set_common_block_data();
		MA28_CppDecl::ma28bd(n,nz,a,licn,ivect,jvect,icn,ikeep,iw,w,iflag);
		get_common_block_data();
	}

	///
	void ma28cd(const f_int& n, f_dbl_prec a[], const f_int& licn, const f_int icn[]
		, const f_int ikeep[], f_dbl_prec rhs[], f_dbl_prec w[], const f_int& mtype)
	{	
		set_common_block_data();
		MA28_CppDecl::ma28cd(n,a,licn,icn,ikeep,rhs,w,mtype);
		get_common_block_data();
	}

	// Common block data setting and retrieval functions

	///
	void		lp(f_int lp)
	{	changed_=true; common_blocks_.ma28ed_.lp=lp; }
	///
	f_int		lp()
	{	return common_blocks_.ma28ed_.lp; }
	///
	void		mp(f_int mp)
	{	changed_=true; common_blocks_.ma28ed_.mp=mp; }
	///
	f_int		mp()
	{	return common_blocks_.ma28ed_.mp; }
	///
	void		lblock(f_logical lblock)
	{	changed_=true; common_blocks_.ma28ed_.lblock=lblock; }
	///
	f_logical	lblock()
	{	return common_blocks_.ma28ed_.lblock; }
	///
	void		grow(f_logical grow)
	{	changed_=true; common_blocks_.ma28ed_.grow=grow; }
	///
	f_logical	grow()
	{	return common_blocks_.ma28ed_.grow; }
	///
	void		eps(f_dbl_prec eps)
	{	changed_=true; common_blocks_.ma28fd_.eps=eps; }
	///
	f_dbl_prec	eps()
	{	return common_blocks_.ma28fd_.eps; }
	///
	void		rmin(f_dbl_prec rmin)
	{	changed_=true; common_blocks_.ma28fd_.rmin=rmin; }
	///
	f_dbl_prec	rmin()
	{	return common_blocks_.ma28fd_.rmin; }
	///
	void		resid(f_dbl_prec resid)
	{	changed_=true; common_blocks_.ma28fd_.resid=resid; }
	///
	f_dbl_prec	resid()
	{	return common_blocks_.ma28fd_.resid; }
	///
	void		irncp(f_int irncp)
	{	changed_=true; common_blocks_.ma28fd_.irncp=irncp; }
	///
	f_int		irncp()
	{	return common_blocks_.ma28fd_.irncp; }
	///
	void		icncp(f_int icncp)
	{	changed_=true; common_blocks_.ma28fd_.icncp=icncp; }
	///
	f_int		icncp()
	{	return common_blocks_.ma28fd_.icncp; }
	///
	void		minirn(f_int minirn)
	{	changed_=true; common_blocks_.ma28fd_.minirn=minirn; }
	///
	f_int		minirn()
	{	return common_blocks_.ma28fd_.minirn; }
	///
	void		minicn(f_int minicn)
	{	changed_=true; common_blocks_.ma28fd_.minicn=minicn; }
	///
	f_int		minicn()
	{	return common_blocks_.ma28fd_.minicn; }
	///
	void		irank(f_int irank)
	{	changed_=true; common_blocks_.ma28fd_.irank=irank; }
	///
	f_int		irank()
	{	return common_blocks_.ma28fd_.irank; }
	///
	void		abort1(f_logical abort1)
	{	changed_=true; common_blocks_.ma28fd_.abort1=abort1; }
	///
	f_logical	abort1()
	{	return common_blocks_.ma28fd_.abort1; }
	///
	void		abort2(f_logical abort2)
	{	changed_=true; common_blocks_.ma28fd_.abort2=abort2; }
	///
	f_logical	abort2()
	{	return common_blocks_.ma28fd_.abort2; }
	///
	void		idisp(f_int val, f_int i)
	{	changed_=true; common_blocks_.ma28gd_.idisp[i]=val; }
	///
	f_int		idisp(f_int i)
	{	return common_blocks_.ma28gd_.idisp[i]; }
	///
	void		tol(f_dbl_prec tol)
	{	changed_=true; common_blocks_.ma28hd_.tol=tol; }
	///
	f_dbl_prec	tol()
	{	return common_blocks_.ma28hd_.tol; }
	///
	void		themax(f_dbl_prec themax)
	{	changed_=true; common_blocks_.ma28hd_.themax=themax; }
	///
	f_dbl_prec	themax()
	{	return common_blocks_.ma28hd_.themax; }
	///
	void		big(f_dbl_prec big)
	{	changed_=true; common_blocks_.ma28hd_.big=big; }
	///
	f_dbl_prec	big()
	{	return common_blocks_.ma28hd_.big; }
	///
	void		dxmax(f_dbl_prec dxmax)
	{	changed_=true; common_blocks_.ma28hd_.dxmax=dxmax; }
	///
	f_dbl_prec	dxmax()
	{	return common_blocks_.ma28hd_.dxmax; }
	///
	void		errmax(f_dbl_prec errmax)
	{	changed_=true; common_blocks_.ma28hd_.errmax=errmax; }
	///
	f_dbl_prec	errmax()
	{	return common_blocks_.ma28hd_.errmax; }
	///
	void		dres(f_dbl_prec dres)
	{	changed_=true; common_blocks_.ma28hd_.dres=dres; }
	///
	f_dbl_prec	dres()
	{	return common_blocks_.ma28hd_.dres; }
	///
	void		cgce(f_dbl_prec cgce)
	{	changed_=true; common_blocks_.ma28hd_.cgce=cgce; }
	///
	f_dbl_prec	cgce()
	{	return common_blocks_.ma28hd_.cgce; }
	///
	void		ndrop(f_int ndrop)
	{	changed_=true; common_blocks_.ma28hd_.ndrop=ndrop; }
	///
	f_int		ndrop()
	{	return common_blocks_.ma28hd_.ndrop; }
	///
	void		maxit(f_int maxit)
	{	changed_=true; common_blocks_.ma28hd_.maxit=maxit; }
	///
	f_int		maxit()
	{	return common_blocks_.ma28hd_.maxit; }
	///
	void		noiter(f_int noiter)
	{	changed_=true; common_blocks_.ma28hd_.noiter=noiter; }
	///
	f_int		noiter()
	{	return common_blocks_.ma28hd_.noiter; }
	///
	void		nsrch(f_int nsrch)
	{	changed_=true; common_blocks_.ma28hd_.nsrch=nsrch; }
	///
	f_int		nsrch()
	{	return common_blocks_.ma28hd_.nsrch; }
	///
	void		istart(f_int istart)
	{	changed_=true; common_blocks_.ma28hd_.istart=istart; }
	///
	f_int		istart()
	{	return common_blocks_.ma28hd_.istart; }
	///
	void		lbig(f_logical lbig)
	{	changed_=true; common_blocks_.ma28hd_.lbig=lbig; }
	///
	f_logical	lbig()
	{	return common_blocks_.ma28hd_.lbig; }

	/// Dump the common block infomation for this solver object.
	void dump_common_blocks(std::ostream& o)
	{	common_blocks_.dump_values(o); }

	/// Copy the state of one solver to another
	MA28Solver& operator=(const MA28Solver& solver)
	{	changed_ = true; common_blocks_ = solver.common_blocks_; return *this; }

	// ///////////////////////////////////
	// Static member functions

	/// Dump the common block infomation for ma28 common blocks
	static void dump_ma28_common_blocks(std::ostream& o)
	{	ma28_common_blocks_.dump_values(o); }

private:

	// ////////////////////////////////////
	// Private member functions

	// Copy the local copy the common block data to MA28 before a MA28 call.
	void set_common_block_data();

	// Retrieve the common block data after a ma28 call.
	void get_common_block_data();
	
	// ///////////////////////////////////
	// Private member data

	// Common block data for this solver object
	 MA28CommonBlockStorage common_blocks_;

	// Flag for if the common bock data has changed
	bool changed_;

	// ///////////////////////////////////
	// Static member data

	// Copies of the default values for the 
	// common block data.
	static MA28CommonBlockStorage default_common_blocks_;

	// References to the MA28 common blocks
	static MA28CommonBlockReferences ma28_common_blocks_;

	// Pointer variable who's purpose it to identify
	// what solver object is the current one.
	static MA28Solver* curr_solver_;
	
};

}	// end namespace MA28_Cpp

#endif // SSP_MA28_SOLVER_H

#endif // SPARSE_SOLVER_PACK_USE_MA28