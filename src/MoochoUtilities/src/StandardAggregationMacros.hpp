// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef STANDARD_AGGRAGATION_MACROS_H
#define STANDARD_AGGRAGATION_MACROS_H

#include "StandardCompositionRelationshipsPack.hpp"

/** \brief \defgroup StandardAggregationMacros_grp Macros that add <<std aggr>> members for an association.
 * \ingroup Misc_grp
 *
 * For example, if you want to include a <<std aggr>> association
 * with an object of type MyClass of the name my_object you
 * would include the macro in the public section of YourClass
 * declaration as follows:
 *
 \verbatim
  class YourClass {
  public:
    STANDARD_AGGREGATION_MEMBERS( MyClass, my_object )
  };
 \endverbatim
 *
 * Note that the macro addes the private member #TYPE* NAME_#
 * to the class declaration and therefore the member NAME_ is
 * available for direct access (in a constructor for example).
 *
 * In order to have a const only association use:
 \verbatim
  class YourClass {
  public:
    STANDARD_CONST_AGGREGATION_MEMBERS( MyClass, my_object )
  };
 \endverbatim
*/
//@{

/// Insert class members for a non-const association
#define STANDARD_AGGREGATION_MEMBERS( TYPE, NAME )						\
public:																	\
  void set_ ## NAME ( TYPE* NAME )									\
  {	NAME ## _ = NAME; }												\
  TYPE* get_ ## NAME()												\
  {	return NAME ## _; }												\
  const TYPE* get_ ## NAME() const									\
  {	return NAME ## _; }												\
  TYPE& NAME()														\
  {																	\
    return StandardCompositionRelationshipsPack::role_name(			\
      NAME ## _, false, " ## NAME ## " );							\
  }																	\
  const TYPE& NAME() const											\
  {																	\
    return StandardCompositionRelationshipsPack::role_name(			\
      NAME ## _, false, " ## NAME ## " );							\
  }																	\
private:																\
  TYPE* NAME ## _;													\
public:

/// Insert class members for a constant association.
#define STANDARD_CONST_AGGREGATION_MEMBERS( TYPE, NAME )				\
public:																	\
  void set_ ## NAME ( const TYPE* NAME )								\
  {	NAME ## _ = NAME; }												\
  const TYPE* get_ ## NAME() const									\
  {	return NAME ## _; }												\
  const TYPE& NAME() const											\
  {																	\
    return StandardCompositionRelationshipsPack::const_role_name(	\
      NAME ## _, false, " ## NAME ## " );							\
  }																	\
private:																\
  const TYPE* NAME ## _;												\
public:
  
#endif	// STANDARD_AGGRAGATION_MACROS_H
