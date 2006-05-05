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

#include <iostream>

#include "IterationPack_TestIterationPack.hpp"
#include "TestingHelperPack_update_success.hpp"

int main() {
  using TestingHelperPack::update_success;
  using namespace IterationPack::TestingPack;

  std::ostream* out = &std::cout;

  bool success = true;
  update_success( TestIterationPack(out), &success );

  if( success )
    std::cerr << "IterationPack seems to check out!\n";
  else
    std::cerr << "Oops! At least one of the tests in IterationPack failed!\n";
    
  return success == true ? 0 : -1;
}
