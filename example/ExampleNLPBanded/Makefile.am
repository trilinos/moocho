# @HEADER
# ************************************************************************
# 
# Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
#                  Copyright (2003) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
# 
# ************************************************************************
# @HEADER

#
# Build the external example for users to see how to do this.
#

if USING_EXPORT_MAKEFILES
EXTERNAL_MAKEFILE_EXAMPLE = external
endif

SUBDIRS = $(EXTERNAL_MAKEFILE_EXAMPLE)

#
# Build the actual test program in the build directory so that
# we can test this.
#

noinst_PROGRAMS = ExampleNLPBanded

ExampleNLPBanded_SOURCES = \
  $(srcdir)/ExampleNLPBandedMain.cpp\
  $(srcdir)/NLPInterfacePack_ExampleNLPBanded.cpp
ExampleNLPBanded_DEPENDENCIES = $(common_dependencies)
ExampleNLPBanded_LDADD = $(common_ldadd)

include $(top_srcdir)/MakefileExeStub.mak

EXTRA_DIST = $(srcdir)/NLPInterfacePack_ExampleNLPBanded.hpp

