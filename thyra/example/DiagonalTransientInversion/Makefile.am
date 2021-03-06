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


all-local: $(DIAGONAL_INPUT_FILES)
	@echo ""
	@echo "moocho/thyra/example/DiagonalTransientInversion built successfully."
	@echo ""

#
# Programs
#

noinst_PROGRAMS = \
  DiagonalTransientInversion

DiagonalTransientInversion_SOURCES = \
  $(common_sources) \
  $(srcdir)/DiagonalTransientInversionMain.cpp

DiagonalTransientInversion_DEPENDENCIES = $(common_dependencies)

DiagonalTransientInversion_LDADD = $(common_ldadd)

#
# Input files
#

DIAGONAL_INPUT_FILES = \
  MoochoDiagTransInvSmall.opt \
  diagTransInvSmall_1.xml \
  MoochoDiagTransInvMedium.opt \
  diagTransInvMedium.xml \
  MoochoDiagTransInvBdfVariable.opt \
  diagTransInvBdfVariable.xml

MoochoDiagTransInvSmall.opt: $(srcdir)/_MoochoDiagTransInvSmall.opt
	cp $< $@

diagTransInvSmall_1.xml: $(srcdir)/_diagTransInvSmall_1.xml
	cp $< $@

MoochoDiagTransInvMedium.opt: $(srcdir)/_MoochoDiagTransInvMedium.opt
	cp $< $@

diagTransInvMedium.xml: $(srcdir)/_diagTransInvMedium.xml
	cp $< $@

MoochoDiagTransInvBdfVariable.opt: $(srcdir)/_MoochoDiagTransInvBdfVariable.opt
	cp $< $@

diagTransInvBdfVariable.xml: $(srcdir)/_diagTransInvBdfVariable.xml
	cp $< $@

#
# General options
#

include $(top_srcdir)/MakefileExeStub.mak

CPPFLAGS += \
  -I$(top_srcdir)/../epetraext/example/model_evaluator/DiagonalTransient

common_sources = \
  $(top_srcdir)/../epetraext/example/model_evaluator/DiagonalTransient/EpetraExt_DiagonalTransientModel.cpp

#
# All extra files that have to be added to the tarball!
#

EXTRA_DIST = \
  _MoochoDiagTransInvSmall.opt \
  _diagTransInvSmall_1.xml \
  _MoochoDiagTransInvMedium.opt \
  _diagTransInvMedium.xml \
  _MoochoDiagTransInvBdfVariable.opt \
  _diagTransInvBdfVariable.xml
