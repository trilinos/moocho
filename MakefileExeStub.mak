# ************************************************************************
# 
# Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
#
#   Copyright (2003) Sandia Corporation
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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER

# Note: This file gets included by automake and these statements are
# interpreted and imported into the create makefiles.  This is very nice
# indeed!

EXEEXT = .exe

@ENABLE_MOOCHO_THYRA_TRUE@MOOCHO_THYRA_LIB = $(top_builddir)/thyra/src/libmoochothyra.a

common_dependencies = \
  $(top_builddir)/src/libmoocho.a \
  $(MOOCHO_THYRA_LIB)

include $(top_builddir)/Makefile.export.moocho

_EXPORT_INCLUDES = $(MOOCHO_INCLUDES) $(RYTHMOS_INCLUDES) $(STRATIMIKOS_INCLUDES) $(EPETRAEXT_INCLUDES)
_EXPORT_LIBS = $(MOOCHO_LIBS) $(RYTHMOS_LIBS) $(STRATIMIKOS_LIBS) $(EPETRAEXT_LIBS)

@USING_GNUMAKE_TRUE@EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(_EXPORT_INCLUDES))
@USING_GNUMAKE_TRUE@EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(_EXPORT_LIBS))
@USING_GNUMAKE_FALSE@EXPORT_INCLUDES = $(_EXPORT_INCLUDES)
@USING_GNUMAKE_FALSE@EXPORT_LIBS = $(_EXPORT_LIBS)

AM_CPPFLAGS = $(EXPORT_INCLUDES)

common_ldadd = $(EXPORT_LIBS)

# This is already added part of MOOCHO_LIBS and therefore automake does not need to add this again!
LIBS =
