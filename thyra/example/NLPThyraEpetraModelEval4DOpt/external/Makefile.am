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

install-data-local:
	$(mkinstalldirs) $(exec_prefix)/example/moocho/thyra/NLPThyraEpetraModelEval4DOpt ; \
	$(INSTALL_DATA) $(srcdir)/../NLPThyraEpetraModelEval4DOptMain.cpp $(exec_prefix)/example/moocho/thyra/NLPThyraEpetraModelEval4DOpt/. ; \
	$(INSTALL_DATA) $(top_srcdir)/../epetraext/example/model_evaluator/4dopt/EpetraModelEval4DOpt.* \
           $(exec_prefix)/example/moocho/thyra/NLPThyraEpetraModelEval4DOpt/. ; \
	$(INSTALL_DATA) $(top_srcdir)/src/MoochoPack/configurations/sample_option_files/Moocho.sample.opt \
           $(exec_prefix)/example/moocho/thyra/NLPThyraEpetraModelEval4DOpt/input_files/. ; \
	$(INSTALL_DATA) Makefile.external.stub $(exec_prefix)/example/moocho/thyra/NLPThyraEpetraModelEval4DOpt/Makefile

sample-makefile:
	cp Makefile.external.stub $(srcdir)/Makefile.external.sample
