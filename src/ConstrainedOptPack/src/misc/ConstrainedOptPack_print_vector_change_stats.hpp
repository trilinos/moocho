// ///////////////////////////////////////////////////////
// print_vector_change_stats.hpp
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

#ifndef PRINT_VECTOR_CHANGE_STATS_H
#define PRINT_VECTOR_CHANGE_STATS_H

#include "ConstrainedOptimizationPackTypes.hpp"

namespace ConstrainedOptimizationPack {

///
/** Compute statistics for change in a vector and output to a stream.
  *
  * Calls the function vector_change_stats(x,d,max_term,max_k,min_term,min_k
	,av_term) then produces the following output to the given stream.
  *
  * max(|d(i)|/(1+|x(i)|)) => |d(max_k)|/(1+|x(max_k)|) = max_term \\
  * min(|d(i)|/(1+|x(i)|)) => |d(min_k)|/(1+|x(min_k)|) = min_term \\ 
  * average(|d(i)|/(1+|x(i)|)) = av_term \\
  *
  * Note that above the names x and d are replaced with their input
  * names #x_name# and #d_name# and max_term, max_k etc.
  * are replaced with their computed values.
  */
void print_vector_change_stats( const VectorSlice& x, const char x_name[]
	, const VectorSlice& d, const char d_name[], std::ostream& out );

}	// end namespace ConstrainedOptimizationPack

#endif	// PRINT_VECTOR_CHANGE_STATS_H
