// ///////////////////////////////////////////////////////
// vector_change_stats.h

#ifndef VECTOR_CHANGE_STATS_H
#define VECTOR_CHANGE_STATS_H

#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

///
/** Compute statistics for change in a vector.
  *
  * Given two vectors x and d where we wish to generate statistics
  * for the update x+d this function computes the following
  * quantitines:
  *
  * max( |d(i)|/(1+|x(i)|), i=1...n )
  *   => #max_k# = k, #max_term# = |d(k)|/(1+|x(k)|) <= 1\\
  * #min( |d(i)|/(1+|x(i)|), i=1...n ) => #min_k# = k, #min_term# = |d(k)|/(1+|x(k)|)#\\
  * #average( |d(i)|/(1+|x(i)|), i=1...10 )# => #av_term#\\
  * 
  * The purpose of generating these satistics is to determine
  * by how much x+d differs from x.
  *
  * If |d(i)|/|x(i)| < mach_eps with x(i) > 0 then we know that d(i) will
  * be lost when added to x(i) so x(i) + d(i) == x(i).
  *
  */
void vector_change_stats( const VectorSlice& x, const VectorSlice& d
	, value_type* max_term, size_type* max_k
	, value_type* min_term, size_type* min_k
	, value_type* av_term );

}	// end namespace ConstrainedOptimizationPack

#endif	// VECTOR_CHANGE_STATS_H
