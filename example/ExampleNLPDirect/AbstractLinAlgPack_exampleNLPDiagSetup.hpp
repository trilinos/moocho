// ////////////////////////////////////////////////////////////////////
// exampleNLPDiagSetup.hpp

#ifndef ALAP_EXPL_NLP_DIAG_SETUP_HPP
#define ALAP_EXPL_NLP_DIAG_SETUP_HPP

#include <assert.h>

#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "AbstractLinAlgPack/src/AbstractLinAlgPackTypes.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace AbstractLinAlgPack {

///
/** Create a vector space given the input arguments argc, argv[] and an MPI communicator
 *
 */
int exampleNLPDiagSetup(
	int argc, char* argv[], MPI_Comm comm
	,Teuchos::RefCountPtr<const VectorSpace> *vec_space
	,size_type *n, value_type *xo, bool *has_bounds, bool *dep_bounded
	);

} // namespace AbstractLinAlgPack

#endif // ALAP_EXPL_NLP_DIAG_SETUP_HPP
