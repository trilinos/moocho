#ifndef GLPAPP_GLPYUEPETRADATAPOOL_H
#define GLPAPP_GLPYUEPETRADATAPOOL_H

//#include "Epetra_config.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_LinearProblem.h"
//#include "AztecOO.h"
//#include "AztecOO_Operator.h"
#include "Epetra_LAPACK.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include "usr_par.h"


#include "GenSQP_DataPool.hpp"
#include "GenSQP_YUEpetraVector.hpp"
#include "includes.h"

namespace GLpApp {

/**
    \brief Implements the GenSQP::DataPool interface module for the parallel
    Ginzburg-Landau (GLp) application.
*/
class GLpYUEpetraDataPool : public GenSQP::DataPool {
public:

  GLpYUEpetraDataPool(
    Teuchos::RefCountPtr<const Epetra_Comm>    const& commptr
    ,const double                              beta
    ,const double                              len_x     // Ignored if myfile is *not* empty
    ,const double                              len_y     // Ignored if myfile is *not* empty
    ,const int                                 local_nx  // Ignored if myfile is *not* empty
    ,const int                                 local_ny  // Ignored if myfile is *not* empty
    ,const char                                myfile[]
    ,const bool                                trace
    );

  /** \brief Calls functions to compute nonlinear quantities and the augmented system matrix.
             These computations are performed after every update of the SQP iterate.
  */
  void computeAll( const GenSQP::Vector &x );

  /** \brief Solves augmented system.*/
  int  solveAugsys( const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsy,
                    const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsu,
                    const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsp,
                    const Teuchos::RefCountPtr<Epetra_MultiVector> & y,
                    const Teuchos::RefCountPtr<Epetra_MultiVector> & u,
                    const Teuchos::RefCountPtr<Epetra_MultiVector> & p,
                    double tol );

  Teuchos::RefCountPtr<const Epetra_Comm> getCommPtr();

  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getA();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getB();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getH();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getR();
  Teuchos::RefCountPtr<Epetra_CrsMatrix> getAugmat();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getNpy();

  Teuchos::RefCountPtr<Epetra_FEVector> getb();
  Teuchos::RefCountPtr<Epetra_FEVector> getq();
  Teuchos::RefCountPtr<Epetra_FEVector> getNy();

  /** \brief Calls the function that computes the nonlinear term.*/
  void computeNy(const Teuchos::RefCountPtr<const Epetra_MultiVector> & y);

  /** \brief Calls the function that computes the Jacobian of the nonlinear term.*/
  void computeNpy(const Teuchos::RefCountPtr<const Epetra_MultiVector> & y);

  /** \brief Assembles the augmented system (KKT-type) matrix.*/
  void computeAugmat();
  
  Teuchos::RefCountPtr<const Epetra_SerialDenseMatrix> getipcoords();
  Teuchos::RefCountPtr<const Epetra_IntSerialDenseVector> getipindx();
  Teuchos::RefCountPtr<const Epetra_SerialDenseMatrix> getpcoords();
  Teuchos::RefCountPtr<const Epetra_IntSerialDenseVector> getpindx();
  Teuchos::RefCountPtr<const Epetra_IntSerialDenseMatrix> gett();
  Teuchos::RefCountPtr<const Epetra_IntSerialDenseMatrix> gete();

  double getbeta();
  
  /** \brief Outputs the solution vector to files.*/
  void PrintVec( const Teuchos::RefCountPtr<const Epetra_Vector> & x );

private:

  Teuchos::RefCountPtr<const Epetra_Comm> commptr_;
          
  /** \brief Coordinates of nodes that are unique to this subdomain.*/
  Teuchos::RefCountPtr<Epetra_SerialDenseMatrix> ipcoords_;
  /** \brief Global nodes (interior, nonoverlapping) in this subdomain.*/
  Teuchos::RefCountPtr<Epetra_IntSerialDenseVector> ipindx_;
  /** \brief Coordinates of all nodes in this subdomain.*/
  Teuchos::RefCountPtr<Epetra_SerialDenseMatrix> pcoords_;
  /** \brief Global nodes (interior + shared, overlapping) in this subdomain.*/
  Teuchos::RefCountPtr<Epetra_IntSerialDenseVector> pindx_;
  /** \brief Elements (this includes all overlapping nodes).*/
  Teuchos::RefCountPtr<Epetra_IntSerialDenseMatrix> t_;
  /** \brief Edges.*/
  Teuchos::RefCountPtr<Epetra_IntSerialDenseMatrix> e_;

  /** \brief Volume stiffness matrix.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> A_;
  /** \brief Control/state mass matrix.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> B_;
  /** \brief Volume mass matrix.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> H_;
  /** \brief Edge mass matrix.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> R_;

  /** \brief Basis matrix for p_bar=B*p.*/
  Teuchos::RefCountPtr<Epetra_MultiVector> B_bar_;

  /** \brief Augmented system matrix:
   [ I  Jac* ]
   [Jac  0   ]
  */
  Teuchos::RefCountPtr<Epetra_CrsMatrix> Augmat_;

  /** \brief Jacobian of the nonlinear term.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> Npy_;

  /** \brief Right-hand side of the PDE.*/
  Teuchos::RefCountPtr<Epetra_FEVector> b_;
  /** \brief The desired state.*/
  Teuchos::RefCountPtr<Epetra_FEVector> q_;

  Teuchos::RefCountPtr<Epetra_FEVector> Ny_;

  /** \brief Regularization parameter.*/
  double beta_;

};

} // namespace GLpApp

#endif
