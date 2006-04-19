#include "GLpApp_AdvDiffReactOptModel.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"

// For orthogonalization of the basis B_bar
#include "sillyModifiedGramSchmidt.hpp" // This is just an example!
#include "Thyra_EpetraThyraWrappers.hpp"

//#define GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF

namespace {

inline double sqr( const double& s ) { return s*s; }

inline double dot( const Epetra_Vector &x, const Epetra_Vector &y )
{
  double dot[1];
  x.Dot(y,dot);
  return dot[0];
}

} // namespace

namespace GLpApp {

AdvDiffReactOptModel::AdvDiffReactOptModel(
  Teuchos::RefCountPtr<GLpApp::GLpYUEpetraDataPool>   const& dat
  ,const int                                                 np
  ,const double                                              x0
  ,const double                                              p0
  ,const double                                              reactionRate
  )
  :dat_(dat),np_(np),reactionRate_(reactionRate)
{

#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  *out << "\nEntering AdvDiffReactOptModel::AdvDiffReactOptModel(...) ...\n";
#endif

  //
  // Get the maps
  //
  const Epetra_SerialComm serialComm;
  const Epetra_Comm &comm = dat_->getA()->OperatorDomainMap().Comm();
  map_x_ = Teuchos::rcp(new Epetra_Map(dat_->getA()->OperatorDomainMap()));
  map_p_bar_ = Teuchos::rcp(new Epetra_Map(dat_->getB()->OperatorDomainMap()));
  map_f_ = Teuchos::rcp(new Epetra_Map(dat_->getA()->OperatorRangeMap()));
  map_g_ = Teuchos::rcp(new Epetra_Map(1,1,0,comm));
  //
  // Initialize the basis matrix for p such that p_bar = B_bar * p
  //
  TEST_FOR_EXCEPTION(
    np_ > map_p_bar_->NumGlobalElements(), std::logic_error
    ,"Error, np="<<np_<<" can not be greater than map_p_bar_->NumGlobalElements()="
    <<map_p_bar_->NumGlobalElements()<<"!"
    );
  if( np_ > 0 ) {
    map_p_ = Teuchos::rcp(new Epetra_Map(np_,np_,0,comm));
    B_bar_ = Teuchos::rcp(new Epetra_MultiVector(*map_p_bar_,np_));
    if( np_ == 1 ) {
      // Just make B_bar a column with ones!
      B_bar_->PutScalar(1.0);
    }
    else if( np_ > 1 ) {
      //
      // Create a random local B_bar that will be the same no matter how the
      // problem is distributed.
      //
      typedef Teuchos::ScalarTraits<double> ST;
      const int numBndyNodes        = map_p_bar_->NumMyElements();
      const int *bndyNodeGlobalIDs  = map_p_bar_->MyGlobalElements();
      for( int i = 0; i < numBndyNodes; ++i ) {
        (*B_bar_)[0][i] = 1.0;
        ST::seedrandom(bndyNodeGlobalIDs[i]);
        for( int j = 1; j < np_; ++j ) {
          (*B_bar_)[j][i] = ST::random();
        }
      }
      //
      // Use modified Gram-Schmidt to create an orthonormal version of B_bar!
      //
      Teuchos::RefCountPtr<Thyra::MultiVectorBase<double> >
        thyra_B_bar = Thyra::create_MPIMultiVectorBase(
          B_bar_
          ,Thyra::create_MPIVectorSpaceBase(Teuchos::rcp(new Epetra_Map(*map_p_bar_)))
          ,Thyra::create_MPIVectorSpaceBase(Teuchos::rcp(new Epetra_Map(*map_p_)))
          ),
        thyra_fact_R;
      sillyModifiedGramSchmidt(&*thyra_B_bar,&thyra_fact_R);
      // We just discard the "R" factory thyra_fact_R
    }
    else {
      TEST_FOR_EXCEPT(true); // Should never get here!
    }
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
    *out << "\nB_bar =\n\n";
    B_bar_->Print(*Teuchos::OSTab(out).getOStream());
#endif
  }
  else {
    // B_bar = I implicitly!
    map_p_ = map_p_bar_;
  }
  //
  // Create vectors
  //
  x0_ = Teuchos::rcp(new Epetra_Vector(*map_x_));
  //xL_ = Teuchos::rcp(new Epetra_Vector(*map_x_));
  //xU_ = Teuchos::rcp(new Epetra_Vector(*map_x_));
  p0_ = Teuchos::rcp(new Epetra_Vector(*map_p_));
  //pL_ = Teuchos::rcp(new Epetra_Vector(*map_p_));
  //pU_ = Teuchos::rcp(new Epetra_Vector(*map_p_));
  //
  // Initialize the vectors
  //
  x0_->PutScalar(x0);
  p0_->PutScalar(p0);
  //
  // Initialize the graph for W
  //
  dat_->computeNpy(x0_);
  //if(dumpAll) { *out << "\nA =\n"; { Teuchos::OSTab tab(out); dat_->getA()->Print(*out); } }
  //if(dumpAll) { *out << "\nNpy =\n"; {  Teuchos::OSTab tab(out); dat_->getNpy()->Print(*out); } }
  W_graph_ = Teuchos::rcp(new Epetra_CrsGraph(dat_->getA()->Graph())); // Assume A and Npy have same graph!
  //
  // Get default objective matching vector q
  //
  q_ = Teuchos::rcp(new Epetra_Vector(*(*dat_->getq())(0))); // From Epetra_FEVector to Epetra_Vector!
  //
  isInitialized_ = true;
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
  *out << "\nLeaving AdvDiffReactOptModel::AdvDiffReactOptModel(...) ...\n";
#endif
}

void AdvDiffReactOptModel::set_q( Teuchos::RefCountPtr<Epetra_Vector> const& q )
{
  q_ = q;
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RefCountPtr<const Epetra_Map>
AdvDiffReactOptModel::get_x_map() const
{
  return map_x_;
}

Teuchos::RefCountPtr<const Epetra_Map>
AdvDiffReactOptModel::get_f_map() const
{
  return map_x_;
}

Teuchos::RefCountPtr<const Epetra_Map>
AdvDiffReactOptModel::get_p_map(int l) const
{
  TEST_FOR_EXCEPT(l!=0);
  return map_p_;
}

Teuchos::RefCountPtr<const Epetra_Map>
AdvDiffReactOptModel::get_g_map(int j) const
{
  TEST_FOR_EXCEPT(j!=0);
  return map_g_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_x_init() const
{
  return x0_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_p_init(int l) const
{
  TEST_FOR_EXCEPT(l!=0);
  return p0_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_x_lower_bounds() const
{
  return xL_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_x_upper_bounds() const
{
  return xU_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_p_lower_bounds(int l) const
{
  TEST_FOR_EXCEPT(l!=0);
  return pL_;
}

Teuchos::RefCountPtr<const Epetra_Vector>
AdvDiffReactOptModel::get_p_upper_bounds(int l) const
{
  TEST_FOR_EXCEPT(l!=0);
  return pU_;
}

Teuchos::RefCountPtr<Epetra_Operator>
AdvDiffReactOptModel::create_W() const
{
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));
}

Teuchos::RefCountPtr<Epetra_Operator>
AdvDiffReactOptModel::create_DfDp_op(int l) const
{
  TEST_FOR_EXCEPT(l!=0);
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy,dat_->getB()->Graph()));
  // See DfDp in evalModel(...) below for details
}

EpetraExt::ModelEvaluator::InArgs
AdvDiffReactOptModel::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(1);
  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
AdvDiffReactOptModel::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1,1);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      reactionRate_!=0.0 ? DERIV_LINEARITY_NONCONST : DERIV_LINEARITY_CONST
      ,DERIV_RANK_FULL
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(
    OUT_ARG_DfDp,0
    ,( np_ > 0
       ? DerivativeSupport(DERIV_MV_BY_COL)
       : DerivativeSupport(DERIV_LINEAR_OP,DERIV_MV_BY_COL)
      )
    );
  outArgs.set_DfDp_properties(
    0,DerivativeProperties(
      DERIV_LINEARITY_CONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DgDx,0,DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDx_properties(
    0,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DgDp,0,0,DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDp_properties(
    0,0,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  return outArgs;
}

void AdvDiffReactOptModel::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
  Teuchos::RefCountPtr<Teuchos::FancyOStream>
    dout = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab dtab(dout);
  *dout << "\nEntering AdvDiffReactOptModel::evalModel(...) ...\n";
#endif
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;
  //
  Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
  const bool trace = ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_LOW) );
  const bool dumpAll = ( static_cast<int>(this->getVerbLevel()) >= static_cast<int>(Teuchos::VERB_EXTREME) );
  //
  Teuchos::OSTab tab(out);
  if(out.get() && trace) *out << "\n*** Entering AdvDiffReactOptModel::evalModel(...) ...\n"; 
  //
  // Get the input arguments
  //
  const Epetra_Vector *p_in = inArgs.get_p(0).get();
  const Epetra_Vector &p = (p_in ? *p_in : *p0_);
  const Epetra_Vector &x = *inArgs.get_x();
  //
  // Get the output arguments
  //
  Epetra_Vector       *f_out = outArgs.get_f().get();
  Epetra_Vector       *g_out = outArgs.get_g(0).get();
  Epetra_Operator     *W_out = outArgs.get_W().get();
  Derivative          DfDp_out = outArgs.get_DfDp(0);
  Epetra_MultiVector  *DgDx_trans_out = get_DgDx_mv(0,outArgs,DERIV_TRANS_MV_BY_ROW).get();
  Epetra_MultiVector  *DgDp_trans_out = get_DgDp_mv(0,0,outArgs,DERIV_TRANS_MV_BY_ROW).get();
  //
  // Precompute some shared quantities
  //
  // p_bar = B_bar*p
  Teuchos::RefCountPtr<const Epetra_Vector> p_bar;
  if(np_ > 0) {
    Teuchos::RefCountPtr<Epetra_Vector> _p_bar;
    _p_bar = Teuchos::rcp(new Epetra_Vector(*map_p_bar_));
    _p_bar->Multiply('N','N',1.0,*B_bar_,p,0.0);
    p_bar = _p_bar;
  }
  else {
    p_bar = Teuchos::rcp(&p,false);
  }
  // R_p_bar = R*p_bar = R*(B_bar*p)
  Teuchos::RefCountPtr<const Epetra_Vector> R_p_bar;
  if( g_out || DgDp_trans_out ) {
      Teuchos::RefCountPtr<Epetra_Vector>
      _R_p_bar = Teuchos::rcp(new Epetra_Vector(*map_p_bar_));
    dat_->getR()->Multiply(false,*p_bar,*_R_p_bar);
    R_p_bar = _R_p_bar;
  }
  //
  // Compute the functions
  //
  if(f_out) {
    //
    // f = A*x + reationRate*Ny(x) + B*(B_bar*p)
    //
    Epetra_Vector &f = *f_out;
    Epetra_Vector Ax(*map_f_);
    dat_->getA()->Multiply(false,x,Ax);
    f.Update(1.0,Ax,0.0);
    if(reactionRate_!=0.0) {
      dat_->computeNy(Teuchos::rcp(&x,false));
      f.Update(reactionRate_,*dat_->getNy(),1.0);
    }
    Epetra_Vector Bp(*map_f_);
    dat_->getB()->Multiply(false,*p_bar,Bp);
    f.Update(1.0,Bp,-1.0,*dat_->getb(),1.0);
  }
  if(g_out) {
    //
    // g = 0.5 * (x-q)'*H*(x-q) + 0.5*regBeta*(B_bar*p)'*R*(B_bar*p)
    //
    Epetra_Vector &g = *g_out;
    Epetra_Vector xq(x);
    xq.Update(-1.0, *q_, 1.0);
    Epetra_Vector Hxq(x);
    dat_->getH()->Multiply(false,xq,Hxq);
    g[0] = 0.5*dot(xq,Hxq) + 0.5*dat_->getbeta()*dot(*p_bar,*R_p_bar);
  }
  if(W_out) {
    //
    // W = A + reationRate*Npy(x)
    //
    Epetra_CrsMatrix &DfDx = dyn_cast<Epetra_CrsMatrix>(*W_out);
    if(reactionRate_!=0.0)
      dat_->computeNpy(Teuchos::rcp(&x,false));
    Teuchos::RefCountPtr<Epetra_CrsMatrix>
      dat_A = dat_->getA(),
      dat_Npy = dat_->getNpy();
    const int numMyRows = dat_A->NumMyRows();
    for( int i = 0; i < numMyRows; ++i ) {
      int dat_A_num_row_entries=0; double *dat_A_row_vals=0; int *dat_A_row_inds=0;
      dat_A->ExtractMyRowView(i,dat_A_num_row_entries,dat_A_row_vals,dat_A_row_inds);
      int DfDx_num_row_entries=0; double *DfDx_row_vals=0; int *DfDx_row_inds=0;
      DfDx.ExtractMyRowView(i,DfDx_num_row_entries,DfDx_row_vals,DfDx_row_inds);
#ifdef _DEBUG
      TEST_FOR_EXCEPT(DfDx_num_row_entries!=dat_A_num_row_entries);
#endif
      if(reactionRate_!=0.0) {
        int dat_Npy_num_row_entries=0; double *dat_Npy_row_vals=0; int *dat_Npy_row_inds=0;
        dat_Npy->ExtractMyRowView(i,dat_Npy_num_row_entries,dat_Npy_row_vals,dat_Npy_row_inds);
#ifdef _DEBUG
        TEST_FOR_EXCEPT(dat_A_num_row_entries!=dat_Npy_num_row_entries);
#endif
        for(int k = 0; k < DfDx_num_row_entries; ++k) {
#ifdef _DEBUG
          TEST_FOR_EXCEPT(dat_A_row_inds[k]!=dat_Npy_row_inds[k]||dat_A_row_inds[k]!=DfDx_row_inds[k]);
#endif
          DfDx_row_vals[k] = dat_A_row_vals[k] + reactionRate_ * dat_Npy_row_vals[k];
        }
      }
      else {
        for(int k = 0; k < DfDx_num_row_entries; ++k) {
#ifdef _DEBUG
          TEST_FOR_EXCEPT(dat_A_row_inds[k]!=DfDx_row_inds[k]);
#endif
          DfDx_row_vals[k] = dat_A_row_vals[k];
        }
      }
    }
  }
  if(!DfDp_out.isEmpty()) {
    if(out.get() && trace) *out << "\nComputing DfDp ...\n"; 
    //
    // DfDp = B*B_bar
    //
    Epetra_CrsMatrix   *DfDp_op = NULL;
    Epetra_MultiVector *DfDp_mv = NULL;
    if(out.get() && dumpAll)
    { *out << "\nB =\n"; { Teuchos::OSTab tab(out); dat_->getB()->Print(*out); } }
    if(DfDp_out.getLinearOp().get()) {
      DfDp_op = &dyn_cast<Epetra_CrsMatrix>(*DfDp_out.getLinearOp());
    }
    else {
      DfDp_mv = &*DfDp_out.getDerivativeMultiVector().getMultiVector();
      DfDp_mv->PutScalar(0.0);
    }
    Teuchos::RefCountPtr<Epetra_CrsMatrix>
      dat_B = dat_->getB();
    if(np_ > 0) {
      //
      // We only support a Multi-vector form when we have a non-idenity basis
      // matrix B_bar for p!
      //
      TEST_FOR_EXCEPT(DfDp_mv==NULL);
      dat_B->Multiply(false,*B_bar_,*DfDp_mv);
    }
    else {
      //
      // Note: We copy from B every time in order to be safe.  Note that since
      // the client knows that B is constant (sense we told them so in
      // createOutArgs()) then it should only compute this matrix once and keep
      // it if it is smart.
      //
      // Note: We support both the CrsMatrix and MultiVector form of this object
      // to make things easier for the client.
      //
      if(DfDp_op) {
        const int numMyRows = dat_B->NumMyRows();
        for( int i = 0; i < numMyRows; ++i ) {
          int dat_B_num_row_entries=0; double *dat_B_row_vals=0; int *dat_B_row_inds=0;
          dat_B->ExtractMyRowView(i,dat_B_num_row_entries,dat_B_row_vals,dat_B_row_inds);
          int DfDp_num_row_entries=0; double *DfDp_row_vals=0; int *DfDp_row_inds=0;
          DfDp_op->ExtractMyRowView(i,DfDp_num_row_entries,DfDp_row_vals,DfDp_row_inds);
#ifdef _DEBUG
          TEST_FOR_EXCEPT(DfDp_num_row_entries!=dat_B_num_row_entries);
#endif
          for(int k = 0; k < DfDp_num_row_entries; ++k) {
#ifdef _DEBUG
            TEST_FOR_EXCEPT(dat_B_row_inds[k]!=DfDp_row_inds[k]);
#endif
            DfDp_row_vals[k] = dat_B_row_vals[k];
          }
          // ToDo: The above code should be put in a utility function called copyValues(...)!
        }
      }
      else if(DfDp_mv) {
        // We must do a mat-vec to get this since we can't just copy out the
        // matrix entries since the domain map may be different from the
        // column map!  I learned this the very very hard way!  I am using
        // Thyra wrappers here since I can't figure out for the life of me how
        // to do this cleanly with Epetra alone!
        Teuchos::RefCountPtr<Epetra_Vector>
          etaVec = Teuchos::rcp(new Epetra_Vector(*map_p_bar_));
        Teuchos::RefCountPtr<const Thyra::MPIVectorSpaceBase<double> >
          space_p_bar = Thyra::create_MPIVectorSpaceBase(Teuchos::rcp(new Epetra_Map(*map_p_bar_)));
        Teuchos::RefCountPtr<Thyra::VectorBase<double> >
          thyra_etaVec = Thyra::create_MPIVectorBase(etaVec,space_p_bar);
        for( int i = 0; i < map_p_bar_->NumGlobalElements(); ++i ) {
          Thyra::assign(&*thyra_etaVec,0.0);
          Thyra::set_ele(i,1.0,&*thyra_etaVec);
          dat_B->Multiply(false,*etaVec,*(*DfDp_mv)(i));
        };
      }
    }
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
    if(DfDp_op) {
      *dout << "\nDfDp_op =\n\n";
      DfDp_op->Print(*Teuchos::OSTab(dout).getOStream());
    }
    if(DfDp_mv) {
      *dout << "\nDfDp_mv =\n\n";
      DfDp_mv->Print(*Teuchos::OSTab(dout).getOStream());
    }
#endif


  }
  if(DgDx_trans_out) {
    //
    // DgDx' = H*(x-q)
    //
    Epetra_Vector &DgDx_trans = *(*DgDx_trans_out)(0);
    Epetra_Vector xq(x);
    xq.Update(-1.0,*q_,1.0);
    dat_->getH()->Multiply(false,xq,DgDx_trans);
  }
  if(DgDp_trans_out) {
    //
    // DgDp' = regBeta*B_bar'*(R*(B_bar*p))
    //
    Epetra_Vector &DgDp_trans = *(*DgDp_trans_out)(0);
    if(np_ > 0) {
      DgDp_trans.Multiply('T','N',dat_->getbeta(),*B_bar_,*R_p_bar,0.0);
    }
    else {
      DgDp_trans.Update(dat_->getbeta(),*R_p_bar,0.0);
    }
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
    *dout << "\nR_p_bar =\n\n";
    R_p_bar->Print(*Teuchos::OSTab(dout).getOStream());
    if(B_bar_.get()) {
      *dout << "\nB_bar =\n\n";
      B_bar_->Print(*Teuchos::OSTab(dout).getOStream());
    }
    *dout << "\nDgDp_trans =\n\n";
    DgDp_trans.Print(*Teuchos::OSTab(dout).getOStream());
#endif
  }
  if(out.get() && trace) *out << "\n*** Leaving AdvDiffReactOptModel::evalModel(...) ...\n"; 
#ifdef GLPAPP_ADVDIFFREACT_OPTMODEL_DUMP_STUFF
  *dout << "\nLeaving AdvDiffReactOptModel::evalModel(...) ...\n";
#endif
}

} // namespace GLpApp
