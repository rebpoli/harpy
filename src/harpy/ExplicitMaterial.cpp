
#include "harpy/ExplicitMaterial.h"
#include "libmesh/explicit_system.h"
#include "util/OutputOperators.h"
#include "postproc/TensorInvariants.h"

/**
 *   Reinits the FE strucutures and assign the element to the material.
 */
void ExplicitMaterial::reinit( const Elem & elem_, uint side ) 
{ 
  UNUSED(side);

  elem = &elem_;
  fe->reinit( elem );  
}

/**
 *    Compute the values in the quadrature points of the material.
 */
void ExplicitMaterial::eval( vector<double> & vals_qp , string vname )
{

  const std::vector<std::vector<Real>> & phi = fe->get_phi();

  uint vid = system.variable_number(vname);
  const DofMap & dof_map = system.get_dof_map();
  std::vector<dof_id_type> dofi;
  dof_map.dof_indices ( elem, dofi, vid);
  
  uint n_dofs = phi.size();

  uint nqp = qrule.n_points();
  vals_qp.resize( nqp, 0.0 );

  for(uint qp=0; qp<nqp; qp++) 
  for(uint B=0; B<n_dofs; B++) 
  {
    double v = (*system.current_local_solution)( dofi[B] );
    vals_qp[qp] += phi[B][qp] * v;
  }
}

/**
 *
 */
double ExplicitMaterial::eval( const Point & pt, string vname )
{
  // Use a local FEBase to avoid poluting the one owned by the object
  uint vid = system.variable_number(vname);
  const DofMap & dof_map = system.get_dof_map();
  std::vector<dof_id_type> dofi;
  dof_map.dof_indices ( elem, dofi, vid);
  
  FEType fetype = system.variable_type( vid );
  unique_ptr<FEBase> _fe = FEBase::build( 3, fetype );

  const std::vector<std::vector<Real>> & phi = _fe->get_phi();

  std::vector<Point> pts_from = { pt };
  std::vector<Point> pts_to;
  FEInterface::inverse_map (3, fetype, elem, pts_from, pts_to, 1E-10);
  _fe->reinit( elem, & pts_to );

  uint n_dofs = phi.size();

  double ret = 0;
  for ( uint B=0; B<n_dofs; B++ )
  {
    double v = (*system.current_local_solution)( dofi[B] );
    ret += phi[B][0] * v;
  }

  return ret;
}

/**
 *    Project the element coupler into the system. 
 */
void ExplicitMaterial::project( vector<double> & vals_qp, string vname )
{
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<Real> & jxw = fe->get_JxW();

  // Update DOF INDICES
  uint vid = system.variable_number(vname);
  const DofMap & dof_map = system.get_dof_map();
  std::vector<dof_id_type> dof_indices;
  dof_map.dof_indices ( elem, dof_indices, vid);

  //  M x = F
  uint n_dofs = phi.size();
  DenseMatrix<Real> MAT(n_dofs, n_dofs);
  DenseVector<Number> F(n_dofs);

  uint nqp = qrule.n_points();
  for(uint qp=0; qp<nqp; qp++) 
  for(uint B=0; B<n_dofs; B++) 
    F(B) += jxw[qp] * vals_qp[qp] * phi[B][qp];

  for(uint qp=0; qp<nqp; qp++) 
  for(uint B=0; B<n_dofs; B++) 
  for(uint M=0; M<n_dofs; M++) 
    MAT(B,M) += jxw[qp] * phi[B][qp] * phi[M][qp];

  for(uint B=0; B<n_dofs; B++) 
    if ( std::abs(MAT(B,B)) < 1e-10 ) MAT(B,B) = 1;

  DenseVector<Number> X;
//  MAT.cholesky_solve(F, X);
  MAT.lu_solve(F, X);

  // Cada processador seta os seus nos apenas
  dof_id_type f = system.solution->first_local_index();
  dof_id_type l = system.solution->last_local_index();
  for(uint B=0; B<n_dofs; B++)
  {
    dof_id_type dof_i = dof_indices[B];
    if ((f <= dof_i) && (dof_i < l )) 
      system.solution->set(dof_i, X(B));
  }
}

/**
 *    Project the element coupler into the system. 
 */
void ExplicitMaterial::project_tensor( vector<RealTensor> & vals_qp, string vname )
{
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<Real> & jxw = fe->get_JxW();

  vector<string> sufvec =         { "XX",  "YY",  "ZZ",  "XY",  "XZ",   "YZ" };
  vector<pair<uint,uint>> ijvec = { {0,0}, {1,1}, {2,2}, {0,1}, {0,2}, {1,2} };

  uint nqp = qrule.n_points();

  if ( vals_qp.size() != nqp ) 
    flog << "These quantities should be the same. Something went wrong (vname:" << vname << ", qp!=qr.n_points(): " << vals_qp.size() << "!=" << qrule.n_points() << ")";
//  dlog(1) << "EC vs QRULE: " << vals_qp.size() << " vs " << qrule.n_points(); 

  for ( uint a=0; a<sufvec.size(); a++ )
  {
    string suf = sufvec[a];
    auto [i,j] = ijvec[a];

    // Update DOF INDICES
    string trg_vname = vname + suf;
    uint vid = system.variable_number(trg_vname);
    const DofMap & dof_map = system.get_dof_map();
    std::vector<dof_id_type> dof_indices;
    dof_map.dof_indices ( elem, dof_indices, vid);

    //  M x = F
    uint n_dofs = phi.size();
    DenseMatrix<Real> MAT(n_dofs, n_dofs);
    DenseVector<Number> F(n_dofs);

    for(uint qp=0; qp<nqp; qp++) 
    for(uint B=0; B<n_dofs; B++) 
      F(B) += jxw[qp] * vals_qp[qp](i,j) * phi[B][qp];

    for(uint qp=0; qp<nqp; qp++) 
    for(uint B=0; B<n_dofs; B++) 
    for(uint M=0; M<n_dofs; M++) 
      MAT(B,M) += jxw[qp] * phi[B][qp] * phi[M][qp];

    for(uint B=0; B<n_dofs; B++) 
      if ( std::abs(MAT(B,B)) < 1e-10 ) MAT(B,B) = 1;

//    dlog(1) << endl << MAT;
    DenseVector<Number> X;
    MAT.cholesky_solve(F, X);
//    MAT.lu_solve(F, X);

    // Cada processador seta os seus nos apenas
    dof_id_type f = system.solution->first_local_index();
    dof_id_type l = system.solution->last_local_index();
    for(uint B=0; B<n_dofs; B++)
    {
      dof_id_type dof_i = dof_indices[B];
      if ((f <= dof_i) && (dof_i < l )) 
        system.solution->set(dof_i, X(B));
    }
  }
}


/**
 *
 *
 */
void ExplicitMaterial::project_tensor_invariants( vector<RealTensor> & vals_qp, string vname )
{
  uint nqp = qrule.n_points();

  /* Initialize the vectors to prepare the projection */
  vector<double> invarQ_qp(nqp), invarP_qp(nqp);
  vector<double> S1_qp(nqp), S2_qp(nqp), S3_qp(nqp);

  /* Principal stresses [i:012=S1,S2,S3 ; j:012=X,Y,Z] */
  vector<double> Sij_qp[3][3];
  for (uint i=0; i<3; i++) 
  for (uint j=0; j<3; j++) 
    Sij_qp[i][j] = vector<double>( nqp );

  /** Do the calculations **/
  for (uint qp=0; qp<nqp; qp++)
  {
    TensorInvariants tinv( vals_qp[qp] );

    // S1_j, S2_j, S3_j
    for (uint i=0; i<3; i++)
    for (uint j=0; j<3; j++) 
      Sij_qp[i][j][qp] = tinv.get_Sij( i, j );
    
    invarP_qp[qp] = tinv.get_P();
    invarQ_qp[qp] = tinv.get_Q();

    // Magnitude of the principal stresses
    S1_qp[qp] = tinv.S1_eval();
    S2_qp[qp] = tinv.S2_eval();
    S3_qp[qp] = tinv.S3_eval();
  }

  /* Project the invariants */
  vector< vector<string> > sij_ = { { vname+"1X", vname+"1Y", vname+"1Z" } , 
                                    { vname+"2X", vname+"2Y", vname+"2Z" } , 
                                    { vname+"3X", vname+"3Y", vname+"3Z" } };
  for (uint i=0; i<3; i++)
  for (uint j=0; j<3; j++) 
    project( Sij_qp[i][j], sij_[i][j] );

  project( invarQ_qp, vname+"_invarQ" );
  project( invarP_qp, vname+"_invarP" );

  project( S1_qp, vname+"1_mag" );
  project( S2_qp, vname+"2_mag" );
  project( S3_qp, vname+"3_mag" );
}

/**
 *
 */
void ExplicitMaterial::init_fem()
{
  SCOPELOG(1);
  system.print_info();
  uint vid = system.variable_number( "T" );
  DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(vid);

  fe = move( FEBase::build(3, fe_type) );

  // Inherits the qrule from the reference material to
  // calculate in the same qp's
  if ( refmat )
    qrule = refmat->qrule;
  else
    qrule = QGauss( 3, fe_type.default_quadrature_order() );

  fe->attach_quadrature_rule (&qrule);

  // Enable calculations calculations
  fe->get_JxW();
  fe->get_phi();
  fe->get_dphi();
  fe->get_xyz();
}

