
#include "harpy/Material.h"
#include "harpy/Solver.h"
#include "config/ModelConfig.h" // MODEL global var
#include "config/SolverConfig.h"
#include "material/ViscoPlasticMaterial.h"

#include "libmesh/elem.h"
#include "libmesh/mesh.h"
#include "libmesh/explicit_system.h"

#include <regex>


// Some local Regex
const string prop_name = R"(([-a-zA-Z_0-9]+))";
const regex RE_PROPERTY ( R"(^)" + prop_name + R"(\.)" + prop_name + R"($)");

/**
 *
 *
 */
Material::Material( suint sid_, const MaterialConfig & config_ ) :
                          config(config_), name(config.name),
                          sid(sid_), qrule(3), elem_coupler(0)
{ 
}

/**
 *
 */
Material::~Material() {};

/**
 *
 *    This function needs the FEBase member to be initialized in the Material to fetch the xyz.
 */
void Material::init_coupler( Elem * elem, ElemCoupler & ec )
{
  // Calc xyz
  const std::vector<Point> & xyz = fe->get_xyz();
  fe->reinit( elem );

  // Feed the coupler
  for ( auto & ctx_pname : required_material_properties )
  {
    smatch match;
    if ( ! regex_search( ctx_pname, match, RE_PROPERTY ) ) flog << "Invalid property name: " << ctx_pname;
    string context = match[1], pname = match[2];

    const MaterialConfig & mconf = config;
    config.get_property( ec.dbl_params[ctx_pname], pname, xyz, context );
  }
}

/**
 *
 */
void Material::get_from_element_coupler( string vname, vector<double> & curr,  vector<double> & old )
{
  if ( ! elem_coupler ) flog << "Element coupler not initialized! Something is wrong.";

  if ( ! elem_coupler->dbl_params.count( vname ) ) flog << "Cannot find variable '" << vname << "' in element coupler. Cannot continue." ;
  if ( ! elem_coupler->dbl_params_old.count( vname ) ) flog << "Cannot find variable '" << vname << "' in element coupler (OLD). Cannot continue." ;

  curr = elem_coupler->dbl_params.at( vname );
  old = elem_coupler->dbl_params_old.at( vname );
}

/**
 *
 */
void Material::get_from_element_coupler( string vname, vector<double> & curr )
{
  if ( ! elem_coupler ) flog << "Element coupler not initialized! Something is wrong.";
  if ( ! elem_coupler->dbl_params.count( vname ) ) {
    dlog(1) << *elem_coupler ;
    flog << "Cannot find variable '" << vname << "' in element coupler. Cannot continue." ;
  }

  curr = elem_coupler->dbl_params.at( vname );
}

/**
 *
 */
void Material::get_from_element_coupler( string vname, vector<RealVectorValue> & curr )
{
  if ( ! elem_coupler ) flog << "Element coupler not initialized! Something is wrong.";
  if ( ! elem_coupler->vector_params.count( vname ) ) flog << "Cannot find variable '" << vname << "' in element coupler. Cannot continue." ;
  curr = elem_coupler->vector_params.at( vname );
}

/**
 *
 */
void Material::get_from_element_coupler( string vname, vector<RealTensor> & curr )
{
  if ( ! elem_coupler ) flog << "Element coupler not initialized! Something is wrong.";
  if ( ! elem_coupler->tensor_params.count( vname ) ) flog << "Cannot find variable '" << vname << "' in element coupler. Cannot continue." ;
  curr = elem_coupler->tensor_params.at( vname );
}

/**
 *    Project the element coupler into the system. 
 */
void MaterialExplicit::project( ElemCoupler & ec, string vname_coupler, string vname_system )
{
  if ( ! vname_system.length() ) vname_system = vname_coupler;

  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<Real> & jxw = fe->get_JxW();
  const vector<double> & vals_qp = ec.dbl_params[vname_coupler];
  uint vid = system.variable_number(vname_system);
  Elem & elem = system.get_mesh().elem_ref( ec.eid );

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

  const DofMap & dof_map = system.get_dof_map();
  std::vector<dof_id_type> dof_indices;
  dof_map.dof_indices (&elem, dof_indices, vid);

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
void MaterialExplicit::project_tensor( ElemCoupler & ec, string vname_coupler, string vname_system )
{
  if ( ! vname_system.length() ) vname_system = vname_coupler;

  vector<RealTensor> & vals_qp = ec.tensor_params[vname_coupler];
  uint nqp = vals_qp.size();

  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<Real> & jxw = fe->get_JxW();
  Elem & elem = system.get_mesh().elem_ref( ec.eid );

  vector<string> sufvec =         { "XX",  "YY",  "ZZ",  "XY",  "XZ",   "YZ" };
  vector<pair<uint,uint>> ijvec = { {0,0}, {1,1}, {2,2}, {0,1}, {0,2}, {1,2} };

  if ( vals_qp.size() != qrule.n_points() ) 
    flog << "These quantities should be the same. Something went wrong (vname:" << vname_coupler << ", qp!=qr.n_points(): " << vals_qp.size() << "!=" << qrule.n_points() << ")";
//  dlog(1) << "EC vs QRULE: " << vals_qp.size() << " vs " << qrule.n_points(); 

  for ( uint a=0; a<sufvec.size(); a++ )
  {
    string suf = sufvec[a];
    auto [i,j] = ijvec[a];

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

    string trg_vname = vname_system + suf;
    uint vid = system.variable_number(trg_vname);

    const DofMap & dof_map = system.get_dof_map();
    std::vector<dof_id_type> dof_indices;
    dof_map.dof_indices (&elem, dof_indices, vid);

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

