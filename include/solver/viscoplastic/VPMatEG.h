#pragma once

#include "harpy/Global.h"

#include "libmesh/transient_system.h"
#include "solver/viscoplastic/EGDataStruct.h"

namespace solver {
namespace viscoplastic {

class ViscoplasticSolver;
using namespace libMesh;

/**
 *
 */
class VPMatEG 
{
public:
  VPMatEG( TransientNonlinearImplicitSystem & sys_,
           ViscoplasticSolver & vpsolver_ );

  void init();                             /// Create datastructures from the mesh

private:
  vector< EGFacePair > gamma_I;        /// Internal skeleton
  vector< EGFace > gamma_H;            /// External boundaries

  TransientNonlinearImplicitSystem & system;
  ViscoplasticSolver & vpsolver; 

  EGFEM fem_p, fem_n;                      /// Shape functions to be reinit'ed

private:
};


}} // ns
