
#pragma once

#include "base/Global.h"
#include "postproc/Probe.h"

/**
 *
 *   This class bridges the viscoplastic solver and the probes.
 *   The probes hold the points that must be evaluated.
 *   The solver holds the mesh and the solution.
 *
 */

class ViscoplasticSolver;

class ViscoplasticReport
{
public:
  ViscoplasticReport( ViscoplasticSolver & solver_ ) : solver(solver_) {}

  void do_export();
  void export_by_point( Probe & probe );
  void export_by_face( GaussProbe & probe );

private:

  ViscoplasticSolver & solver;
  ProbeCol probes;
};
