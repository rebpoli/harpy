
#pragma once

#include "harpy/Global.h"
#include "postproc/probes/Probe.h"
#include "util/CsvFile.h"
#include "util/NetCDFWriter.h"

namespace solver { namespace viscoplastic { class ViscoplasticSolver; } }


namespace postproc {
namespace report {

/**
 *
 *   This class bridges the viscoplastic solver and the probes.
 *   The probes hold the points that must be evaluated.
 *   The solver holds the mesh and the solution.
 *
 */

using postproc::probes::Probe;
using postproc::probes::GaussProbe;
using postproc::probes::ProbeCol;
using solver::viscoplastic::ViscoplasticSolver;

class ViscoplasticReport
{
public:
  ViscoplasticReport( ViscoplasticSolver & solver_ );

  void init();
  void init_material( Probe & probe  );
  void do_export();
  void export_by_point( Probe & probe );
  void export_by_face( GaussProbe & probe );
  void export_scalars();

private:

  ViscoplasticSolver & solver;
  ProbeCol probes;

  string scalars_fn;
};

}} // ns
