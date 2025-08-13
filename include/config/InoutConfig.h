#pragma once

#include "base/Global.h"

#include "libmesh/point.h"
#include <optional>

/**
 *   Holdes the probe information, to be used in the
 *   actual probe class (postproc/Probe)
 */
class ProbeConfig 
{
public:
  ProbeConfig( string & name_ ) : name(name_) {};

  string name;
  optional<string> type;

  // Planar
  optional<libMesh::Point> p0, p1, p2;

  // Linear
  optional<libMesh::Point> from, to;

  // Radial
  optional<libMesh::Point> center;
  optional<libMesh::Point> normal;
  optional<double> radius, dtheta;

  // Multiple
  optional<uint> npts;
  
  // Gauss
  vector<string>   boundaries;
  optional<string> order;
};

/**
 *
 *   Holds the model information for input and output.
 *
 */

struct InoutConfig 
{
    InoutConfig();
    ~InoutConfig();

    vector<ProbeConfig> probes;

    string model_dir;
};

/** ** ** ** **/
ostream& operator<<(ostream& os, const InoutConfig & m);
ostream& operator<<(ostream& os, const ProbeConfig & m);
