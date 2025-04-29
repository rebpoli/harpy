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

  // Radial
  optional<libMesh::Point> center, normal;
  optional<double> radius, dtheta;

  // Linear
  optional<libMesh::Point> from, to;
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
