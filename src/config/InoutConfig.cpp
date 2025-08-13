
#include "config/InoutConfig.h"
#include "config/reader/InoutReader.h"

#include "config/ModelConfig.h"

#include "util/OutputOperators.h"
#include "util/String.h"

using namespace harpy_string;


/**
 *   Calls the reader to fill in this object.
 */
InoutConfig::InoutConfig() : model_dir(MODEL->model_dir) 
{
  InoutReader reader(*this);
}


/**
 *
 */
InoutConfig::~InoutConfig() { }

/**
 *
 */
ostream& operator<<(ostream& os, const InoutConfig & m)
{
  os << "InoutConfig" << endl;
  for ( auto &  v : m.probes ) 
    os << right << setw(10) << "Probe '" + v.name + "'" << " : " << endl << v;
  return os;
}

/**
 *
 */
ostream& operator<<(ostream& os, const ProbeConfig & m)
{
  os << right << setw(20) << "Type"       << setw(20) << " " << m.type << endl;

  if (! m.type ) return os;
  string type = *m.type;

  if ( iequals(type, "planar") ) 
  {
    os << right << setw(20) << "P0"     << setw(20) << " " << m.p0 << endl;
    os << right << setw(20) << "P1"     << setw(20) << " " << m.p1 << endl;
    os << right << setw(20) << "P2"     << setw(20) << " " << m.p2 << endl;
    os << right << setw(20) << "Npts"   << setw(20) << " " << m.npts << endl;
  }

  if ( iequals(type, "radial") ) 
  {
    os << right << setw(20) << "Center"     << setw(20) << " " << m.center << endl;
    os << right << setw(20) << "Normal"     << setw(20) << " " << m.normal << endl;
    os << right << setw(20) << "Radius"     << setw(20) << " " << m.radius << endl;
    os << right << setw(20) << "DTheta"     << setw(20) << " " << m.dtheta << endl;
  }

  if ( iequals(type, "linear") ) 
  {
    os << right << setw(20) << "From"       << setw(20) << " " << m.from << endl;
    os << right << setw(20) << "To"         << setw(20) << " " << m.to << endl;
    os << right << setw(20) << "Npts"       << setw(20) << " " << m.npts << endl;
  }

  if ( iequals(type, "gauss") ) 
  {
    os << right << setw(20) << "Order"      << setw(20) << " " << m.order << endl;
    os << right << setw(20) << "Boundaries" << setw(20) << " " << m.boundaries << endl;
  }

  return os;
}
