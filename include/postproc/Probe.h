
#pragma once

#include "base/Global.h"
#include "util/CsvFile.h"
#include "libmesh/system.h"

#include "config/InoutConfig.h"
#include "harpy/Solver.h"

using namespace libMesh;

/**
 *
 *  Stores the collection of probes for this model
 *
 */

class Probe
{
public:
  Probe( ProbeConfig & config );
  virtual ~Probe() = default;

  virtual void print( ostream & os ) const { UNUSED(os); };
  static Probe * Factory( vector<Probe *> & ret, ProbeConfig & config );

  virtual void eval( Solver & solver );

  virtual bool is_gauss() { return false; }

  string name, filename;
  string header;

  vector<Point> points;
};


/**
 *
 */
struct RadialProbe : public Probe
{ 
  RadialProbe( ProbeConfig & config ); 
  ~RadialProbe() = default;

  void print( ostream & os ) const;
};

/**
 *
 */
struct LinearProbe : public Probe
{
  LinearProbe( ProbeConfig & config ); 
  ~LinearProbe() = default;

  void print( ostream & os ) const;
};

/**
 *
 */
struct GaussProbe : public Probe
{
  GaussProbe( ProbeConfig & config );
  ~GaussProbe() = default;
//  void eval( System & sys, vector<string> vars = {} );

  void print( ostream & os ) const;
  virtual bool is_gauss() { return true; }

  vector<string> boundaries;
  libMesh::Order order;
};

/**
 *  Collection of probes. Constructs the structure from the configuration.
 */
struct ProbeCol : public vector<Probe *>
{ 
  ProbeCol(); 
  ~ProbeCol() { for ( auto p : *this ) delete(p); this->clear(); }
};

/**
 *
 *
 */
ostream& operator<<(ostream& os, const Probe & m);
ostream& operator<<(ostream& os, const ProbeCol & m);
