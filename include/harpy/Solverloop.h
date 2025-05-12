
#pragma once

/**
 *
 * This is an abstract class.
 *
 */

class Timestep;

class Solverloop 
{

public:
  Solverloop( Timestep & ts_ ) : ts(ts_) {}
  void solve();
  void export_results();

  Timestep & ts;        // Owned by Timeloop

};
