
#pragma once

/**
 *
 * This is an abstract class.
 *
 */

class Solverloop 
{

public:
  Solverloop( const Timestep & ts_ ) : ts(ts_) {}
  void solve();
  void export_results();

  const Timestep & ts;        // Owned by Timeloop

};
