
#pragma once

/**
 *
 * This is an abstract class.
 *
 */

namespace timeloop { class Timestep; }

namespace solverloop {

using timeloop::Timestep;

class Solverloop 
{
public:
  Solverloop( Timestep & ts_ ) : ts(ts_) {}
  void solve();
  void export_results();

  Timestep & ts;        // Owned by Timeloop

};

} // ns
