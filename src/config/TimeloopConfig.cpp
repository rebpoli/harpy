#include "config/Config.h"
#include "config/TimeloopConfig.h"
#include "util/OutputOperators.h"

using namespace rapidjson;

/**
 *
 *
 */
TimeloopConfig::TimeloopConfig() {
  dt = 0;
  CFG.uinteger("timeloop","n_steps",  n_steps, 10);
  CFG.dbl    ("timeloop","dt_k",     dt_k,    1);
  CFG.dbl    ("timeloop","dt_max",   dt_max,  1000);
  CFG.dbl    ("timeloop","t_max",   t_max,  -1);
  CFG.dbl    ("timeloop","dt_min",   dt_min,  1e-5);

  CFG.dblvec( "timeloop", "force_ts", force_ts );

  build_tsvec();
}

/**
 * Carrega timesteps configurados nas condicoes de contorno 
 * importante ter eles na lista de timesteps
 */
void TimeloopConfig::build_tsvec()
{
  set<double> tsset;

  vector<string> sysnames = {"poroelastic", "thermal" };
  for ( string & sname : sysnames ) {
    if ( ! CFG.exists( sname, "boundary_conditions" ) ) continue;
    const Value& dbc_config = CFG._rj[sname.c_str()]["boundary_conditions"];
    for ( auto & V : dbc_config.GetArray() ) {
      if ( ! V.HasMember("time") ) continue;
      double t = V["time"].GetDouble();
      if ( t > 0 ) // evita ts inicial (-1)
      tsset.insert(t);
    }
  }

  // Push the forced timesteps into the queue
  for ( double t : force_ts )
    tsset.insert(t);

  for ( double t : tsset )
    tsqueue.push(t);
}

