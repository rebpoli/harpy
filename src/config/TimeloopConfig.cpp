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

/**
 *
 *
 */
bool Timestep::test_end() const {
  if ( t_step > cfg.n_steps ) return true;

  if ( cfg.t_max > 0 )
    if ( time > cfg.t_max ) return true;

  return false;
}

/**
 *
 *
 */
double Timestep::next() {
  SCOPELOG1(1);
  t_step++;

  // Primeiro timestep após a inicializacao, roda com dt_min
  if ( ! time ) { 
    dt = cfg.dt_min; 
    time += dt;
    update();
    return dt; 
  }

  // Controle normal dos timesteps ( dt = progressao geometrica )
  dt = dt * cfg.dt_k;
  if ( dt > cfg.dt_max ) dt = cfg.dt_max;
  if ( dt < cfg.dt_min ) dt = cfg.dt_min;

  double prev_time = time;

  // Verifica se nao estamos bypassando um tsvec
  queue<double> & tsqueue = cfg.tsqueue;
  while(true) {
    if ( tsqueue.empty() ) break;
    double _t = tsqueue.front();
    double _dt = _t - prev_time;

    if ( _dt < 0 ) flog << "Failed to control timestep -- _dt<0? (_t="<<_t<<", prev_time="<<prev_time<<")";

    // nao estamos interessados no tempo zero, esse sempre existira
    if ( ! _t ) 
    {
      tsqueue.pop();
      continue; 
    } 

    // nao estamos interessados em passos de tempo muito pequenos
    if ( _dt < cfg.dt_min ) 
    {
      tsqueue.pop();
      continue; 
    } 

    // Update only if _dt is significantly larger than the current ts.
    if ( _dt < dt ) 
    {
      dt = _dt ; tsqueue.pop(); break; 
    }

//    Acredito que nao precisa disso, porque o usuario esta escolhendo o passo de 
//    tempo. Mas coloquei um warning caso de problema porque teoricamente o dt pode
//    ficar criticamente pequeno.
//    if ( _dt < cfg.dt_min )
//      _dt = cfg.dt_min;
    if ( _dt < cfg.dt_min )
      wlog << "Timestep less than minimum (_dt="<<_dt<<" < dt_min="<<cfg.dt_min << ")";

    break; // _dt >= dt ... good to go.
  }

  time += dt;

  dlog(1) << "Timestep::next: t_step=" << t_step << ", dt="<<dt<<", time="<< scientific << setw(20) << time << ")";
  update();
  return dt;
}

/*
 *
 */
ostream& operator<<(ostream& os, const Timestep & ts)
{
  os << "TIMESTEP " << ts.t_step << " : dt("<<ts.dt<<") time("<<ts.time<<")";
  return os;
}
