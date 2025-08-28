#include "timeloop/Timestep.h"
#include <iomanip>

#include "config/ModelConfig.h"
#include "config/TimestepConfig.h"

namespace timeloop {

using config::MODEL;
using config::BCConfig;

/**
 *
 */
Timestep::Timestep() :
       config(MODEL->timestep),
       t_step(0), dt(config.dt0), 
       time(-1) , cut_next(0)
{

  // Build tsqueue
  set<double> tsset;
  BCConfig & bcc = MODEL->boundary_config;
  for ( auto & [k,v] : bcc.entry_by_time ) 
  if ( k > 0 ) tsset.insert(k);

  for ( double t : tsset ) 
  {
    dlog(1) << "--- build tsqueue: " << t;
    tsqueue.push(t);
  }
}


/**
 * Returns true if we reached the last timestep of the run.
 */
bool Timestep::test_end() const {
  if ( t_step > config.max_steps ) return true;

  if ( config.t_max > 0 )
  if ( time > config.t_max ) return true;

  return false;
}

/**
 *  Advances the timestep.
 */
double Timestep::next() {
  SCOPELOG1(1);

  if ( cut_next )  cut_next--;

  // This means that we should not increment t_step and we should cut dt in half
  if ( cut_next == 4 ) 
  {
    time -= dt; 
    dt /= config.dtk; 
    time += dt; 
    if ( time < 0 ) { dt = config.dt_min; time=dt; }

    dlog(1) << "Cutting timestep (time: " << time << ", dt= " << dt << ")";
    if ( dt < config.dt_min ) flog << "Cannot continue: dt(" << dt << ")<dt_min(" << config.dt_min << ")";
    return dt;
  }

  t_step++;

  // Primeiro timestep após a inicializacao, roda com dt0
  if ( time < 0 ) 
  { 
    dt = config.dt0; 
    time = dt;
    return dt; 
  }

  // Controle normal dos timesteps ( dt = progressao geometrica )
  if ( ! cut_next ) dt = dt * config.dtk; // Increment dt only if cut_next=0
  // Constrain dt
  if ( dt > config.dt_max ) dt = config.dt_max;
  if ( dt < config.dt0 ) dt = config.dt0; // This happens after using an entry from tsqueue, near the previous TS

  double prev_time = time;

  // Verifica se nao estamos bypassando um tsvec
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
    if ( _dt < config.dt_min ) 
    {
      tsqueue.pop();
      continue; 
    } 

    // Update only if _dt is significantly larger than the current ts.
    if ( _dt < dt ) 
    {
      dt = _dt ; tsqueue.pop(); break; 
    }

    if ( _dt < config.dt_min )
      wlog << "Timestep less than minimum (_dt="<<_dt<<" < dt_min="<<config.dt_min << ")";

    break; // _dt >= dt ... good to go.
  }

  time += dt;

  dlog(1) << "Timestep::next: t_step=" << t_step << ", dt="<<dt<<", time="<< scientific << setw(20) << time << ")";
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

} // ns
