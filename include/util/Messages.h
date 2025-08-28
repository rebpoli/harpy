#pragma once

#include "harpy/Global.h"

#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/transient_system.h"

#include <iomanip>

/**
 * Utilities, como conversoes, etc.
 *
 */

namespace timeloop { class Timestep; }

namespace util {
using timeloop::Timestep;

string fmt_i( int num, uint w=4, bool zerofill=true );
string fmt_d( double num );
string fmt_sci( double num );

/*
 * MENSAGENS DE FLUXO   
 */
string MSG_THERMAL_INIT_TIMESTEP( uint t_step, double time );
string MSG_POROEL_INIT_TIMESTEP( Timestep & ts );
string MSG_POROEL_WRITE( Timestep & ts );
string MSG_POROEL_EVAL_PROBE( string pname );
string MSG_CONVERGED_REASON( const libMesh::TransientLinearImplicitSystem & sys );

template<class T>
inline const char * to_cstr(const T & v) {
  // Por algum motivo esta linha é necessaria - não podemos retornar direto (??)
  const char *cs = to_str(v).c_str(); 
  return cs;
}

template<typename F>
inline F trunc_decs(const F& f, int decs)
{
  F ff = f * pow(10,decs);
  ff = trunc(ff);
  ff /= pow(10,decs);
  if ( ff == 0 ) ff = 0; // avoid "-0"
  return ff;
}

} // ns
