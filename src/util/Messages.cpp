#include "util/Messages.h"
#include "config/TimeloopConfig.h"
#include "libmesh/linear_solver.h"

using namespace libMesh;

/**
 *
 * Funcao para formatar um inteiro...
 *
 */
string fmt_i( int num, uint w, bool zerofill ) {
  std::ostringstream out;
  out << std::setw(w);
  if ( zerofill ) out << std::setfill('0');
  out << num;
  return out.str();
}
/**
 *
 *
 */
string fmt_d( double num ) {
  std::ostringstream out;
  out << num;
  return out.str();
}

/**
 *
 *
 */
string fmt_sci( double num ) {
  std::ostringstream out;
  out << std::scientific << num;
  return out.str();
}

/**
 *
 * Principais mensagens
 * 
 */

/**
 *
 */
string MSG_THERMAL_INIT_TIMESTEP( uint t_step, double time ) {
  std::ostringstream out;

  out << "Solving THERMAL time step ";
  out << std::setw(2) << std::right << t_step
    << ", time=" << std::fixed << std::setw(6)
    << std::setprecision(3) << std::setfill('0') << std::left
    << time <<  "...";

  return out.str();
}

/**
 *
 */
string MSG_POROEL_INIT_TIMESTEP( Timestep & ts ) {
  std::ostringstream out;
  out << "Solving POROELASTIC time step "
    << std::setw(2) << std::right << ts.t_step
    << ", time=" << std::fixed << std::setw(6)
    << std::setprecision(3) << std::setfill('0') << std::left
    << ts.time <<  " (dt:"<< ts.dt;

  return out.str();
}

/**
 *
 */
string MSG_POROEL_WRITE( Timestep & ts ) {
  UNUSED(ts);
  std::ostringstream out;
  out << "Writing poroelastic timestep ...";
  return out.str();
}

/**
 *
 */
string MSG_POROEL_EVAL_PROBE( string pname ) {
  std::ostringstream out;
  out << "Evaluating poroelastic probe '"<< pname << "' ...";
  return out.str();
}


/**
 *
 *
 */
string MSG_CONVERGED_REASON( const libMesh::TransientLinearImplicitSystem & sys ) {
  std::ostringstream out;
  uint nit = sys.n_linear_iterations();
  out << "LINEAR solver convergence/divergence reason: ";
  out << libMesh::Utility::enum_to_string(sys.get_linear_solver()->get_converged_reason());
  out << " after " << nit << " LINEAR iterations"; 
  return out.str();
}

