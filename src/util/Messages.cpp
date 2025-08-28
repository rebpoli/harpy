#include "util/Messages.h"
#include "timeloop/Timestep.h"
#include "libmesh/linear_solver.h"

namespace util {

using namespace libMesh;
using namespace std;

/**
 *
 * Funcao para formatar um inteiro...
 *
 */
string fmt_i( int num, uint w, bool zerofill ) {
  ostringstream out;
  out << setw(w);
  if ( zerofill ) out << setfill('0');
  out << num;
  return out.str();
}
/**
 *
 *
 */
string fmt_d( double num ) {
  ostringstream out;
  out << num;
  return out.str();
}

/**
 *
 *
 */
string fmt_sci( double num ) {
  ostringstream out;
  out << scientific << num;
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
  ostringstream out;

  out << "Solving THERMAL time step ";
  out << setw(2) << right << t_step
    << ", time=" << fixed << setw(6)
    << setprecision(3) << setfill('0') << left
    << time <<  "...";

  return out.str();
}

/**
 *
 */
string MSG_POROEL_INIT_TIMESTEP( Timestep & ts ) {
  ostringstream out;
  out << "Solving POROELASTIC time step "
    << setw(2) << right << ts.t_step
    << ", time=" << fixed << setw(6)
    << setprecision(3) << setfill('0') << left
    << ts.time <<  " (dt:"<< ts.dt;

  return out.str();
}

/**
 *
 */
string MSG_POROEL_WRITE( Timestep & ts ) {
  UNUSED(ts);
  ostringstream out;
  out << "Writing poroelastic timestep ...";
  return out.str();
}

/**
 *
 */
string MSG_POROEL_EVAL_PROBE( string pname ) {
  ostringstream out;
  out << "Evaluating poroelastic probe '"<< pname << "' ...";
  return out.str();
}


/**
 *
 *
 */
string MSG_CONVERGED_REASON( const libMesh::TransientLinearImplicitSystem & sys ) {
  ostringstream out;
  uint nit = sys.n_linear_iterations();
  out << "LINEAR solver convergence/divergence reason: ";
  out << libMesh::Utility::enum_to_string(sys.get_linear_solver()->get_converged_reason());
  out << " after " << nit << " LINEAR iterations"; 
  return out.str();
}

} // ns
