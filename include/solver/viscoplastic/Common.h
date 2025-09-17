#pragma once

#include "harpy/Global.h"

namespace solver {
namespace viscoplastic {

/** Useful enums */

enum VARIABLE { TEMPERATURE, INITIAL_TEMPERATURE, PRESSURE, INITIAL_PRESSURE };

inline string v_to_str( VARIABLE v ) 
{
  if ( v == TEMPERATURE )         return "Temperature";
  if ( v == INITIAL_TEMPERATURE ) return "Init. Temperature";
  if ( v == PRESSURE )            return "Pressure";
  if ( v == INITIAL_PRESSURE )    return "Init. Pressure";
  flog << "Should not be here.";
  return string("???");
}



}} // ns
