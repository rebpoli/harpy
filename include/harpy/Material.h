#pragma once

#include "base/Global.h" 

/**
 *
 * This is an abstract class to be the common inteface 
 * for different materials being passed around in the code.
 *
 */

class MaterialConfig;

class Material 
{
  public:
    Material( MaterialConfig & mat_conf );

    static Material * Factory( uint sid );
};
