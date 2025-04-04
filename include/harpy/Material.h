#pragma once

/**
 *
 * This is an abstract class to be the common inteface 
 * for different materials being passed around in the code.
 *
 */

class Material {
  public:
    Material();

    static Material * Factory( uint sid );
};
