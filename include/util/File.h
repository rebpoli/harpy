#pragma once

#include "base/Global.h"

#include <boost/filesystem.hpp>

#include <fstream>

/**
 *
 */
inline bool file_exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

/**
 *
 */
inline bool dir_exists (const std::string& name) {
    return boost::filesystem::is_directory(name);
}

/**
 *
 */
inline string file_canonical( const string & fn ) {
  namespace fs = boost::filesystem;
  if ( ! file_exists(fn) )
    flog << "[file_canonical] File '"<< fn <<"' nao existe!";
  auto p = fs::path(fn);
  p = fs::canonical(p);
  return p.native();

}

