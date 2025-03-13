#pragma once

#include "base/Global.h"
#include "util/Json.h"

#include "libmesh/system.h"

#include <string>

using namespace std;

const double GRAVITY_FORCE = 9.80665;
class ConfigValidate;

/**
 *
 * Classe de configuração generica para fazer interface com
 * o json
 *
 */
class Config {
  public:
    Config();
    void init();

    /**
     *  Existe o namespace?
     */
    bool exists( string n ) {
      if ( !_rj.HasMember( n.c_str() ) ) return false;
      return true;
    }

    /**
     *
     */
    bool exists( string n, string p ) {
      if ( !_rj.HasMember( n.c_str() ) ) return false;
      return _rj[n.c_str()].HasMember(p.c_str()); 
    }

    /**
     *
     */
    bool str(string n, string p, string & ret, string def="") {
      if ( ! exists(n,p) ) { ret = def; return false; }

      ret = _rj[n.c_str()][p.c_str()].GetString();
      return true;
    }

    /**
     *
     */
    bool dbl(string n, string p, double & ret, double def=0) {
      if ( ! exists(n,p) ) { ret = def; return false; }
      ret = _rj[n.c_str()][p.c_str()].GetDouble();
      return true;
    }

    /**
     *
     */
    bool bln(string n, string p, bool & ret, bool def=0) {
      if ( ! exists(n,p) ) { ret = def; return false; }
      ret = _rj[n.c_str()][p.c_str()].GetBool();
      return true;
    }

    /**
     *
     */
    bool integer(string n, string p, int & ret, int def=0) {
      if ( ! exists(n,p) ) { ret = def; return false; }
      ret = _rj[n.c_str()][p.c_str()].GetInt();
      return true;
    }

    /**
     *
     */
    bool uinteger(string n, string p, uint & ret, uint def=0) {
      if ( ! exists(n,p) ) { ret = def; return false; }
      ret = _rj[n.c_str()][p.c_str()].GetUint();
      return true;
    }

    /**
     *
     */
    bool dblvec(string n, string p, vector<double> & ret, vector<double> def = {} ) {
      if ( ! exists(n,p) ) { ret=def; return false; }

      using Value = rapidjson::Value;
      const Value& vv = _rj[n.c_str()][p.c_str()];
      for ( auto & v : vv.GetArray() )
        ret.push_back( v.GetDouble() );

      return true;
    }

    /**
     *
     */
    bool strvec(string n, string p, vector<string> & ret, vector<string> def = {} ) {
      if ( ! exists(n,p) ) { ret=def; return false; }

      using Value = rapidjson::Value;
      const Value& vv = _rj[n.c_str()][p.c_str()];
      for ( auto & v : vv.GetArray() )
        ret.push_back( v.GetString() );

      return true;
    }

    /**
     *
     */
    bool strset(string n, string p, set<string> & ret, set<string> def = {} ) {
      if ( ! exists(n,p) ) { ret=def; return false; }

      using Value = rapidjson::Value;
      const Value& vv = _rj[n.c_str()][p.c_str()];
      for ( auto & v : vv.GetArray() )
        ret.insert( v.GetString() );

      return true;
    }

    /**
     *
     */
    bool strmap(string n, string p, map<string,string> & ret, map<string,string> def = {} ) {
      if ( ! exists(n,p) ) { ret=def; return false; }

      using Value = rapidjson::Value;
      Value & X = _rj[n.c_str()][p.c_str()];
      if ( ! X.IsObject() ) flog << "JSon error @ <"<<n<<"."<<p<<">: not an object.";

      for ( auto const & p : X.GetObject() )
        ret[p.name.GetString()] = p.value.GetString();

      return true;
    }

    /**
     *
     */
    string str(string n, string p) {
      if ( ! exists(n,p) ) { return ""; }
      return _rj[n.c_str()][p.c_str()].GetString();
    }

    /**
     *
     */
    double dbl(string n, string p) {
      if ( ! exists(n,p) ) { return 0; }
      return _rj[n.c_str()][p.c_str()].GetDouble();
    }

    /**
     *
     */
    static std::string getenv( const std::string & key );
    
    bool validate();

    /// The database
    rapidjson::Document _rj;

  private:
    friend ConfigValidate;
};
/**
 *  Parametro global, carregado antes que qualquer coisa no
 *  sistema.
 *
 *  Cuidado com o que se faz no contrutor do Config, deve ter
 *  poucas dependencias porque a ordem de carga das coisas pode
 *  variar.
 */
extern Config CFG;

/**
 *
 *
 */
class SolverConfig {
  public:
    SolverConfig();

    double dtol, atol, rtol, maxits;
    string ksptype, pctype;
};
