#pragma once

#include "base/Global.h"
#include "util/Json.h"

#include <list>

class Config;

class ConfigValidate {

  public:
    ConfigValidate( Config & cfg_ );
  
  private:
    void validate( RJValue & trg, RJValue & gram);
    RJValue * scan( RJValue & gram, string & key);
    RJValue * scan( RJValue & gram, uint i);

    bool check_namelist( const string & key, RJValue & gram );

    void validate_string( RJValue & trg, RJValue & gram);
    void validate_type( RJValue & trg, RJValue & gram);
    void validate_range( RJValue & trg, RJValue & gram);
    void validate_file( RJValue & trg, RJValue & gram);
    void validate_matrix( RJValue & trg, RJValue & gram);
    void validate_point( RJValue & trg, RJValue & gram);

    void help( RJValue & gram );

    string get_type( RJValue & gram );

    string str( RJValue & gram, string key );
    double dbl( RJValue & gram, string key );

    Config & cfg;

    class Context : public list<string> {
      public:
        class ContextScope {
          public: 
            ContextScope( Context & context_, const string & key ) : context(context_) { context.push_back(key); }
            ~ContextScope() { context.pop_back(); }
          private : Context & context;
        };
        ContextScope scope( const string key ) { return ContextScope( *this, key ); }
        string str() const { string ret; for ( auto & s : *this ) ret += "/" + s; return ret; }
    } context;
    friend ostream& operator<<(ostream& os, const ConfigValidate::Context & m);
};
ostream& operator<<(ostream& os, const ConfigValidate::Context & m);
