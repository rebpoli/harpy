#include "config/ConfigValidate.h"
#include "config/Config.h"
#include "util/Messages.h"
#include "util/File.h"

#include <fstream>

using namespace rapidjson;

/**
 *
 *
 */
ConfigValidate::ConfigValidate( Config & cfg_ ) :
  cfg( cfg_ ) 
{
  string json_fn = CHIMAS_HOME+"/share/etc/grammar.json";
//  Processos nao inicializados ainda, vai imprimir em todos eles.
//  ilog1 << "Lendo gramatica json @ '" << json_fn << "'...";
  RJDoc gram_doc;
  read_json( gram_doc, json_fn );

  auto & rj = CFG._rj;
  validate( rj, gram_doc );
}

/**
 *
 *
 */
string ConfigValidate::str( const Value & gram, string key )
{
  if ( ! gram.IsObject() ) flog << context<<" entrada da gramatica mal formatada! Deve ser um dicionario.";
  if ( ! gram.HasMember( key.c_str() ) ) { help(gram) ; flog << "Nao achei '" << key << "' na gramatica @ " << context; }
  return gram[ key.c_str() ].GetString();
}
/**
 *
 *
 */
double ConfigValidate::dbl( const Value & gram, string key )
{
  if ( ! gram.IsObject() ) flog << context<<" Entrada da gramatica mal formatada! Deve ser um dicionario.";
  if ( ! gram.HasMember( key.c_str() ) ) { help(gram) ; flog << "Nao achei '" << key << "' na gramatica @ " << context; }
  return gram[ key.c_str() ].GetDouble();
}

/**
 *
 *
 *
 */
void ConfigValidate::validate_type( RJValue & trg, RJValue & gram ) 
{
  string type = get_type( gram );
  if ( type == "dict" )   if ( ! trg.IsObject() )                  { help(gram) ; flog << context << " Tipo 'dict' não respeitado.";}
  if ( type == "double" ) if ( ! trg.IsDouble() and !trg.IsInt() ) { help(gram) ; flog << context << " Tipo 'double' não respeitado.";}
  if ( type == "uint" )   if ( !trg.IsUint() )                     { help(gram) ; flog << context << " Tipo 'uint' não respeitado.";}
  if ( type == "string" ) if ( ! trg.IsString() )                  { help(gram) ; flog << context << " Tipo 'string' não respeitado.";}
  if ( type == "list" )   if ( ! trg.IsArray() )                   { help(gram) ; flog << context << " Tipo 'list' não respeitado.";}
  if ( type == "point" )   if ( ! trg.IsArray() )                   { help(gram) ; flog << context << " Tipo 'point' não respeitado.";}
  if ( type == "bool" )   if ( ! trg.IsBool() )                    { help(gram) ; flog << context << " Tipo 'bool' não respeitado.";}

  // melhorar isso: o valor de uma variavvel pode ser um numero ou uma string apontando para 
  // um escalar.
  if ( type == "var_value" )   if ( ! trg.IsNumber() and !trg.IsString()) { help(gram) ; flog << context << " Tipo 'var_value' não respeitado.";}
}

/**
 *
 *
 *
 */
void ConfigValidate::validate_string( RJValue & trg, RJValue & gram ) 
{
  if ( ! gram.IsObject() ) return;
  if ( ! trg.IsString() ) return;

  string key = trg.GetString();
  if ( ! check_namelist(key,gram) ) {
    help(gram); 
    flog << context << " String '" << key << "' nao eh valida - nao esta na namelist.";
  }
}

/**
 *
 *
 *
 */
void ConfigValidate::validate_file( RJValue & trg, RJValue & gram ) 
{
  if ( ! gram.IsObject() ) return;
  if ( ! trg.IsString() ) return;

  string fn = trg.GetString();
  if ( gram.HasMember("_file_exists") ) {
    if ( gram.HasMember("_dir") )
        fn = gram["_dir"].GetString() + string("/") + fn;

    if ( ! file_exists( fn ) ) {
      help(gram); flog << context << " Arquivo '" << fn << "' nao abriu para leitura;";
    }
  }
}

/**
 *
 *
 */
void ConfigValidate::help( RJValue & gram )
{
  if ( ! gram.IsObject() ) return;
  if ( gram.HasMember("_help") ) elog << gram["_help"].GetString();
}

/**
 *
 *
 *
 */
void ConfigValidate::validate_range( RJValue & trg, RJValue & gram ) 
{
  if ( ! trg.IsNumber() ) return;
  if ( ! gram.IsObject() ) return;

  double val = trg.GetDouble();

  if ( gram.HasMember("_min") ) {
    double min = dbl(gram, "_min");
    if ( val < min ) { help(gram) ; flog << context << " Valor menor que o mínimo ( " << val << " < " << min << " )";}
  }

  if ( gram.HasMember("_max") ) {
    double max = dbl(gram, "_max");
    if ( val > max ) { help(gram) ; flog << context << " Valor maior que o máximo ( " << val << " > " << max << " )";}
  }
}

/**
 *
 *
 */
void ConfigValidate::validate_point( RJValue & trg, RJValue & gram )
{
  if ( ! trg.IsArray() ) return;
  string type = get_type(gram);
  if ( type != "point" ) return;

  uint l = 0; // elems
  for ( auto& v : trg.GetArray() ) {
    if ( l++ > 3 ) break;
    if ( ! v.IsNumber() ) flog << context << " Cada componente de um ponto deve ser um numero.";
  }
  if ( l != 3 ) flog << context << " Um ponto no espaço tridimensional deve ter 3 componentes.";
}

/**
 *
 *
 */
void ConfigValidate::validate_matrix( RJValue & trg, RJValue & gram )
{
  if ( ! trg.IsArray() ) return;
  string type = get_type(gram);
  if ( type != "matrix" ) return;

  uint l = 0; // linhas
  uint c = 0; // colunas
  for ( auto& v : trg.GetArray() ) {
    l++;
    if ( v.IsArray() ) {
      uint _c=0;
      for ( auto &vv : v.GetArray() ) {
        if ( ! vv.IsNumber() ) 
          flog << context << " Matriz tem elementos nao numericos?";
        _c++;
      }

      if ( ! c ) c = _c;
      if ( c != _c ) flog << context << " Matriz nao tem linhas e colunas do mesmo tamanho! Revisar.";
    }
  }

  dlog(3) << "Lida matriz com dimensoes: (l,c) = (" << l << "," << c << ").";
  if ( gram.IsObject() ) return;
  if ( gram.HasMember("_dim") ) {
    auto  a = gram["_dim"].GetArray();
    uint _l = a[0].GetInt();
    uint _c = a[1].GetInt();
    string dim = "(" + fmt_i(_l,0) + "," + fmt_i(_c,0) + ")";
    if ( (_l != l) || (_c != c) ) flog << context << " Matriz nao tem dimensoes adequadas ("<<l<<","<<c<<") != " << "["<<dim<<"].";
  }
}

/**
 *
 * Funcao recursiva para validacao do json
 *
 */
void ConfigValidate::validate( const Value & trg, const Value & gram )
{

  // Le na gramatica o que esse cara deve ser
  string type = get_type(gram);
  dlog(3) << "Validando '"<< context << "' ... type:" << type;
  
  validate_type ( trg, gram );
  validate_string( trg, gram );
  validate_range( trg, gram );
  validate_file ( trg, gram );
  validate_matrix ( trg, gram );
  validate_point ( trg, gram );

  // Estamos diante de um objeto? ==> itera e aprofunda (é um dict)
  if ( trg.IsObject() )
  for (Value::ConstMemberIterator itr = trg.MemberBegin(); itr != trg.MemberEnd(); ++itr) 
  {
    string tkey = itr->name.GetString();
    auto CS = this->context.scope( tkey );
    const Value * gval = scan( gram, tkey );
    if ( gval ) validate( itr->value, *gval );
  }

  // Estamos diante de um array? ==> itera e aprofunda 
  if ( type == "list" ) {
    uint i=0;
    for (auto & V : trg.GetArray() )
    {
      auto CS = this->context.scope(fmt_i(i,0));
      const Value * gval = scan( gram, i );
      if ( gval ) validate( V, *gval );
      ++i;
    }
  }
  
}

string ConfigValidate::get_type( RJValue & gram )
{
  string ret;
  if ( gram.IsObject() ) ret = str( gram, "_type" );
  else if ( ! gram.IsString() ) flog << "[context]" << " Nao consegui determinar o tipo.";
  else ret = gram.GetString();
  
  set<string> types = {"dict","double","uint","string","list","named_dict","matrix","bool","point", "var_value"};
  if ( ! types.count(ret) ) flog << "[context]" << " Tipo '"<< ret <<"' nao existe.";
  return ret;
}

/**
 *
 *
 */
const Value * ConfigValidate::scan( RJValue & gram, uint i )
{
  UNUSED(i);
  string type = get_type(gram);
  if ( type != "list" ) flog << context << " Essa funcao deveria ser chamada apenas para listas? Habemos bug?";

  dlog(3) << "Looking for key '" << context << "' ...";

  // Nao tem descricao de entrada -- skip.
  if ( ! gram.IsObject() ) return 0;

  if ( gram.HasMember( "_entry" ) ) return & gram["_entry"];
  flog << context << " Não encontrei '_entry' na gramatica para a lista.";
  return 0;
}

/**
 *
 * Retorna verdadeiro se a gramatica tem o namelist e a chave esta la.
 *
 */
bool ConfigValidate::check_namelist( const string & key, RJValue & gram )
{
  if ( ! gram.IsObject() ) return true;
  if ( ! gram.HasMember( "_namelist" ) ) return true;

  if (!gram["_namelist"].IsArray() ) flog << "[context] _namelist deve ser uma lista!";

  for ( auto & v : gram["_namelist"].GetArray() )
    if (v.GetString() == key ) 
      return true;

  return false;
}

/**
 *
 *
 */
const Value * ConfigValidate::scan( RJValue & gram, string & key )
{
  string type = get_type(gram);

  dlog(3) << "Looking for key '" << context << "' ...";

  if ( gram.HasMember( key.c_str() ) ) return & gram[key.c_str()];

  if ( check_namelist( key, gram ) )
  if ( gram.HasMember( "_entry" ) ) 
    return & gram["_entry"];
  
  if ( type == "named_dict" ) { help(gram) ; flog << context << " Nao encontrei chave '" << key  << "' na gramatica. O 'named_dict' só aceita chaves conhecidas.";}
  
// Nao tem descricao de entrada -- skip.
//  if ( ! gram.IsObject() ) return 0;

  if ( gram.HasMember( "_entry" ) ) return & gram["_entry"];
  flog << context << " Não encontrei '_entry' na gramatica.";
  return 0;
}


/**
 *
 *
 */
ostream& operator<<(ostream& os, const ConfigValidate::Context & m)
{ os << "[" << m.str() << "]"; return os; }

