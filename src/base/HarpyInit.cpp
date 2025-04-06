#include "base/Global.h"
#include "base/HarpyInit.h"
#include "config/Config.h"
#include "util/Stopwatch.h"
#include "util/File.h"

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

// Namespaces
using namespace std;
using namespace libMesh;

// Init with null pointer. HarpyInit should initialize with the right stuff.
libMesh::Parallel::Communicator * LIBMESH_COMMUNICATOR = 0;

/**
 *
 * Inicializacao global do simulador.
 * Deve estar instanciado no inicio do main e o destrutor sera chamado
 * automaticamente no final do main.
 *
 * A primeira linha do main deve ser:
 *
 *          HarpyInit init(argc, argv);
 *
 * Nota: o MyLogInit tem que rodar antes do lm_init para
 * conseguir encaminhar as mensagens do petsc. Logo, a ordem
 * da inicialização nao pode ser mudada.
 *
 */
HarpyInit::HarpyInit( int argc, char ** argv ) :
 pl(), lm_init() 
{
// Soh coloca esse codigo se, em tempo de COMPILACAO
// o CHDEBUG estiver setado. Caso afirmativo, em tempo de 
// EXECUCAO toma esse parametro como o nivel de depuracao
#ifdef CHDEBUG
  string dl = Config::getenv("CHDEBUG");
//  ilog1 << "Setting debug level: " << dl;
  if ( ! dl.length() ) debug_level = 10;
  else debug_level = std::stoi( dl );
#endif

  // Inicializa o mais elementar. Isso tem que ser antes de qualquer coisa
    // Le o config.json (nota: essa funcao nao pode contar com nada configurado, evitar prints de depuracao etc.,
    // que usam o rank do proceso, que ainda nao esta inicializado)
    CHIMAS_HOME = find_home();
//    CFG.init();  // DEPRECATED -- ALL IN config/* and MODEL
    // Aqui temos que ajustar as configuracoes de depuracao, para a inicializacao
    // fazer direito -- como elas sao definidas no json, eh importante manter a ordem
//    PENumericalConfig pec; 
//    pec.apply();
    // Inicializacao da estrutura. O paralelismo etc
    lm_init = new libMesh::LibMeshInit( argc, argv );

    // Variavel global
    LIBMESH_COMMUNICATOR = &(lm_init->comm());
    RANK = lm_init->comm().rank();

  // Cria as pastas de saida
    fs::create_directory("run");
    fs::create_directory("run/log");
    fs::create_directory("run/exo");
    fs::create_directory("run/csv");
    fs::create_directory("run/msh");

  // Inicializacoes finais
    Stopwatch::init();
    Log::init(RANK);
}

HarpyInit::~HarpyInit() { delete(lm_init); }

/**
 *
 *
 */
string HarpyInit::find_home()
{
  bool ok = true;
  string ret = Config::getenv("CHIMAS_HOME");
  if ( ! ret.length() ) {
    ok = false;
    ret = "./";
    for ( uint i = 0; i < 10 ; ++i ) {
      if ( ! file_exists( ret + ".chimas_top" ) ) { ret += "../"; continue; }
      ret = file_canonical(ret);
      ok = true;
      break;
    }
  }

  // Garante que o negocio existe
  if ( !file_exists(ret) || ! fs::is_directory(ret) || ! ok ) 
    flog << "Nao encontrei o CHIMAS_HOME. Utilize a variavel de ambiente ou rode dentro de um repositorio git do Chimas.";

  ret = file_canonical(ret);

  return ret;
}
