#include "base/Global.h"
#include "config/Config.h"
#include "util/TeeBuff.h"

#include <fstream>
#include <stdio.h>
#include <boost/stacktrace.hpp>


unsigned int Log::_curr_indent = 0;

std::string Log::LogColor = RST;
bool Log::show_msg_lines = 1;
int Log::msg_lines_col = 80;

map<tuple<string,string,int>,uint> msg_count;

/**
 *
 */
Log::Log(FILE * out, bool fail,const char * file, int line, uint maxtimes, string tag, bool log1) :
  _out(out), _fail(fail), _file(file), _line(line),
  _newline(1), _maxtimes(maxtimes), _tag(tag), _log1(log1)
{ }

/**
 *
 * Redireciona as mensagens do petsc.
 *
 */

/**
 *
 * Arquivos de log, por processador.
 *
 */
std::ofstream * Log::_dfile = 0;
void Log::init( uint rank ) {
  string filename("run/log/proc_out");
  filename = filename + string("-");
  filename = filename + std::to_string(rank);
  filename = filename + string(".log");
  _dfile = new ofstream(filename);

  // Duplica os buffer do libMesh para cout e cerr
  teebuf * tb = new teebuf( _dfile->rdbuf(), cerr.rdbuf() );
  libMesh::out.rdbuf( tb );
  libMesh::err.rdbuf( tb );
}

/**
 * Retorna um stringstream onde serao armazenadas as mensagens de log.
 *
 * \param level define o tipo de mensgem:
 *
 *            0: informativo (# )
 *           -1: aviso (# Warning)
 *           -2: erro (# Error)
 *           -3: log do petsc
 *           >0: mensagem de depuração (# Debug)
 *
 * Se a flag \a _fail estiver setada, é uma falha não tolerável (Failure)
 *
*/
ostringstream & Log::get( int level )
{
  // Color only to the terminal (not to the log files)
  if (_fail) color_red();
  else if (level == 0)  color_green();
  else if (level == -1) color_yellow();
  else if (level == -2) color_red();
  else if (level == -3) color_rst();
  else color_cyan();

  // Prefix to the lines
  if (_fail) _os << "# Fail[p" << RANK <<"]: ";
  else if (level == 0) _os <<  "# Info:     ";
  else if (level == -1) _os <<  "# Warn:     ";
  else if (level == -2) _os  << "# Err:      ";
  else if (level == -3) {
    _newline = 0;
    _os << "# PETSC: ";
  }
  else _os << "# Deb[" << level << ":" << RANK << "]: ";

  if ( _curr_indent ) {
    for ( uint i=0;i<_curr_indent;i++) _os << "_";
    _os << " ";
  }

//  _os.setf( std::ios_base::scientific, std::ios_base::floatfield );
  return _os;
}


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
std::string Log::now() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
  return buf;
}

/**
 * Destrutor. Escreve no output configurado no construtor ("_out") o que
 * foi armazenado no ostringsteam "os".
 *
 * Coloca uma nova linha no final. Evita que o usuario fique sempre colocando
 * "endl" no final de cada log.
*/
Log::~Log()
{
  bool screen = true;
  if ( _log1 && RANK >= 1 ) screen = false;

  // Aborta se ja foi exibido mais do que o solicitado
  if ( _maxtimes ) {
    string f(_file);
    tuple<string,string,int> p = { f, _tag, _line };
    if ( ! msg_count.count(p) ) {
      msg_count[p] = 1;
    } else {
      uint c = msg_count.at(p);
      msg_count[p] = ( c+1 );
      if ( c >= _maxtimes ) return;
    }
  }

  if ( screen ) {
    if ( LogColor != RST ) fprintf(_out, "%s", LogColor.c_str());
    fprintf(_out, "%s", _os.str().c_str());
    if ( LogColor != RST ) fprintf(_out, "%s", string(RST).c_str());
  }

  if ( _dfile ) * _dfile << _os.str();


  if ( Log::show_msg_lines ) {
    if ( _file ) { 
      // Tenta colocar os @ alinhados - nao funcionara em multilinha
      uint ssize = _os.str().size();
      if ( ssize>80 ) ssize = 80;
      int spaces = Log::msg_lines_col - ssize;
      if ( spaces < 5 ) spaces = 5;

      if ( screen ) {
        if ( LogColor != RST ) fprintf(_out, "%s", LogColor.c_str());
        for ( int i=0;i<spaces;i++) fprintf(_out, " ");
        fprintf( _out, "@ %s", _file ); 
        if ( _line > 0 ) { fprintf(_out, ":%d", _line); }
        if ( LogColor != RST ) fprintf(_out, "%s", string(RST).c_str());
      }

      if ( _dfile ) {
        for ( int i=0;i<spaces;i++) *_dfile << " ";
        *_dfile << "@ " << _file;
        if ( _line > 0 ) *_dfile << ":" << _line;
      }
    }
  }

  if ( _dfile ) {
    if ( _newline ) *_dfile << endl;
    _dfile->flush();
  }

  if ( screen ) {
    if ( _newline ) fprintf(_out, "\n");
    fflush(_out);
  }

  if (_fail) {
    boost::stacktrace::stacktrace st;
    ostringstream os;
    os << st;

    fprintf(_out, "%s", os.str().c_str());
    if ( _dfile ) *_dfile << st;
    breakpoint();
    exit(-1);
  }
}

/**
 *
 * Acerta o indent para o escopo.
 *
 */
ScopeLog::ScopeLog(const string & scope, const string & file, int line, uint dlevel, bool log1) :
    _scope(scope), _dlevel(dlevel), _log1(log1) { 

  if ( _log1 )
    { // importante o escopo para o objeto ser distruido logo.
      dlog1(_dlevel)<<"ScopeLog @ " << _scope << " @ " << file << ":" << line;
    }
  else
    { // importante o escopo para o objeto ser distruido logo.
      dlog(_dlevel)<<"ScopeLog @ " << _scope << " @ " << file << ":" << line;
    }

  Log::indent_inc(); 
}

/**
 *
 * Volta o indent.
 *
 */
ScopeLog::~ScopeLog() {
  Log::indent_dec(); 
  { // importante o escopo para o objeto ser distruido logo.
    dlog(_dlevel)<<"~ScopeLog @ " << _scope;
  }
}

//
// Funcao dummy para ter um breakpoint intuitivo
// (interface com o gdb).
//
// Chamada do gdb:
//   gdb --ex 'set breakpoint pending on' --ex 'b MPI_Abort' --ex 'b breakpoint' --ex 'run' --args $run $progargs
//
//
void breakpoint() {}





