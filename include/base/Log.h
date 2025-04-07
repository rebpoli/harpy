// Isso eh um include guard
#ifndef __LOG_H
#define __LOG_H
/**
 * \file
 * Elementos de log.
 *
 * @param debug_level: nivel de depuracao.
 * Se maior que zero, esta em debug. 
 *
 */

#include "base/Global.h"
#include "util/TermColors.h"

#include <string>
#include <sstream>
#include <streambuf>

using namespace std;

#define __METHOD_NAME__ methodName(__PRETTY_FUNCTION__)

#define flog Log(stderr,true,__FILE__,__LINE__).get(-1)
#define ilog Log(stdout,false,__FILE__,__LINE__).get(0)
#define wlog Log(stderr,false,__FILE__,__LINE__).get(-1)
#define elog Log(stderr,false,__FILE__,__LINE__).get(-2)
#define dlog(l) if (debug_level>=l) Log(stderr,false,__FILE__,__LINE__).get(l) 
#define SCOPELOG(l) ScopeLog _______SLOG(__METHOD_NAME__, __FILE__, __LINE__, l); 
#define SCOPELOG1(l) ScopeLog _______SLOG(__METHOD_NAME__, __FILE__, __LINE__, l, true); 
#define petsc_log Log(stdout,false,__FILE__,__LINE__).get(-3)

// Apenas em um processador
#define dlog1(l) if (debug_level>=l) Log(stderr,false,__FILE__,__LINE__, 0, "", true).get(l) 
#define wlog1 Log(stderr,false,__FILE__,__LINE__, 0, "", true).get(-1)
#define elog1 Log(stderr,false,__FILE__,__LINE__, 0, "", true).get(-2)
#define ilog1 Log(stdout,false,__FILE__,__LINE__, 0, "", true).get(0)

#define dlog1_mt(l,t,tg) if (debug_level>=l) if (RANK<1) Log(stderr,false,__FILE__,__LINE__,t, tg).get(l) 

#define __PRINTLINE__ dlog(1) << "PL @ " << __FILE__ << ":" << __LINE__;

/**
 * 
 * Gerencia o logging do software.
 *
 * Com alteracoes simples pode se enviar todas as mensagens para um arquivo, por exemplo
 * ou ter apenas algum nivel de depuracao, conforme a conveniencia.
 *
 * Nao aceita concorrencia.
 *
 * Cada chamada de ilog,elog,dlog de fato cria uma instancia da classe Log, sem
 * manter sua referencia. A instancia armazena as mensagens, fazendo um processamento
 * de formatacao e destina ao terminal quando a instancia morre. Como nao se 
 * mantem a referencia da instancia, isto ocorre na linha seguinte.
 *
 * \param _os ostringstream que armazena as mensagens ate a destruicao do objeto
 * \param _out FILE * de saida. Tipicamente inicializado com stdout ou stderr, mas poderia
 *             ser um arquivo qualquer. A ideia eh ter uma funcao com argumento "string filename" para
 *             enviar o log para um arquivo qualquer. Neste caso, eh importante manter o arquivo
 *             aberto, sem fechar para cada instancia.
*/
class Log {
private:
  ostringstream _os;
  FILE * _out;   
  bool _fail;
  const char * _file;
  int _line;
  bool _newline;
  uint _maxtimes;
  string _tag;
  bool _log1;

  static unsigned int _curr_indent;

public:
  Log(FILE * out = stdout, bool fail=false,const char * file=0, int line=-1, uint maxtimes=0, string tag="", bool log1=false);
  ostringstream & get(int level);
  virtual ~Log();

  static void init( uint rank );
  static int msg_lines_col;
  static bool show_msg_lines;
  static std::string LogColor;
  static void color_rst() { LogColor = RST; }
  static void color_blue() { LogColor = std::string(KFGBLU) + KBGBLK; }
  static void color_cyan() { LogColor = std::string(KFGCYN) + KBGBLK; }
  static void color_green() { LogColor = std::string(KFGGRN) + KBGBLK; }
  static void color_red() { LogColor = std::string(KFGRED) + KBGBLK; }
  static void color_yellow() { LogColor = std::string(KFGYEL) + KBGBLK; }

  static void indent_inc() { _curr_indent += 2; };
  static void indent_dec() { _curr_indent -= 2; };
  static std::string now();

  static std::ofstream * _dfile;

  friend class ScopeLog;
};

/**
 *
 * Lida com o incremento e decremento de indentação para as mensagens do sistema
 *
 */
class ScopeLog {
  private:
      string _scope;
      uint _dlevel;
      bool _log1;
  public:
    ScopeLog( const string & scope, const string & file, int line, uint dlevel, bool log1=false ); 
    ~ScopeLog();
};

/**
 *
 * Torna as mensagens mais bonitas - especialmente as de escopo
 *
 */
inline std::string methodName(const std::string& prettyFunction)
{
  size_t colons = prettyFunction.find("::");
  size_t begin = prettyFunction.substr(0,colons).rfind(" ") + 1;
  size_t end = prettyFunction.rfind("(") - begin;

  return prettyFunction.substr(begin,end) + "()";
}

void breakpoint();

#endif
