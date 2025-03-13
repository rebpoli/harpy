#include "util/Stopwatch.h"
#include "util/CsvFile.h"

#include "libmesh/perf_log.h"

#include <iostream>
#include <ctime>
#include <iostream>     // std::cout, std::fixed
#include <iomanip>


uint Stopwatch::_debug_level = 5;
libMesh::PerfLog perf ("GLOBAL_PROFILE");

/*
  Criacao do objeto. Atribui o nome e registra o horario de inicio.
  
 @param name nome para identificar a instancia
*/
Stopwatch::Stopwatch(string name, bool lap, double * out_time) : 
  _start(0), _name(name), _parent(0), _dt(0), _lap(lap), _out_time(out_time), _st_elapsed(0)
{
  dlog(_debug_level) << "[" << Log::now() << "] Stopwatch \"" << _name << "\": start... " << (lap ? "(lap)" : "") ;
  _start = clock(); 
  _st_elapsed = _start;
  perf.push(_name);
}

Stopwatch Stopwatch::lap( ) { 
  return Stopwatch( this ); 
}

/**
 *
 *
 */
void Stopwatch::init() { 
  CsvFile ofile(csv_fn(),"\t", false);
  ofile << "Tag" << "Elapsed (ms)" << endrow;
}

/**
 *
 * Retorna o tempo percorrido em ms ate o momento
 *
 */
double Stopwatch::elapsed_ms() {
  double ms = (double) (clock()-_start) / CLOCKS_PER_SEC * 1000;
  return ms;
}
double Stopwatch::delta_elapsed_ms() {
  double now = clock();
  double ms = (double) (now-_st_elapsed) / CLOCKS_PER_SEC * 1000;
  _st_elapsed = now;
  return ms;
}

/**
 *
 * Retorna o tempo percorrido em s ate o momento
 *
 */
double Stopwatch::elapsed_s() {
  double s = (double) (clock()-_start) / CLOCKS_PER_SEC;
  return s;
}

/*
  Finaliza a medicao, impriminto o tempo que passou em milisegundos.
  A mensagem eh impressa como informacao.
*/
Stopwatch::~Stopwatch() {
  string now = Log::now();
  double ms = (double) (clock()-_start) / CLOCKS_PER_SEC * 1000;
  if ( _parent ) {
    _parent->_dt += ms;
    dlog(_debug_level) << "SW/" << _parent->_name << " (lap): " << ms;
  } else {
    if ( _lap ) ms = _dt;
    double s = ms / 1000;
    dlog(_debug_level) << "["<<now<<"] Stopwatch \"" << _name <<"\": " << std::setprecision(6) << s << " s";
  }

  CsvFile ofile(csv_fn());
  ofile << _name << ms << endrow;

  // Callback - coloca o tempo do SW em uma variavel externa
  if ( _out_time ) *_out_time = ms;
  perf.pop(_name);
}
