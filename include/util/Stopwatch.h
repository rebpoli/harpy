#pragma once

#include "harpy/Global.h"
#include "util/Messages.h"

namespace util {

/*
 * Cronometro para perfilagem de desempenho.
 * Na criacao ele recebe um nome e grava o horario de inicio.
 * Quando eh destruido, imprime o tempo que passou.
 *
 * Uso dentro de um escopo (loop ou funcao, por exemplo):
 *
 * while (xxxx) {
 *    StopWatch sw("Iteracao");
 *    . . .
 * }
 *
 * Neste caso, ao final de cada iteracao, sera impresso o tempo que se passou.
 *
 *
 * @param _start o horario inicial, na criacao do objeto
 * @param _name nome para identificar a instancia.
 *
 */

class Stopwatch {
  private:
    clock_t _start;
    string _name;
    Stopwatch * _parent;
    double _dt;
    bool _lap;
    double * _out_time;
  public:
    bool info_log; // if 1, the object also logs as info.

    Stopwatch( Stopwatch * parent ) : Stopwatch("", false) { _parent = parent; };
    Stopwatch(string name, bool lap=false, double * out_time=0);
    ~Stopwatch();

    static string csv_fn() { return string("run/csv/profile-") + fmt_i(RANK) + string(".csv"); }
    static void init();

    double _st_elapsed;
    double delta_elapsed_ms();
    double elapsed_ms();
    double elapsed_s();
    Stopwatch lap( );

    static uint _debug_level;
};

} // ns
