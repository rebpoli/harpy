#pragma once

#include "base/Global.h"
#include <memory>
#include <vector>

#include "libmesh/fe_base.h"
#include "libmesh/explicit_system.h"
#include "libmesh/equation_systems.h"

/** \file 
 *
 * Projeta a solução de um sistema em outro, para pegar a solução que foi
 * feita em uma malha e jogar para outra malha, mais (ou menos) granular.
 *
 */

namespace libMesh { 
  class Elem;
  class QGauss;
}

using namespace libMesh;
using namespace std;

class ElemProjection {
  public:
    /// Inicializa a engine
    ElemProjection( ExplicitSystem & sys_, unique_ptr<FEBase> & fe_, QGauss * qrule_ );

    /// Finaliza o sistema. Faz isso apos finalizar todas as projecoes desejadas,
    /// fora do loop de elemento.
    void close_system();

    /// reinicializa esta classe, seta o elemento alvo
    void reinit( Elem * e, bool fe_reinit = 0, int side = -1 );

    /// Faz a projecao propriamente dita
    void eval( uint var, const vector<double> & vals_qp );

    // Adaptadores
    void eval( string vname, const vector<double> & vals_qp );

    const libMesh::Parallel::Communicator & comm() { return sys.get_equation_systems().comm(); }

  private:
    ExplicitSystem & sys;
    unique_ptr<FEBase> & fe;
    QGauss * qrule;
    Elem * elem;
};


