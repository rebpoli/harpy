// Isso eh um include guard
#pragma once

// Para melhorar a eficiencia, zerar esse
// parametro em tempo de compilacao. O compilador
// ira remover todas as chamadas de mensagens de debug.
//
// O PARAMETRO DEBUG é usado tanto da lincagem quanto na
// compilacao. Importante setálo no ambiente e nao aqui dentro
// senao pode dar erro de lincagem por inconsistencia
//

#define UNUSED(expr) do { (void)(expr); } while (0)

// NOTA: O NUMERO 1E-8 EH UM NUMERO MAGICO.
// FOI UTILIZADO BASEADO NA EXPERIENCIA COM O LIBMESH E O GMSH, QUE VEZ POR OUTRA
// GERAM COORDENADAS NO MESMO PLANO, MAS COM ATÉ 1E-10 DE IMPRECISAO NUMERICA

/// Compara dois libMesh::Point para ver se estao no mesmo local 
#define PT_CMP(p1,p2) ( (p1-p2).norm() < 1E-8 )
/// Compara dois doubles utilizando um epsilon de tolerancia
#define DCMP(v1,v2) ( std::abs(v1-v2) < 1E-8 )
/// LESS THAN de dois doubles utilizando um epsilon de tolerancia
#define DLT(v1,v2) ( v1 - v2 + 1E-8 < 0)
/// GREATER THAN de dois doubles utilizando um epsilon de tolerancia
#define DGT(v1,v2) ( v1 - v2 - 1E-8 > 0)
/// LESS OR EQUAL de dois doubles utilizando um epsilon de tolerancia
#define DLEQ(v1,v2) ( DCMP(v1,v2) || DLT(v1,v2) )
/// GREATER THAN OR EQUAL de dois doubles utilizando um epsilon de tolerancia
#define DGEQ(v1,v2) ( DCMP(v1,v2) || DGT(v1,v2) )



#include <string>
using namespace std;
extern double PI;
typedef unsigned int uint;
typedef short unsigned int suint;
typedef short int sint;

extern uint debug_level ;
extern string CHIMAS_HOME;
extern uint RANK ;

#define SINGLEPROC (RANK<1)

#include "harpy/Log.h"

class Math {
  public:
    inline static double kronecker_delta(unsigned int i, unsigned int j) { return i == j ? 1. : 0.; }
};


class Tester;

// Universal constants
const double GRAVITY_FORCE = 9.80665;

// Enables assert tests (can harm performance)
#define ENABLE_ASSERTS
#include "harpy/Assert.h"
