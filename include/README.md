
## Organization of soiurce files

| Directory  | Content                                                                                                                                |
|------------|----------------------------------------------------------------------------------------------------------------------------------------|
| base       | Logger and global utilities. Most source files will use these components.                                                               |
| config     | Configuration files that read JSON and provide processed (validated with default values) getters.                                       |
| harpy      | The heart of Harpy: abstract interface definitions for multiple modules.                                                               |
| material   | Material definitions sharing the harpy/Material interface.                                                                             |
| solver     | Solver implementations sharing the harpy/Solver interface.                                                                             |
| solverloop | SolverLoop implementations sharing the harpy/SolverLoop interface (coupling between multiple solvers).                                 |
| timeloop   | Timeloop implementations sharing the harpy/Timeloop interface. Handles different methods of time advancement.                          |
| util       | Peripheral code for various purposes. Includes file handling, string handling, output operators, etc.                                   |

```plantuml

@startuml

' Class definitions with descriptions
class Base {
  +Logger
  +GlobalSettings
  +init()
  +shutdown()
}

class Config {
  -jsonData
  +loadFromFile(path: string)
  +getConfigValue(key: string): any
  +validate(): boolean
}

package "Harpy Core" {
  interface Material {
    +getProperties(): Properties
    +updateState(conditions: Conditions)
  }

  interface Solver {
    +initialize()
    +solve(timeStep: double): Results
    +getConvergenceInfo(): ConvergenceInfo
  }

  interface SolverLoop {
    +addSolver(solver: Solver)
    +iterate(maxIterations: int): boolean
  }

  interface TimeLoop {
    +advance(): boolean
    +getCurrentTime(): double
    +setTimeStep(dt: double)
  }
}

' Implementations
class ConcreteTimeLoop implements TimeLoop {
  -currentTime: double
  -timeStep: double
  +advance(): boolean
  +getCurrentTime(): double
  +setTimeStep(dt: double)
}

class ConcreteSolver implements Solver {
  -material: Material
  -tolerance: double
  +initialize()
  +solve(timeStep: double): Results
  +getConvergenceInfo(): ConvergenceInfo
}

' Relationships
Base <-- Config: uses
Base <-- Harpy Core: uses
Material <-- Solver: uses
Solver <-- SolverLoop: manages
TimeLoop <-- SolverLoop: controls

@enduml
