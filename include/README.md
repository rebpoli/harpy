
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

```mermaid

classDiagram
    class Base {
        +Logger
        +GlobalSettings
        +init()
        +shutdown()
    }

    class Config {
        -jsonData
        +loadFromFile(path) string
        +getConfigValue(key) any
        +validate() boolean
    }

    class Material {
        <<interface>>
        +getProperties() Properties
        +updateState(conditions) void
    }

    class Solver {
        <<interface>>
        +initialize() void
        +solve(timeStep) Results
        +getConvergenceInfo() ConvergenceInfo
    }

    class SolverLoop {
        <<interface>>
        +addSolver(solver) void
        +iterate(maxIterations) boolean
    }

    class TimeLoop {
        <<interface>>
        +advance() boolean
        +getCurrentTime() double
        +setTimeStep(dt) void
    }

    class ConcreteTimeLoop {
        -currentTime double
        -timeStep double
        +advance() boolean
        +getCurrentTime() double
        +setTimeStep(dt) void
    }

    class ConcreteSolver {
        -material Material
        -tolerance double
        +initialize() void
        +solve(timeStep) Results
        +getConvergenceInfo() ConvergenceInfo
    }

    Config ..> Base : uses
    Base <.. "Harpy Core" : uses
    Material <.. Solver : uses
    Solver <.. SolverLoop : manages
    TimeLoop <.. SolverLoop : controls

    TimeLoop <|.. ConcreteTimeLoop : implements
    Solver <|.. ConcreteSolver : implements
