
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
graph TD
    Timeloop --> SolverLoop --> Solver --> Material

classDiagram
    class TimeLoop {
        <<interface>>
        +advance() boolean
        +getCurrentTime() double
        +setTimeStep(dt) void
        +reset() void
    }
    
    class ConcreteTimeLoop {
        -currentTime double
        -timeStep double
        -initialTime double
        +advance() boolean
        +getCurrentTime() double
        +setTimeStep(dt) void
        +reset() void
    }
    
    class AdaptiveTimeLoop {
        -currentTime double
        -timeStep double
        -minTimeStep double
        -maxTimeStep double
        -errorTolerance double
        +advance() boolean
        +getCurrentTime() double
        +setTimeStep(dt) void
        +reset() void
        -calculateNextTimeStep() double
    }
