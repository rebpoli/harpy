
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
%% Interfaces
    class Timeloop { 
        solve()
    }

    class Solverloop {
        solve()
    }

    class Solver {
        solve()
        -jacobian()
        -residual()
    }
    note for Solver "Multiplexes element materials."


    class Material {
        +jacobian()
        +residual()
    }

%% Helper
    class Timestep { }

%% Implementations
    class SolverloopBasic {
        FEM EquationSystems
        Coupling beween solvers
    }

    class MatViscoPlastic {
        FEM Element
    }

%% Dependencies
    Timeloop ..> Timestep
    Timeloop ..> Solverloop

    note for SolverloopSingleNR "Bypasses to the NR solver."
    SolverloopSingleNR ..> SolverNR

%% Interfaces << implementations
    Timeloop <|.. TimeloopBasic
    Material <|.. MatPoroelastic
    Material <|.. MatViscoPlastic
    Solverloop <|.. SolverloopSingleNR
    Solver <|.. SolverNR

    Solverloop ..> Solver
    Solver..> Material

%%
%% classA <|-- classB    // Inheritance (B inherits from A)
%% classA <|.. classB    // Implementation (B implements A)
%% classA *-- classB     // Composition (B is part of A)
%% classA o-- classB     // Aggregation (B is used by A)
%% classA --> classB     // Association (A references B)
%% classA ..> classB     // Dependency (A depends on B)
%% classA ..|> classB    // Realization (B realizes A)
%% classA <--> classB    // Bidirectional association

