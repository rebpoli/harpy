
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
        +project_from(Solver, vars)
        +solve()
        -jacobian()
        -residual()
    }


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

    class MatViscoPlastic { FEM Element }

    class MatPoroelastic { FEM Element }

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
    note for SolverNR "Multiplexes element materials."

    Solverloop ..> Solver
    Solver..> Material

    SolverNR ..> MatViscoPlastic
    SolverNR ..> MatPoroelastic

%%
%% classA <|-- classB    // Inheritance (B inherits from A)
%% classA <|.. classB    // Implementation (B implements A)
%% classA *-- classB     // Composition (B is part of A)
%% classA o-- classB     // Aggregation (B is used by A)
%% classA --> classB     // Association (A references B)
%% classA ..> classB     // Dependency (A depends on B)
%% classA ..|> classB    // Realization (B realizes A)
%% classA <--> classB    // Bidirectional association

```
2. **function** Merge(A, low, mid, high):
   - Create temporary arrays L and R
   - Copy A[low...mid] to L
   - Copy A[mid+1...high] to R
   - i = 0, j = 0, k = low
   - **while** i < length(L) **and** j < length(R) **do**
      1. **if** L[i] ≤ R[j] **then**
         1. A[k] = L[i]
         2. i = i + 1
      2. **else**
         1. A[k] = R[j]
         2. j = j + 1
      3. k = k + 1
