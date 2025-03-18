
## Organization of source files - include/ and src/

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

## Class diagram

#### Interfaces


```mermaid
%% Interfaces
    class Timeloop { 
        BoundaryConditions bc
        Timestep ts
        solve()
    }
    BoundaryConditions <-- Timeloop

    class Solverloop {
        solve()
    }

    class Solver {
        map< var, Calculator * > calculators
        map< sid, Material * > material_by_subdomain
        ~Solver: delete materials and calculators
        project_from(Solver, vars)
        solve()
        jacobian()
        residual()
    }
    class SolverConst { A solver that returns const values.  }


    class Material {
        BoundaryConditions *bc
        +jacobian()
        +residual()
    }
    BoundaryConditions <.. Material

    class CalcEntry {
        vector[qp] val
        vector[qp] grad
    }
    class Calculator {
        map< eid, vector< CalcEntry > > entries_by_eid
        set_material( Material )
        eval( vector<Point> )
    }
    Material <.. Calculator
    CalcEntry <-- Calculator
    Solverloop <|.. SolverloopBasic
    Calculator <-- Solver

```

#### Implementation


```mermaid
classDiagram

    class Timestep {
        Controls the time.
        Register a callback to BoundaryConditions.
    }
    class BoundaryConditions { 
        Holds current BCs
        Knows the current delta_t
        Updated when Timestep is updated by callback
    }
    BoundaryConditions ..> Timestep

%% Implementations
    class SolverloopBasic {
        FEM EquationSystems
        Coupling beween solvers
    }

    class MatViscoPlastic { FEM Element }

    class MatPoroelastic {
        FEM Element 
        FE ShapeFunc
        QRule
        reinit(Elem)
    }

    class SolverloopSingleNR { Bypasses to the NR solver. } 
    class SolverNR { Multiplexes element materials. }

%% Dependencies
    Timeloop --> Timestep
    Timeloop ..> Solverloop

    SolverloopSingleNR ..> SolverNR

%% Interfaces << implementations
    Timeloop <|.. TimeloopBasic
    Material <|.. MatPoroelastic
    Material <|.. MatViscoPlastic
    Solverloop <|.. SolverloopSingleNR
    Solver <|.. SolverNR
    Solver <|.. SolverConst

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

```

##  Algorithms

#### Timeloop::Timeloop()

<pre>
    Create Timestep structure
    Add TS callback to BoundaryConditions (created during config)
</pre>

#### Timeloop::Timeloop()
<pre>
    while ( true ) 
        Solverloop.solve()
        Solverloop.export()
        Timestep.next()
        if ( Timestep > max ) : break
</pre>

#### Calculator::eval( Solver S, Elem TE, var )
Evaluates the values at the quadrature points of the target materials.
Register in the entries_by_eid. Only in the processor that owns the element.
<pre>
    Mat TM = get_mat( TE )                 // Target material (shape funcs)
    qpxyz = TM.get_qpxyz()                 // The points of each quadrature point in the target
    foreach (Point pt) in (qpxyz)
        entries = entries_by_eid[TM.eid]
        Elem SE = S.find_elem( pt )        // Colective task! All processors in sync
        Mat SM = S.get_mat( SE )           // Source material (shape funcs)
        if ( curr_proc )
            CalcEntry ce = SM.calc()       // Now only the right processor does the calc
            entries.push( ce )
</pre>
        
#### Solver::project_from(Solver S, vars)
Creates a fully calculated structure in each integration point of the target.
Remember: if we need to find elements by point, this is a collective task (need to iterate in all
elements of the mesh in all processors, in sync).
<pre>
    foreach (Elem E) in (this)
        foreach (var) in (vars)
            calc = calculators[var]
            calc.eval( S, E, var )  // CalcEntry holds the information at the list of points
</pre>

#### Solverloop::solve()
Integrates many solvers.
This workflow should be implemented in the child classes.
<pre>
    Solver S1, S2                   // Instantiate the desired types
    while ( true )
        S1.project_from( S2, vars )
        S1.solve()

        S2.project_from( S1, vars )
        S2.solve()

        if (converged) : break
        if (maxit) : break

        // export intermediate results for debugging
</pre>

#### Solver::get_mat( Elem E )   && get_mat( Elem E, Side S )
Retrieves the material for the element.
Creates a new one if does not exist (lazy worker).
The material holds the FE shape functions and quadrature points for integration.
It is responsible for building the element matrix and RHS.

Two versions: one for the element, another for a side of the element.
The FE structures are the same, but the dimensions are different.

<pre>
    sid = E.subdomain()
    if not material_by_subdomain[sid] :
        material_by_subdomain[sid] = Material::Factory(sid)

    material_by_subdomain[sid].reinit(E)
</pre>

#### Solver::Solver
Creates all materials of the local processor upfront.

<pre>
    for (Elem E) in (this.local) : get_mat(E)   // Create all materials.  
</pre>
    
#### **static** Material::factory( sid )
Multiplexes the material from the configuration.
Instantiates the right material for the element.
Should only be called if it hasnt been created befor (see Solver::get_mat)
<pre>
    matid = get_mat_id(sid) 
    if matid == VISCOPLASTIC :
        return new MatViscoPlastic( E, BC ) 
    if matid == POROELASTIC :
        return new MatViscoPlastic( E, BC )
</pre>


#### Material::Material( EquationSystem, System )
The constructor should be able to query the configuration structure, build the 
FE structures for the material, the element matrix and RHS structure etc.

This constructor fetches all the information needed from EquationSystems,
like the volume of parameters (permeability, porosity etc). Thhe best is to
map these parameters into an input_system to be queried in the jaobian and
residual functions. 

#### Material::reinit( Elem E )
Reinitializes the FE shape function and quadrature for the element.

#### Material::jacobian( solution, K )  &&  jacobian_bc( solution, K )
Fills the global matrix K with the contributions of the current element.
**jacobian_bc** adds the boundary conditions equations to the matrix.

#### Solver::jacobian( solution, K )
Multiplexes the material and calls the jacobian in the material.

<pre>
    // Continua
    for (Elem E) in (local.this)
        M = get_mat(E)
        M.jacobian( solution, K )

    // Boundary conditions
    for (Elem E) in (local.this)
    for (Side S) in (E)
        if ( not has_bc ) continue
        M = get_mat( E, S )
        M.jacobian_bc( solution, K )
</pre>
