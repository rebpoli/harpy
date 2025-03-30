cmake_minimum_required(VERSION 3.28)

project( harpy )

# if(DEFINED CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose build type")
# else()
#   set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose build type")
# endif()

#set(CMAKE_BUILD_TYPE Release)

message(-----------------)
message(METHOD: "$ENV{METHOD}")

set(CMAKE_BUILD_TYPE Debug)

# LIBMESH
set(LIBMESH_DIR $ENV{LIBMESH_DIR})
if ( "${LIBMESH_DIR}" STREQUAL "" )
   set (LIBMESH_DIR "/usr/local" )
endif()
set(LIBMESH_CONFIG "${LIBMESH_DIR}/bin/libmesh-config")

# Remover isso para ficar mais rapido . Nao usar a flag DEBUG, que eh usada
# pelo libmesh e pode gerar conflitos
add_definitions(-DCHDEBUG)
add_definitions(-DBOOST_ALLOW_DEPRECATED_HEADERS)

# Compiladores MPI
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

# execute_process(COMMAND ${LIBMESH_CONFIG} --cc  OUTPUT_VARIABLE CMAKE_C_COMPILER OUTPUT_STRIP_TRAILING_WHITESPACE)
# execute_process(COMMAND ${LIBMESH_CONFIG} --cxx OUTPUT_VARIABLE CMAKE_CXX_COMPILER OUTPUT_STRIP_TRAILING_WHITESPACE)

# gmsh

find_library(BOOST_FS boost_filesystem)
find_library(GMSH_LIB gmsh)
find_path(GMSH_INC gmsh.h)
include_directories(${GMSH_INC})

# libMesh

# Get the libMesh's flags, includes, libs ...
execute_process(COMMAND ${LIBMESH_CONFIG} --cppflags OUTPUT_VARIABLE LIBMESH_CPPFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LIBMESH_CONFIG} --cxxflags OUTPUT_VARIABLE LIBMESH_CXXFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LIBMESH_CONFIG} --include OUTPUT_VARIABLE LIBMESH_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${LIBMESH_CONFIG} --libs OUTPUT_VARIABLE LIBMESH_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND ${LIBMESH_CONFIG} --version OUTPUT_VARIABLE LIBMESH_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

# libmesh compiler flags
set( CMAKE_CXX_FLAGS " -std=gnu++17 -D BOOST_STACKTRACE_LINK ${CMAKE_CXX_FLAGS} ${LIBMESH_CPPFLAGS} ${LIBMESH_CXXFLAGS} ${LIBMESH_INCLUDE}" )

message(CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS})
message(LIBMESH_VERSION: ${LIBMESH_VERSION})
message(LIBMESH_CPPFLAGS: ${LIBMESH_CPPFLAGS})
message(LIBMESH_CXXFLAGS: ${LIBMESH_CXXFLAGS})
message(LIBMESH_LIBS: ${LIBMESH_LIBS})

# Headers
include_directories(include)

# Compilar arquivos-fonte dentro da biblioteca estática libModules.a
ADD_LIBRARY(Modules
  src/harpy/Timeloop.cpp
  src/harpy/Timestep.cpp
  src/harpy/MeshInit.cpp
  src/config/TimeloopConfig.cpp
  src/config/MeshConfig.cpp
  src/config/ConfigValidate.cpp
  src/config/Config.cpp
  src/config/BCConfig.cpp
  src/config/ModelConfig.cpp
  src/config/MaterialConfig.cpp
  src/config/SolverConfig.cpp
  src/config/reader/ModelReader.cpp
  src/config/reader/SolverReader.cpp
  src/config/reader/SpatialDataReader.cpp
  src/util/Json.cpp
  src/util/Stopwatch.cpp
  src/util/OutputOperators.cpp
  src/util/TeeBuff.cpp
  src/util/Messages.cpp
  src/base/Global.cpp
  src/base/HarpyInit.cpp
  src/base/Log.cpp
  src/base/PetscLogger.cpp
  src/main/harpy.cpp
  src/solver/BC.cpp
  src/solver/SolverTHM.cpp
  src/solverloop/SolverloopTHM.cpp



#   src/StringUtil.cpp
#   src/Solver.cpp 
#   src/Stress.cpp 
#   src/StressInvariants.cpp 
#   src/Poroelastic.cpp 
#   src/Probe.cpp 
#   src/RadialAnalytic.cpp 
#   src/Terzaghi.cpp
#   src/Util.cpp
#   src/DGConfig.cpp
#   src/BoundaryConditions.cpp
#   src/PoroelasticConfig.cpp
#   src/ProbeConfig.cpp
#   src/ThermalConfig.cpp
#   src/Thermal.cpp
#   src/FECalc.cpp
#   src/MeshProjection.cpp
#   src/ElemProjection.cpp
#   src/PENumericalConfig.cpp
#   src/FractureBoundary.cpp
#   src/DamageConfig.cpp
#   src/InjectionPoint.cpp
#   src/Report.cpp
#   src/PEDebugSys.cpp
#   src/WellMirrors.cpp
#   src/Volumetric.cpp
#   src/MeshFracture.cpp
)

# ADD_LIBRARY(DEMO_EQUIL_MEC
#   demos/equil_mecanico/elasticity_system.cpp
# )

# Indicar bibliotecas compartilhadas no sistema
# target_link_libraries(Modules ${GMSH_LIB})
target_link_libraries(Modules ${BOOST_FS})
target_link_libraries(Modules ${LIBMESH_LIBS})

# To enable stack traces on fail using boost
target_link_libraries(Modules "-lboost_stacktrace_backtrace -ldl -lbacktrace" )


# Compilar e lincar executável final
add_executable( harpy src/main/harpy.cpp  )
target_link_libraries( harpy Modules )

set(RT "${PROJECT_SOURCE_DIR}") ## proj root

# Tests
add_executable( ut-timestep demos/UNIT-TEST/01-TIMESTEP/timestep-main.cpp  )
target_link_libraries( ut-timestep Modules )

# add_executable( ut-bcconfig demos/UNIT-TEST/02-BCCONFIG/bcconfig-main.cpp  )
# target_link_libraries( ut-bcconfig Modules )

add_executable( ut-param-interp demos/UNIT-TEST/04-PARAM-INTERPOLATOR/main.cpp  )
target_link_libraries( ut-param-interp Modules )

add_executable( ut-model-config demos/UNIT-TEST/05-MODEL-CONFIG//main.cpp  )
target_link_libraries( ut-model-config Modules )

# if(EXISTS ${RT}/demos/UnitTests/FractureNodeSplit)
#   add_executable( FractureNodeSplit demos/UnitTests/FractureNodeSplit/main.cpp  )
#   target_link_libraries( FractureNodeSplit Modules )
# endif()

# DEBUG ## se precisar ver as variaveis
# get_cmake_property(_variableNames VARIABLES)
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
#     message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()

message(-----------------)
