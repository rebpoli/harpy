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
set(METHOD $ENV{METHOD})

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
set(CMAKE_CXX_COMPILER /usr/bin/mpicxx)

# execute_process(COMMAND ${LIBMESH_CONFIG} --cc  OUTPUT_VARIABLE CMAKE_C_COMPILER OUTPUT_STRIP_TRAILING_WHITESPACE)
# execute_process(COMMAND ${LIBMESH_CONFIG} --cxx OUTPUT_VARIABLE CMAKE_CXX_COMPILER OUTPUT_STRIP_TRAILING_WHITESPACE)


# gmsh

find_library(BOOST boost_filesystem)
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

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# libmesh compiler flags
set( CMAKE_CXX_FLAGS " -D BOOST_STACKTRACE_LINK ${CMAKE_CXX_FLAGS} ${LIBMESH_CPPFLAGS} ${LIBMESH_CXXFLAGS} ${LIBMESH_INCLUDE}" )

message(CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS})
message(LIBMESH_VERSION: ${LIBMESH_VERSION})
message(LIBMESH_CPPFLAGS: ${LIBMESH_CPPFLAGS})
message(LIBMESH_CXXFLAGS: ${LIBMESH_CXXFLAGS})
message(LIBMESH_LIBS: ${LIBMESH_LIBS})

# Setup Autodiff
find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
    message(STATUS "Eigen library has been found.")
    add_definitions(-DAUTODIFF_EIGEN_FOUND)
else()
    message(STATUS "Eigen was not found.")
endif()

# Headers
include_directories(include)

add_compile_options(-Wno-parentheses)   # Avoid stupid warnings

# Compilar arquivos-fonte dentro da biblioteca estática libModules.a
ADD_LIBRARY(Modules
  src/config/BCConfig.cpp
  src/config/Config.cpp
  src/config/ConfigValidate.cpp
  src/config/InoutConfig.cpp
  src/config/MaterialConfig.cpp
  src/config/ModelConfig.cpp
  src/config/SolverConfig.cpp
  src/config/reader/InoutReader.cpp
  src/config/reader/MaterialReader.cpp
  src/config/reader/ModelReader.cpp
  src/config/reader/SolverReader.cpp
  src/config/reader/SpatialDataReader.cpp
  src/harpy/Global.cpp
  src/harpy/HarpyInit.cpp
  src/harpy/Log.cpp
  src/harpy/PetscLogger.cpp
  src/postproc/probes/Probe.cpp
  src/postproc/report/VPReport.cpp
  src/postproc/stress/StressPostProc.cpp
  src/postproc/stress/TensorInvariants.cpp
  src/restart/File.cpp
  src/restart/Header.cpp
  src/restart/Util.cpp
  src/solver/common/ExplicitMaterial.cpp
  src/solver/common/Solver.cpp
  src/solver/pressure/PSolverFile.cpp
  src/solver/thermal/TIfc.cpp
  src/solver/thermal/TSolverCnt.cpp
  src/solver/thermal/TSolverFile.cpp
  src/solver/viscoplastic/VPBC.cpp
  src/solver/viscoplastic/VPIfc.cpp
  src/solver/viscoplastic/VPMaterial.cpp
  src/solver/viscoplastic/VPProbe.cpp
  src/solver/viscoplastic/VPProps.cpp
  src/solver/viscoplastic/VPSolver.cpp
  src/solverloop/SLViscoplastic.cpp
  src/timeloop/Timeloop.cpp
  src/timeloop/Timestep.cpp
  src/util/DirManager.cpp
  src/util/GridFile.cpp
  src/util/Json.cpp
  src/util/MeshUtils.cpp
  src/util/Messages.cpp
  src/util/OutputOperators.cpp
  src/util/Stopwatch.cpp
  src/util/TeeBuff.cpp
)

# Common: export symbols so backtrace resolves function names
# (-rdynamic is key on Linux)
add_link_options(-Wl,--no-as-needed)

if(METHOD STREQUAL "devel")
  message(STATUS "METHOD=devel: enabling debug symbols and nice stack traces")
  add_compile_options(-O0 -g3 -ggdb3 -fno-omit-frame-pointer)
  add_link_options(-rdynamic)
elseif(METHOD STREQUAL "asan")
  message(STATUS "METHOD=asan: AddressSanitizer + symbols")
  add_compile_options(-O1 -g3 -ggdb3 -fno-omit-frame-pointer -fsanitize=address)
  add_link_options(-rdynamic -fsanitize=address)
else() # release
  message(STATUS "METHOD=release: optimized build (no symbol export)")
  add_compile_options(-O3)
  # If you *do* strip in release, do it explicitly; devel will not strip.
  # add_link_options(-s)
endif()



# Indicar bibliotecas compartilhadas no sistema
# target_link_libraries(Modules ${GMSH_LIB})
target_link_libraries(Modules ${BOOST_FS})
target_link_libraries(Modules ${LIBMESH_LIBS})

# To enable stack traces on fail using boost
target_link_libraries(Modules "-lboost_stacktrace_backtrace -ldl -lbacktrace -lboost_serialization -lboost_mpi -lboost_filesystem" )

# Compilar e lincar executável final
add_executable( harpy src/main/harpy.cpp  )
target_link_libraries( harpy Modules )

add_executable( msh2exo src/main/msh2exo.cpp  )
target_link_libraries( msh2exo Modules )

set(RT "${PROJECT_SOURCE_DIR}") ## proj root

# Tests
add_executable( ut-timestep demos/UNIT-TEST/01-TIMESTEP/timestep-main.cpp  )
target_link_libraries( ut-timestep Modules )

add_executable( ut-autodiff-jac demos/UNIT-TEST/07-JACOBIAN-AUTODIFF/autodiff-jac.cpp )
target_link_libraries( ut-autodiff-jac Modules )

add_executable( ut-bcconfig demos/UNIT-TEST/02-BCCONFIG/bcconfig-main.cpp  )
target_link_libraries( ut-bcconfig Modules )

add_executable( ut-material-setup demos/UNIT-TEST/03-MATERIAL-SETUP/main.cpp  )
target_link_libraries( ut-material-setup Modules )

add_executable( ut-param-interp demos/UNIT-TEST/04-PARAM-INTERPOLATOR/main.cpp  )
target_link_libraries( ut-param-interp Modules )

add_executable( ut-model-config demos/UNIT-TEST/05-MODEL-CONFIG//main.cpp  )
target_link_libraries( ut-model-config Modules )

add_executable( ut-viscoplastic demos/UNIT-TEST/06-VISCOPLASTIC/main.cpp  )
target_link_libraries( ut-viscoplastic Modules )

add_executable( map-mpi-localize demos/UNIT-TEST/08-MAP-MPI-LOCALIZE/map-mpi-localize.cpp  )
target_link_libraries( map-mpi-localize Modules )

add_executable( templated-namespaces demos/UNIT-TEST/09-TEMPLATED-NAMESPACES/templated-namespaces.cpp )
target_link_libraries( templated-namespaces Modules )

add_executable( ut-radial-temperature demos/UNIT-TEST/10-RADIAL-TEMPERATURE/main.cpp  )
target_link_libraries( ut-radial-temperature Modules )

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

