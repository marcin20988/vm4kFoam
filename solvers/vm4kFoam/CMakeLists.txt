add_executable(vm4kFoam vm4kFoam.C)

# twoPhaseEuler libraries
include_directories ("$ENV{FOAM_SOLVERS}/multiphase/twoPhaseEulerFoam/twoPhaseSystem/lnInclude")
include_directories ("$ENV{FOAM_SOLVERS}/multiphase/twoPhaseEulerFoam/interfacialModels/lnInclude")
include_directories ("$ENV{FOAM_SOLVERS}/multiphase/twoPhaseEulerFoam/phaseIncompressibleTurbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SOLVERS}/multiphase/twoPhaseEulerFoam/")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/phaseIncompressible/PhaseIncompressibleTurbulenceModel")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/phaseCompressible/PhaseCompressibleTurbulenceModel")
# OF libraries
include_directories ("$ENV{FOAM_SRC}/finiteVolume/lnInclude")
include_directories ("$ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude")
include_directories ("$ENV{FOAM_SRC}/transportModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/incompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/phaseIncompressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/TurbulenceModels/compressible/lnInclude")
include_directories ("$ENV{FOAM_SRC}/fvOptions/lnInclude")
include_directories ("$ENV{FOAM_SRC}/meshTools/lnInclude")
include_directories ("$ENV{FOAM_SRC}/sampling/lnInclude")

include_directories ("${PROJECT_SOURCE_DIR}/src/kineticModel/")




target_link_openfoam_libraries(vm4kFoam
    fluidThermophysicalModels
    specie
    turbulenceModels
    incompressibleTurbulenceModels
    incompressibleTransportModels
    compressibleEulerianInterfacialModels
    meshTools
    finiteVolume
    compressibleTwoPhaseSystem
    phaseCompressibleTurbulenceModels
    compressibleTurbulenceModels
    fvOptions
    OpenFOAM)

  target_link_libraries(vm4kFoam vm4kModel)

set_target_properties(vm4kFoam PROPERTIES RUNTIME_OUTPUT_DIRECTORY "$ENV{FOAM_USER_APPBIN}")
