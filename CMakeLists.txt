cmake_minimum_required (VERSION 2.6)
project (vm4kSolvers)

#The build is based on OpenFOAM environmental varialbles Need to check 
#if they are loaded and if not TODO search and load an appropriate file?
if(DEFINED ENV{WM_PROJECT})
  message("\nOpenFOAM $ENV{WM_PROJECT_VERSION} found in $ENV{WM_PROJECT_DIR}.\n")
else()
  message("WM_PROJECT is not defined. You must set OpenFOAM specific environment
  variables for CMake to configure the build.")
  # exit early 
  return()
endif()

#TODO Is there a better way to break lines?
set(CMAKE_CXX_FLAGS "-m64 -Dlinux64 -DWM_DP -ldl ") 
#g++ will inform you of unused typedefs and we want to suppress this as there
#are typedefs like this in OF

if (CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNoRepository -ftemplate-depth-100")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Xlinker --add-needed -Xlinker --no-as-needed")
else ()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_HAS_TR1_TUPLE=0 -wd327,654,819,1125,1476,1505,1572 -xSSE3 -O2 -no-prec-div  -DNoRepository")
endif()

#define a macro to simplify adding OpenFOAM libraries
#TODO  add for each
macro (TARGET_LINK_OPENFOAM_LIBRARIES target)
  foreach(arg ${ARGN})
    if(NOT DEFINED ${arg}_LIB)
      find_library(${arg}_LIB lib${arg}.so $ENV{FOAM_LIBBIN})
    endif()
    target_link_libraries(${target} ${${arg}_LIB})
  endforeach(arg ${ARGN})
endmacro (TARGET_LINK_OPENFOAM_LIBRARIES)

#TODO  add for each
macro (TARGET_LINK_OPENFOAM_USERLIB arg1 arg2)
  if(NOT DEFINED ${arg2}_LIB)
    find_library(${arg2}_LIB lib${arg2}.so $ENV{FOAM_USER_LIBBIN})
  endif()
  #message("Using ${${arg2}_LIB}")
  target_link_libraries(${arg1} ${${arg2}_LIB})
endmacro (TARGET_LINK_OPENFOAM_USERLIB)


macro (ADD_MY_TEST testname executable)
  add_test(
      NAME ${testname}
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${testname}"
      COMMAND ${executable}
  )

  add_custom_command(
      TARGET ${executable}
      PRE_BUILD
      COMMAND mkdir -p ${testname}
      COMMAND cp -r ${PROJECT_SOURCE_DIR}/testing/${testname}/* ${testname}
      COMMENT "Creating and copying directories for ${testname}"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_custom_command(
      TARGET ${executable}
      POST_BUILD
      COMMAND blockMesh &> log.blockMesh
      COMMENT "Running blockMesh in ${CMAKE_CURRENT_BINARY_DIR}/${testname}"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${testname}
  )
endmacro (ADD_MY_TEST)

#macro to list all subdirectories
MACRO(SUBDIRLIST result curdir)
FILE(GLOB children ${curdir} ${curdir}/*)
SET(dirlist "")
FOREACH(child ${children})
	IF(IS_DIRECTORY ${child})
		SET(dirlist ${dirlist} ${child})
	ENDIF()
ENDFOREACH()
SET(${result} ${dirlist})
ENDMACRO()


#option (UseOldLibs "Use old libraries" OFF)
#General includes required by all code
include_directories ("$ENV{FOAM_SRC}/OpenFOAM/lnInclude")
include_directories ("$ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude")

#Less general but seems to be used by everything too?
include_directories ("$ENV{FOAM_SRC}/transportModels")
include_directories ("$ENV{FOAM_SRC}/finiteVolume/lnInclude")

#We need to tell CMake to use use C files as CXX (OpenFOAM standard)
set_source_files_properties(filename.c PROPERTIES LANGUAGE CXX )

#Library
add_subdirectory(src/kineticModel)

#Solvers
add_subdirectory(solvers/vm4kFoam)

