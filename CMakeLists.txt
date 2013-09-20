CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

# Include all macros. Some of them are copied brutally from deal.II
FILE(GLOB _cmake_files "cmake/macros/*.cmake")
FOREACH(_file ${_cmake_files})
  INCLUDE(${_file})
ENDFOREACH()

# These ones are for GSL, OCE and SUNDIALS
LIST(APPEND CMAKE_MODULE_PATH
  $ENV{MATHLAB_SRC}/cmake/
  $ENV{MATHLAB_SRC}/cmake/modules/
  )



FIND_PACKAGE(deal.II 8.0 REQUIRED HINTS $ENV{DEAL_II_DIR})
FIND_PACKAGE(OpenCASCADE MODULE REQUIRED)


SET(TARGET "waveBem")
set(CMAKE_BUILD_TYPE "Release")

file(GLOB TARGET_SRC source/*cc)
include_directories(
	include 
	$ENV{MATHLAB_SRC}/include 
	${DEAL_II_INCLUDE_DIRS}
	${OpenCASCADE_INCLUDE_DIRS}
	$ENV{SUNDIALS_DIR}/include
	$ENV{GSL_DIR}/include
)

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()

STRING(TOLOWER ${CMAKE_BUILD_TYPE} _build)
FIND_LIBRARY(mathlab mathlab.${_build} 
	     HINTS $ENV{MATHLAB_SRC}/lib)
MESSAGE('mathlab.${_build}')
MESSAGE(${mathlab})
target_link_libraries(
	${TARGET} 
	${OpenCASCADE_LIBRARIES}
	${mathlab}
	gsl gslcblas
	-L$ENV{GSL_DIR}/lib
)
