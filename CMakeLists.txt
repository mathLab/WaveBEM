##
# CMake script for large project, with one executable
# and a library of objects depending on deal2lkit.
#
# The structure of the directory is assumed to be: 
# ./source: containing implementations and one file with "int main()"
# ./include: containing all class declarations
# ./tests: containing pairs of files: test.cc/test.output to test
#
# If you set the environemnt variable DEAL2LKIT_DIR or D2K_DIR, 
# everything will work out of the box
##

# If your application follows the structure above, you don't need to 
# specify anything else. 
SET(TARGET WaveBEM)
SET(_main source/main.cc)

# Declare all source files the target consists of:
file(GLOB _files source/*cc)

# Don't compile the main file into the library
IF(NOT "${_main}" STREQUAL "") 
  LIST(REMOVE_ITEM _files
    ${CMAKE_CURRENT_SOURCE_DIR}/${_main}
    )
ENDIF()


CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.3 REQUIRED
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES()

INCLUDE_DIRECTORIES(include ${DEAL_II_DIR}/include/deal.II)


PROJECT(${TARGET})
FIND_PACKAGE(deal2lkit 1.0 REQUIRED
  HINTS ${D2K_DIR} $ENV{D2K_DIR} $ENV{DEAL2LKIT_DIR}
  )
D2K_INITIALIZE_CACHED_VARIABLES()
OPTION(ENABLE_APP_TESTING "Enable tests in this application." ON)
OPTION(ENABLE_DOCUMENTATION "Enable Doxygen documentation." ON)

# We add one library and one target for each type of deal.II library
# we found. If you compiled deal.II with both Release and Debug
# mode, this will generate both Release and Debug programs for you.
# The debug library and program are postfixed with ".g"
SET(_d2_build_types "Release;Debug")
SET(Release_postfix "")
SET(Debug_postfix ".g")

FOREACH(_build_type ${_d2_build_types})
  # Postfix to use everywhere
  SET(_p "${${_build_type}_postfix}")
  # Only build this type, if deal.II was compiled with it.
  IF(CMAKE_BUILD_TYPE MATCHES "${_build_type}" AND
      DEAL_II_BUILD_TYPE MATCHES "${_build_type}"  AND
      D2K_BUILD_TYPE MATCHES "${_build_type}")

    MESSAGE("-- Found ${_build_type} version of deal.II.")
    MESSAGE("-- Found ${_build_type} version of deal2lkit.")

    SET(_lib "${TARGET}-lib${_p}")
    MESSAGE("-- Configuring library ${_lib}")

    STRING(TOUPPER "${_build_type}" _BUILD_TYPE)

    ADD_LIBRARY(${_lib} SHARED ${_files})
    SET(_exe "${TARGET}${${_build_type}_postfix}")
    MESSAGE("-- Configuring executable ${_exe}")
    ADD_EXECUTABLE(${_exe} ${_main})
    TARGET_LINK_LIBRARIES(${_exe} ${_lib})
    D2K_SETUP_TARGET(${_exe} ${_BUILD_TYPE})

    SET(TEST_LIBRARIES_${_BUILD_TYPE} ${_lib})

    D2K_SETUP_TARGET(${_lib} ${_BUILD_TYPE})
    SET_TARGET_PROPERTIES(${_lib} PROPERTIES MACOSX_RPATH OFF
	BUILD_WITH_INSTALL_RPATH OFF
	INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib
    )
  ENDIF()
ENDFOREACH()


if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests AND ${ENABLE_APP_TESTING})
  ADD_SUBDIRECTORY(tests)
  ENABLE_TESTING()
endif()


IF(ENABLE_DOCUMENTATION)
  ADD_SUBDIRECTORY(doxygen)
ENDIF()


ADD_CUSTOM_TARGET(indent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./scripts/indent
)

