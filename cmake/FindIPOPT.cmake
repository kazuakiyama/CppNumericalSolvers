# - Try to find IPOPT
# Once done this will define
#  IPOPT_FOUND - System has IpOpt
#  IPOPT_INCLUDE_DIRS - The IpOpt include directories
#  IPOPT_LIBRARY_DIRS - The library directories needed to use IpOpt
#  IPOPT_LIBRARIES    - The libraries needed to use IpOpt

if (IPOPT_INCLUDE_DIR)
  # in cache already
  set (IPOPT_FIND_QUIETLY TRUE)
endif (IPOPT_INCLUDE_DIR)

find_package(PkgConfig QUIET)

if (PKG_CONFIG_FOUND)
  if (IPOPT_DIR)
    set (_pkgconfig_path_old $ENV{PKG_CONFIG_PATH})
    set (ENV{PKG_CONFIG_PATH} "${_pkgconfig_path_old}:${IPOPT_DIR}/lib/pkgconfig")
  endif (IPOPT_DIR)
  if (IPOPT_FIND_VERSION)
    if (IPOPT_FIND_VERSION_EXACT)
      pkg_check_modules(_PC_IPOPT QUIET ipopt=${IPOPT_FIND_VERSION})
    else (IPOPT_FIND_VERSION_EXACT)
      pkg_check_modules(_PC_IPOPT QUIET ipopt>=${IPOPT_FIND_VERSION})
    endif (IPOPT_FIND_VERSION_EXACT)
  else (IPOPT_FIND_VERSION)
    pkg_check_modules(_PC_IPOPT QUIET ipopt)
  endif (IPOPT_FIND_VERSION)
  set (ENV{PKG_CONFIG_PATH} ${_pkgconfig_path_old})

  if (_PC_IPOPT_FOUND)
    set (IPOPT_INCLUDE_DIRS ${_PC_IPOPT_INCLUDE_DIRS} CACHE PATH "IPOPT include directory")
    set (IPOPT_DEFINITIONS ${_PC_IPOPT_CFLAGS_OTHER} CACHE STRING "Additional compiler flags for IPOPT")
    set (IPOPT_LIBRARIES "" CACHE STRING "IPOPT libraries" FORCE)
    foreach(_LIBRARY IN ITEMS ${_PC_IPOPT_LIBRARIES})
      find_library(${_LIBRARY}_PATH
        NAMES ${_LIBRARY}
        PATHS ${_PC_IPOPT_LIBRARY_DIRS})
      list (APPEND IPOPT_LIBRARIES ${${_LIBRARY}_PATH})
    endforeach(_LIBRARY)
  endif (_PC_IPOPT_FOUND)

  include (FindPackageHandleStandardArgs)

  # handle the QUIETLY and REQUIRED arguments and set LIBCPLEX_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args (IPOPT DEFAULT_MSG IPOPT_LIBRARIES IPOPT_INCLUDE_DIRS)
  mark_as_advanced(IPOPT_INCLUDE_DIRS IPOPT_LIBRARIES IPOPT_DEFINITIONS)

else (PKG_CONFIG_FOUND)

  if (WIN32)
    find_path(IPOPT_INCLUDE_DIR NAMES IpNLP.hpp
      PATHS
      "C:\\libs\\Ipopt-3.8.2\\include\\coin"
      ${IPOPT_DIR}/include
      )

    if (IPOPT_INCLUDE_DIR)
      find_library( IPOPT_LIBRARY_RELEASE 
        Ipopt
        PATHS "C:\\libs\\Ipopt-3.8.2\\lib\\win32\\release" )
      find_library( IPOPT_LIBRARY_DEBUG
        Ipopt
        PATHS "C:\\libs\\Ipopt-3.8.2\\lib\\win32\\debug" )

      set (IPOPT_LIBRARY "optimized;${IPOPT_LIBRARY_RELEASE};debug;${IPOPT_LIBRARY_DEBUG}" 
        CACHE  STRING "IPOPT Libraries" )

      set (IPOPT_FOUND TRUE)
      set (IPOPT_INCLUDE_DIR ${IPOPT_INCLUDE_DIR})
      # Todo, set right version depending on build type (debug/release)
      #GET_FILENAME_COMPONENT( IPOPT_LIBRARY_DIR ${GLEW_LIBRARY} PATH )
    else (IPOPT_INCLUDE_DIR)
      set (IPOPT_FOUND FALSE)
      set (IPOPT_INCLUDE_DIR ${IPOPT_INCLUDE_DIR})
    endif (IPOPT_INCLUDE_DIR)

  else ( WIN32 )
    find_path(IPOPT_INCLUDE_DIR NAMES IpNLP.hpp
      PATHS  "$ENV{IPOPT_HOME}/include/coin"
      "/usr/include/coin"
      "${IPOPT_DIR}/include/coin"
      )

    find_library (IPOPT_LIBRARY 
      ipopt
      PATHS "$ENV{IPOPT_HOME}/lib"
      "/usr/lib"
      "${IPOPT_DIR}/lib"
      )
    
    #wrong config under Debian workaround
    #add_definitions( -DHAVE_CSTDDEF )
    
    # set optional path to HSL Solver
    find_path(IPOPT_HSL_LIBRARY_DIR 
      NAMES libhsl.so
      libhsl.dylib
      PATHS "$ENV{IPOPT_HSL_LIBRARY_PATH}"
      "$ENV{HOME}/opt/HSL/lib"
      )
    
    if ( IPOPT_HSL_LIBRARY_DIR)
      if ( NOT IPOPT_FIND_QUIETLY )
        message ( "IPOPT_HSL_LIBRARY_DIR found at ${IPOPT_HSL_LIBRARY_DIR} ")
      endif ()
      set (IPOPT_LIBRARY_DIR ${IPOPT_HSL_LIBRARY_DIR})
      list ( APPEND IPOPT_LIBRARY_DIRS "${IPOPT_HSL_LIBRARY_DIR}")
    endif (IPOPT_HSL_LIBRARY_DIR)
    
    
    set (IPOPT_INCLUDE_DIRS "${IPOPT_INCLUDE_DIR}" )
    set (IPOPT_LIBRARIES "${IPOPT_LIBRARY}" )

    include (FindPackageHandleStandardArgs)

    # handle the QUIETLY and REQUIRED arguments and set LIBCPLEX_FOUND to TRUE
    # if all listed variables are TRUE
    find_package_handle_standard_args (IPOPT DEFAULT_MSG IPOPT_LIBRARIES IPOPT_INCLUDE_DIRS)

    mark_as_advanced(IPOPT_INCLUDE_DIRS IPOPT_LIBRARIES )
    
  endif (WIN32)

endif (PKG_CONFIG_FOUND)
