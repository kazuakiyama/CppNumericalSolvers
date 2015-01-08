# - Try to find Eigen3
# Once done this will define
#  EIGEN3_FOUND - System has Eigen3
#  EIGEN3_INCLUDE_DIRS - The Eigen3 include directories
#  EIGEN3_LIBRARY_DIRS - The library directories needed to use Eigen3
#  EIGEN3_LIBRARIES    - The libraries needed to use Eigen3

find_package(PkgConfig QUIET)

if (PKG_CONFIG_FOUND)
  if (EIGEN3_FIND_VERSION)
    if (EIGEN3_FIND_VERSION_EXACT)
      pkg_check_modules(_PC_EIGEN3 QUIET eigen3=${EIGEN3_FIND_VERSION})
    else (EIGEN3_FIND_VERSION_EXACT)
      pkg_check_modules(_PC_EIGEN3 QUIET eigen3>=${EIGEN3_FIND_VERSION})
    endif (EIGEN3_FIND_VERSION_EXACT)
  else (EIGEN3_FIND_VERSION)
    pkg_check_modules(_PC_EIGEN3 QUIET eigen3)
  endif (EIGEN3_FIND_VERSION)

  if (_PC_EIGEN3_FOUND)
    set (EIGEN3_INCLUDE_DIRS ${_PC_EIGEN3_INCLUDE_DIRS} CACHE PATH "EIGEN3 include directory")
    set (EIGEN3_DEFINITIONS ${_PC_EIGEN3_CFLAGS_OTHER} CACHE STRING "Additional compiler flags for EIGEN3")
  endif (_PC_EIGEN3_FOUND)

  include (FindPackageHandleStandardArgs)

  # handle the QUIETLY and REQUIRED arguments and set LIBCPLEX_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args (EIGEN3 DEFAULT_MSG  EIGEN3_INCLUDE_DIRS)

  mark_as_advanced(EIGEN3_INCLUDE_DIRS EIGEN3_DEFINITIONS)

endif (PKG_CONFIG_FOUND)
