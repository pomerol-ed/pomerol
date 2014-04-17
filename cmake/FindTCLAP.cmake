#  Copyright Olivier Parcollet 2010, Andrey Antipov 2013.
#
# This module looks for tclap command option parser.
# It sets up : TCLAP_INCLUDE_DIR
# 

find_package(PkgConfig)
pkg_check_modules(PC_TCLAP QUIET tclap)

SET(TRIAL_PATHS
 $ENV{TCLAP_ROOT}/include
 ${TCLAP_ROOT}/include
 /usr/include
 /usr/local/include
 /opt/local/include
 /sw/include
 )

find_path(TCLAP_INCLUDE_DIR tclap/CmdLine.h
          HINTS ${PC_TCLAP_INCLUDEDIR} ${TRIAL_PATHS}
          DOC "Include for TCLAP"
         )

set(TCLAP_INCLUDE_DIRS ${TCLAP_INCLUDE_DIR} )

mark_as_advanced(TCLAP_INCLUDE_DIR)

find_package_handle_standard_args(TCLAP DEFAULT_MSG TCLAP_INCLUDE_DIR)

