#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([pomerol],[0.8.0],[pomerol@googlecode.com])
AC_COPYRIGHT(GNU General Public License v3)
#AC_CONFIG_SRCDIR([src/main.cpp])

PKG_CHECK_MODULES([eigen3], [eigen3 >= 3.0.0])
PKG_CHECK_MODULES([hdf5], [hdf5 >= 1.8.0])

AM_INIT_AUTOMAKE
AC_PROG_RANLIB

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.

# Checks for header files.
#AC_CHECK_HEADERS([limits.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
#AC_C_INLINE
#AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([strchr strpbrk])

AC_PROG_MAKE_SET

AC_CONFIG_FILES([Makefile
                 jsoncpp/Makefile
                 src/Makefile
                 jsoncpp/src/Makefile
               ])
                 #jsoncpp/example/Makefile
AC_OUTPUT
