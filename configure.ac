#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([pomerol],[0.8.0],[pomerol@googlecode.com])
AC_COPYRIGHT(GNU General Public License v3)
AC_CONFIG_SRCDIR([src/main.cpp])


# Checks for programs.
AC_PROG_CXX

AM_INIT_AUTOMAKE
AC_PROG_RANLIB
# Checks for libraries.
# start: Check for HDF5 support
AC_ARG_WITH([hdf5],
            [AS_HELP_STRING([--with-hdf5=path], [location of the hdf5 library])],[
                hdf5_LIBS="-L$withval -lhdf5 -lhdf5_cpp"
                hdf5_CFLAGS="-I$withval"
                AC_SUBST([hdf5_CFLAGS])
                AC_SUBST([hdf5_LIBS])
                ],
            []
           )

AC_ARG_WITH([hdf5-include],
            [AS_HELP_STRING([--with-hdf5-include=path], [location of the hdf5 headers, defaults to /usr/include])],[ 
                hdf5_CFLAGS="-I$withval"
                AC_SUBST([hdf5_CFLAGS])
                ],
            [])

AC_ARG_WITH([hdf5-lib],
            [AS_HELP_STRING([--with-hdf5-lib=path], [location of the hdf5 libraries])], [
                hdf5_LIBS="-L$withval -lhdf5 -lhdf5_cpp"
                AC_SUBST([hdf5_LIBS])
            ],
            [])


PKG_CHECK_MODULES([hdf5], [hdf5 >= 1.8.0],[
                    AC_SUBST([hdf5_CFLAGS])
                    AC_SUBST([hdf5_LIBS])
                    ],
                    [
                    hdf5_LIBS="-L/usr/lib -lhdf5 -lhdf5_cpp"
                    hdf5_CFLAGS="-I/usr/include"
                    AC_MSG_WARN(["No hdf5 >= 1.8.0 found automatically. Check your \$PKG_CONFIG_PATH or set --with-hdf5 flag. Assuming hdf5 is stored in /usr "])
                    ]
                  )
AC_CHECK_HEADERS([hdf5.h],,AC_MSG_ERROR(cannot find hdf5.h in hdf5 include path))
AC_CHECK_LIB(hdf5, H5Fflush,,AC_MSG_ERROR([invalid hdf5 in hdf5 path]))
AC_CHECK_LIB(hdf5_hl, H5DSis_scale,,AC_MSG_ERROR(invalid hdf5 high-level library in hdf5 path))

#End :: Check HDF5 suport
#Check Eigen3 via pkg-config 
PKG_CHECK_MODULES([eigen3], [eigen3 >= 3.0.0])
#End




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
