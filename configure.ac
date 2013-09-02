#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([fusionseq], [0.7.2], [david.z.chen@outlook.com])
AC_CONFIG_SRCDIR([src/geneFusions.c])
AC_CONFIG_SRCDIR([src/cgi/geneFusions_cgi.c])
AC_CONFIG_SRCDIR([src/test/bp2fasta.c])
AC_CONFIG_SRCDIR([src/optional/bp2wig.c])
AC_CONFIG_SRCDIR([scripts/Confp.pm])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

#------------------------------------------------------------------------------
# Initialize automake
#------------------------------------------------------------------------------
AM_INIT_AUTOMAKE([-Wall -Werror -Wno-portability subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

#------------------------------------------------------------------------------
# Initialize libtool
#------------------------------------------------------------------------------
LT_INIT
AC_SUBST([FUSIONSEQ_SO_VERSION], [2:7:0])

#------------------------------------------------------------------------------
# Options for disabling features.
#------------------------------------------------------------------------------

AC_ARG_ENABLE([cgi], 
    AS_HELP_STRING([--disable-cgi], [Disabling building CGIs]))

AM_CONDITIONAL([BUILD_CGI], [test "x$enable_cgi" != "xno"])

AC_ARG_ENABLE([tests], 
    AS_HELP_STRING([--disable-tests], [Disabling building test programs.]))

AM_CONDITIONAL([BUILD_TESTS], [test "x$enable_tests" != "xno"])

AC_ARG_ENABLE([optional], 
    AS_HELP_STRING([--disable-optional], [Disabling building optional programs.]))

AM_CONDITIONAL([BUILD_OPTIONAL], [test "x$enable_optional" != "xno"])

#------------------------------------------------------------------------------
# Checks for programs.
#------------------------------------------------------------------------------
AC_PROG_SED
AC_PROG_MKDIR_P
AC_PROG_CC
AC_PROG_CXX
ROOT_PATH([], [], [AC_MSG_ERROR([Cannot find [ROOT]])])
AC_PROG_LIBTOOL
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN"; then
	AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi
AM_CONDITIONAL([USE_DOXYGEN], [test -n "$DOXYGEN"])

#------------------------------------------------------------------------------
# Checks for libraries.
#------------------------------------------------------------------------------
AC_CHECK_LIB([m], [log], [], [AC_MSG_ERROR([Cannot find standard math library])])
AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [AC_MSG_ERROR([Cannot find cblas library])])
AC_CHECK_LIB([gsl], [gsl_ran_hypergeometric_pdf], [], [AC_MSG_ERROR([Cannot find gsl library])])
AC_CHECK_LIB([bios], [needMem], [], [AC_MSG_ERROR([Cannot find bios library])])
AC_CHECK_LIB([mrf], [mrf_init], [], [AC_MSG_ERROR([Cannot find mrf library])])
AX_CHECK_GD

#------------------------------------------------------------------------------
# Checks for header files.
#------------------------------------------------------------------------------
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])
AC_CHECK_HEADERS([gd.h])

#------------------------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
#------------------------------------------------------------------------------
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

#------------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------------
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([getcwd pathconf rint setenv strchr strdup strstr])

#------------------------------------------------------------------------------
# Output
#------------------------------------------------------------------------------
AM_COND_IF([USE_DOXYGEN], [AC_CONFIG_FILES([doc/Doxyfile])])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
