AC_INIT([dash],[0.5.10.2])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_MACRO_DIRS([m4])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES(yes)])

dnl Checks for programs.
AC_PROG_CC
AC_USE_SYSTEM_EXTENSIONS
AC_PROG_YACC
dnl MMG 2018-09-26 support building the library
AM_PROG_AR


AC_MSG_CHECKING([for build system compiler])
if test "$cross_compiling" = yes; then
	CC_FOR_BUILD=${CC_FOR_BUILD-cc}
else
	CC_FOR_BUILD=${CC}
fi
AC_MSG_RESULT(${CC_FOR_BUILD})
AC_SUBST(CC_FOR_BUILD)

AC_MSG_CHECKING([for __attribute__((__alias__()))])
dash_cv_have_attribute_alias=no
AC_LINK_IFELSE([AC_LANG_PROGRAM([void t() {}
                                 void a() __attribute__((__alias__("t")));],
                                [a();])],
               [dash_cv_have_attribute_alias=yes])
AC_MSG_RESULT($dash_cv_have_attribute_alias)
if test "x$dash_cv_have_attribute_alias" = xyes; then
  AC_DEFINE([HAVE_ALIAS_ATTRIBUTE], 1,
            [Define if __attribute__((__alias__())) is supported])
fi

AC_ARG_ENABLE(static, AS_HELP_STRING(--enable-static, \
				     [Build statical linked program]))
if test "$enable_static" = "yes"; then
	export LDFLAGS="-static -Wl,--fatal-warnings"
fi

AC_ARG_ENABLE(fnmatch, AS_HELP_STRING(--enable-fnmatch, \
				      [Use fnmatch(3) from libc]))
AC_ARG_ENABLE(glob, AS_HELP_STRING(--enable-glob, [Use glob(3) from libc]))

dnl Checks for libraries.

dnl Checks for header files.
AC_CHECK_HEADERS(alloca.h paths.h)

dnl Check for declarations
AC_CHECK_DECL([_PATH_BSHELL],,AC_DEFINE_UNQUOTED([_PATH_BSHELL], "/bin/sh", [Define to system shell path]),[
#ifdef HAVE_PATHS_H
#include <paths.h>
#endif
])
AC_CHECK_DECL([_PATH_DEVNULL],,AC_DEFINE_UNQUOTED([_PATH_DEVNULL], "/dev/null", [Define to devnull device node path]),[
#ifdef HAVE_PATHS_H
#include <paths.h>
#endif
])
AC_CHECK_DECL([_PATH_TTY],,AC_DEFINE_UNQUOTED([_PATH_TTY], "/dev/tty", [Define to tty device node path]),[
#ifdef HAVE_PATHS_H
#include <paths.h>
#endif
])

dnl Some systems lack isblank
AC_CHECK_DECLS([isblank],,,[#include <ctype.h>])

dnl Check for sizes of types
AC_CHECK_SIZEOF([intmax_t])
AC_CHECK_SIZEOF([long long int])

dnl Select a fallback format string for intmax_t in case we don't find PRIdMAX
if test "x$ac_cv_sizeof_intmax_t" = "x$ac_cv_sizeof_long_long_int"; then
  intmax_fstr="lld"
else
  intmax_fstr="jd"
fi

dnl Check for PRIdMAX and define it to a fallback if not found
AC_CHECK_DECL([PRIdMAX],,
	[AC_DEFINE_UNQUOTED([PRIdMAX], "$intmax_fstr",
				       [Define to printf format string for intmax_t])],
        [
#include <inttypes.h>
])

dnl Checks for library functions.
AC_CHECK_FUNCS(bsearch faccessat getpwnam getrlimit isalpha killpg \
	       mempcpy \
	       sigsetmask stpcpy strchrnul strsignal strtod strtoimax \
	       strtoumax sysconf)

dnl Check whether it's worth working around FreeBSD PR kern/125009.
dnl The traditional behavior of access/faccessat is crazy, but
dnl POSIX.1-2008 explicitly allows those functions to misbehave.
dnl
dnl Unaffected kernels:
dnl
dnl - all versions of Linux
dnl - NetBSD sys/kern/vfs_subr.c 1.64, 1997-04-23
dnl - FreeBSD 9 (r212002), 2010-09-10
dnl - OpenBSD sys/kern/vfs_subr.c 1.166, 2008-06-09
dnl
dnl Also worked around in Debian's libc0.1 2.13-19 when using
dnl kFreeBSD 8.

AC_ARG_ENABLE(test-workaround, AS_HELP_STRING(--enable-test-workaround, \
	[Guard against faccessat(2) that tells root all files are executable]),,
	[enable_test_workaround=auto])

if test "enable_test_workaround" = "auto" &&
   test "$ac_cv_func_faccessat" = yes; then
	case `uname -s 2>/dev/null` in
	GNU/kFreeBSD | \
	FreeBSD)
		enable_test_workaround=yes
	esac
fi
if test "$enable_test_workaround" = "yes"; then
	AC_DEFINE([HAVE_TRADITIONAL_FACCESSAT], [1],
		[Define if your faccessat tells root all files are executable])
fi

if test "$enable_fnmatch" = yes; then
	use_fnmatch=
	AC_CHECK_FUNCS(fnmatch, use_fnmatch=yes)
fi

if test "$use_fnmatch" = yes && test "$enable_glob" = yes; then
	AC_CHECK_FUNCS(glob)
fi

dnl Check for klibc signal.
AC_CHECK_FUNC(signal)
if test "$ac_cv_func_signal" != yes; then
	AC_CHECK_FUNC(bsd_signal,
		      [AC_DEFINE(signal, bsd_signal,
				 [klibc has bsd_signal instead of signal])])
fi

dnl Check for stat64 (dietlibc/klibc).
AC_CHECK_FUNC(stat64,, [
	AC_DEFINE(fstat64, fstat, [64-bit operations are the same as 32-bit])
	AC_DEFINE(lstat64, lstat, [64-bit operations are the same as 32-bit])
	AC_DEFINE(stat64, stat, [64-bit operations are the same as 32-bit])
])

AC_CHECK_FUNC(open64,, [
	AC_DEFINE(open64, open, [64-bit operations are the same as 32-bit])
])

dnl Check if struct stat has st_mtim.
AC_MSG_CHECKING(for stat::st_mtim)
AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>],
[struct stat foo; return sizeof(foo.st_mtim.tv_sec)])],
have_st_mtim=yes, have_st_mtim=no)
AC_MSG_RESULT($have_st_mtim)
if test "$have_st_mtim" = "yes"; then
	AC_DEFINE([HAVE_ST_MTIM], [1],
		[Define if your `struct stat' has `st_mtim'])
fi

AC_ARG_WITH(libedit, AS_HELP_STRING(--with-libedit, [Compile with libedit support]))
use_libedit=
if test "$with_libedit" = "yes"; then
	AC_CHECK_LIB(edit, history_init, [
		AC_CHECK_HEADER([histedit.h], [use_libedit="yes"],
				AC_MSG_ERROR(
					[Can't find required header files.]))])
fi
if test "$use_libedit" != "yes"; then
	AC_DEFINE([SMALL], 1, [Define if you build with -DSMALL])
else
	export LIBS="$LIBS -ledit"
fi
AC_ARG_ENABLE(lineno, AS_HELP_STRING(--disable-lineno, \
				     [Disable LINENO support]))
if test "$enable_lineno" != "no"; then
	AC_DEFINE([WITH_LINENO], 1, [Define if you build with -DWITH_LINENO])
fi

dnl MMG 2018-09-27 support building the shared library
AC_PROG_LIBTOOL
LT_INIT

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
