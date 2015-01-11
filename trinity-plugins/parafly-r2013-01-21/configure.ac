AC_INIT([Parafly], [0.1], [bug-report@address])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_LANG([C++])
AC_PROG_CXX
#AC_OPENMP # requires autoconf >= 2.62
AC_SUBST([AM_CXXFLAGS], [-m64])
case $CXX in
  g++*) AC_SUBST([AM_CXXFLAGS],["-pedantic -fopenmp -Wall -Wextra -Wno-long-long -Wno-deprecated $AM_CXXFLAGS"]);;
  sunCC*) AC_SUBST([AM_CXXFLAGS], ["-library=stlport4 -xopenmp -xvpara -fast $AM_CXXFLAGS"]) ;;
  icpc*) AC_SUBST([AM_CXXFLAGS], ["-Wall -openmp $AM_CXXFLAGS"]) ;;
esac
AC_SEARCH_LIBS([cos], [m])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
