#--------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
#--------------------------------------------------------------------

AC_PREREQ([2.59])
AC_INIT([Chombo], [3.2], [chombousers-request@hpcrdm.lbl.gov], [Chombo.tar.gz], [https://seesar.lbl.gov/anag/chombo/index.html])
AC_CONFIG_SRCDIR([lib/src/BaseTools/REAL.H])

#--------------------------------------------------------------------
# Created files
#--------------------------------------------------------------------

AC_CONFIG_HEADERS([./lib/src/BaseTools/CH_config.H])
AC_CONFIG_FILES([./lib/mk/Make.defs.defaults_ac])

#--------------------------------------------------------------------
# System information
#--------------------------------------------------------------------

AC_CANONICAL_HOST
ch_int_systemptrsize="unknown"
AC_CHECK_SIZEOF([void *])
if test "$ac_cv_sizeof_void_p" -eq "8"; then
  ch_int_systemptrsize="64"
elif test "$ac_cv_sizeof_void_p" -eq "4"; then
  ch_int_systemptrsize="32"
fi

#--------------------------------------------------------------------
# Custom tests use C++
#--------------------------------------------------------------------

AC_PROG_CXX()
AC_LANG([C++])

# Try to determine the make of the compiler.  Many others than gcc set GXX to
# yes so we have to be a bit more clever.
AC_MSG_CHECKING([the make of the C++ compiler])
ch_str_cxxmake=unknown
case $host in
  *86*)
    if test "x`$CXX --version 2>&1 | sed '1d;2{s/.*Free Software Foundation.*/g++/;q;}'`" = "xg++"; then
      ch_str_cxxmake=gcc
    elif test "x`$CXX -V 2>&1 | sed '1{s/.*Intel.*/icpc/;q;}'`" = "xicpc"; then
      ch_str_cxxmake=intel
    elif test "x`$CXX -V 2>&1 | sed '1d;2d;3{s/.*Portland.*/pgCC/;q;}'`" = "xpgCC"; then
      ch_str_cxxmake=portland
    elif test "x$GXX" = "xyes"; then
      # Pathscale will usually come in here.
      ch_str_cxxmake=gcc
    fi
    AC_MSG_RESULT([$ch_str_cxxmake])
    if test "x$ch_str_cxxmake" = "xunknown"; then
      AC_MSG_WARN([Unknown compiler: testing of x86 hardware may fail simply because the right flags are not being used])
    fi
  ;;
  *    )
    AC_MSG_RESULT([$ch_str_cxxmake])
  ;;
esac

#--------------------------------------------------------------------
# Check for memalign
#--------------------------------------------------------------------

AC_MSG_CHECKING([how to allocate aligned memory])
ch_stat_memalign=unknown

# Check for posix memalign
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#define _XOPEN_SOURCE 600
#include <stdlib.h>]],
[[void* addr;
int m_alignment = 64;
int ier = posix_memalign(&addr, m_alignment, 64);]])],
                  [ch_stat_memalign=posix_memalign])

# Results
AC_MSG_RESULT([$ch_stat_memalign])
if test "x$ch_stat_memalign" = "xposix_memalign"; then
  AC_DEFINE([CHDEF_SYSTEM_HAVE_ALIGNEDALLOC])
  AC_DEFINE([CHDEF_SYSTEM_HAVE_POSIXMEMALIGN])
fi
if test "x$ch_stat_memalign" = "xunknown"; then
  AC_MSG_WARN([No method found for allocating aligned memory - this has
  possible performance implications])
fi

#--------------------------------------------------------------------
# Check for cache line size (unable if cross-compiling)
#--------------------------------------------------------------------

ch_stat_cacheinfo=unknown
ch_int_cachelinesize=unknown
ch_file_cacheinfo=unknown

# Preliminaries
AC_CHECK_FILE([/sys/devices/system/cpu/cpu0/cache/index0])
AC_MSG_CHECKING([cache line size])

# Try system files
if test "x$ac_cv_file__sys_devices_system_cpu_cpu0_cache_index0" = "xyes"; then
  ch_file_cacheinfo=/sys/devices/system/cpu/cpu0/cache
  for d in 0 1; do
    if test -r "$ch_file_cacheinfo/index$d/type"; then
      ch_str_cachetype=`cat $ch_file_cacheinfo/index$d/type`
    fi
    if test "x$ch_str_cachetype" = "xData"; then
      if test -r "$ch_file_cacheinfo/index$d/coherency_line_size"; then
        ch_int_cachelinesize=`cat $ch_file_cacheinfo/index$d/coherency_line_size`
        ch_stat_cacheinfo=file
        break
      fi
    fi
  done
fi

# Try sysconf
if test "x$ch_stat_cacheinfo" = "xunknown"; then
  ch_stat_sysconfcache=no
  AC_RUN_IFELSE([AC_LANG_PROGRAM(
[[#include <unistd.h>
#include <fstream>]],
[[std::ofstream ofile("conftest.data");
ofile << sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
ofile.close();]])],
                [ch_stat_sysconfcache=yes])
  if test "x$ch_stat_sysconfcache" = "xyes" && test -r "conftest.data"; then
    ch_int_cachelinesize=`cat conftest.data`
    ch_stat_cacheinfo=sysconf
  fi
fi

# Results
AC_MSG_RESULT([$ch_int_cachelinesize])
if test "x$ch_stat_cacheinfo" = "xunknown"; then
  AC_MSG_WARN([Unable to find cache line size - using default from
  CH_config.H.in])
else
  AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CACHELINESIZE], [$ch_int_cachelinesize])
  if test "x$ch_stat_cacheinfo" = "xfile"; then
    AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_HAVE_CACHEINFOFILE], [$ch_file_cacheinfo])
  elif test "x$ch_stat_cacheinfo" = "xsysconf"; then
    AC_DEFINE([CHDEF_SYSTEM_HAVE_CACHESYSCONF])
  fi
fi

#--------------------------------------------------------------------
# Check for x86 vector extensions and intrinsics
#--------------------------------------------------------------------

ch_str_veccpu_feature=unknown
ch_str_veccxx_feature=unknown
ch_int_veccpu_featurebits=0
ch_int_veccxx_featurebits=0
ch_str_sse_h=unknown
ch_str_sse2_h=unknown
ch_str_sse3_h=unknown
ch_str_ssse3_h=unknown
ch_str_sse4d1_h=unknown
ch_str_sse4d2_h=unknown
ch_str_sse4a_h=unknown
ch_str_avx_h=unknown
ch_str_fma_h=unknown
ch_str_avx2_h=unknown
ch_str_tsxhle_h=unknown
ch_str_tsxrtm_h=unknown
ch_str_avx512f_h=unknown
ch_str_avx512cd_h=unknown
ch_str_avx512dq_h=unknown
ch_str_avx512bw_h=unknown
ch_str_avx512fma_h=unknown
ch_str_avx512vbmi_h=unknown
ch_str_intelvec_h=unknown
# Vector operators are tested for sse, sse2, avx, and avx2
ch_stat_vecopp=unknown

AC_MSG_CHECKING([for x86 vector extensions])

case $host_cpu in
  *86*)

#--Test hardware

    ch_stat_cpuid=no
    # Avoid using an optimization flags -- the assembly can cause problems.
    # We will also gradually add sse flags
    ch_str_savecxxflags="$CXXFLAGS"
    CXXFLAGS=""
    AC_RUN_IFELSE([AC_LANG_PROGRAM(
[[#include <fstream>
#include <string>]],
[[unsigned bx, cx, dx;
unsigned volatile res = 0;
std::string str("");
__asm__ __volatile__ ("cpuid" : "=c"(cx), "=d"(dx) : "a"(0x1));
if (dx & (1 << 25)) { res |= (1 <<  0); str += "sse"; }
if (dx & (1 << 26)) { res |= (1 <<  1); str += " sse2"; }
if (cx & (1 <<  0)) { res |= (1 <<  2); str += " sse3"; }
if (cx & (1 <<  9)) { res |= (1 <<  3); str += " ssse3"; }
if (cx & (1 << 19)) { res |= (1 <<  4); str += " sse4.1"; }
if (cx & (1 << 20)) { res |= (1 <<  5); str += " sse4.2"; }
if (cx & (1 << 28)) { res |= (1 <<  7); str += " avx1"; }
if (cx & (1 << 12)) { res |= (1 <<  8); str += " fma"; }
__asm__ __volatile__ ("cpuid" : "=c"(cx) : "a"(0x80000001));
if (cx & (1 <<  6)) { res |= (1 <<  6); str += " sse4a"; }
__asm__ __volatile__ ("cpuid" : "=b"(bx) : "a"(0x7), "c"(0x0));
if (bx & (1 <<  5)) { res |= (1 <<  9); str += " avx2"; }
if (bx & (1 <<  4)) { res |= (1 << 10); str += " tsxhle"; }
if (bx & (1 << 11)) { res |= (1 << 11); str += " tsxrtm"; }
if (bx & (1 << 16)) { res |= (1 << 12); str += " avx512f"; }
if (bx & (1 << 28)) { res |= (1 << 13); str += " avx512cd"; }
if (bx & (1 << 17)) { res |= (1 << 14); str += " avx512dq"; }
if (bx & (1 << 30)) { res |= (1 << 15); str += " avx512bw"; }
if (bx & (1 << 21)) { res |= (1 << 16); str += " avx512fma"; }
if (bx & (1 <<  1)) { res |= (1 << 17); str += " avx512vbmi"; }
std::ofstream ofile("conftest.data");
ofile << str << " featurebits=" << res;
ofile.close();]])],
                  [ch_stat_cpuid=yes])
    if test "x$ch_stat_cpuid" = "xyes" && test -r "conftest.data"; then
      ch_str_veccpu_feature=`cat conftest.data`
      ch_str_have_cpusse=`expr "x$ch_str_veccpu_feature" : 'x.*\(sse\) '`
      ch_str_have_cpusse2=`expr "x$ch_str_veccpu_feature" : 'x.*\(sse2\) '`
      ch_str_have_cpusse3=`expr "x$ch_str_veccpu_feature" : 'x.* \(sse3\) '`
      ch_str_have_cpussse3=`expr "x$ch_str_veccpu_feature" : 'x.*\(ssse3\) '`
      ch_str_have_cpusse4d1=`expr "x$ch_str_veccpu_feature" : 'x.*\(sse4\.1\) '`
      ch_str_have_cpusse4d2=`expr "x$ch_str_veccpu_feature" : 'x.*\(sse4\.2\) '`
      ch_str_have_cpusse4a=`expr "x$ch_str_veccpu_feature" : 'x.*\(sse4a\) '`
      ch_str_have_cpuavx1=`expr "x$ch_str_veccpu_feature" : 'x.*\(avx1\) '`
      ch_str_have_cpufma=`expr "x$ch_str_veccpu_feature" : 'x.*\(fma\) '`
      ch_str_have_cpuavx2=`expr "x$ch_str_veccpu_feature" : 'x.*\(avx2\) '`
      ch_str_have_cputsxhle=`expr "x$ch_str_veccpu_feature" : 'x.*\(tsxhle\) '`
      ch_str_have_cputsxrtm=`expr "x$ch_str_veccpu_feature" : 'x.*\(tsxrtm\) '`
      ch_str_have_cpuavx512f=`expr "x$ch_str_veccpu_feature" : 'x.*\(avx512f\) '`
      ch_str_have_cpuavx512cd=`expr "x$ch_str_veccpu_feature" : 'x.*\(avx512cd\) '`
      ch_str_have_cpuavx512dq=`expr "x$ch_str_veccpu_feature" : 'x.*\(avx512dq\) '`
      ch_str_have_cpuavx512bw=`expr "x$ch_str_veccpu_feature" : 'x.*\(avx512bw\) '`
      ch_str_have_cpuavx512fma=`expr "x$ch_str_veccpu_feature" : 'x.*\(avx512fma\) '`
      ch_str_have_cpuavx512vbmi=`expr "x$ch_str_veccpu_feature" : 'x.*\(avx512vbmi\) '`
      [ch_int_veccpu_featurebits=`expr "x$ch_str_veccpu_feature" : 'x.*featurebits=\([0-9]*\)'`]
#DBG      echo "$ch_str_have_cpusse x $ch_str_have_cpusse2 x $ch_str_have_cpusse3 x $ch_str_have_cpussse3 x $ch_str_have_cpusse4d1 x $ch_str_have_cpusse4d2 x $ch_str_have_cpusse4a x $ch_str_have_cpuavx1 x $ch_str_have_cpufma x $ch_str_have_cpuavx2 x $ch_int_veccpu_featurebits"
      ch_str_veccpu_feature=`expr "x$ch_str_veccpu_feature" : 'x\(.*\) featurebits=.*'`

#--Test compiler
#--SSE

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-msse" ;;
      esac
      if test "x$ch_str_have_cpusse" = "xsse"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <xmmintrin.h>  // SSE]],
[[__m128 test_v4sf = _mm_set1_ps(1.f) + _mm_setr_ps(1.f, 2.f, 3.f, 4.f);]])],
                          [ch_str_sse_h="xmmintrin.h"])
        if test "x$ch_str_sse_h" = "xxmmintrin.h"; then
          if test "x$ch_stat_vecopp" != "xno"; then
            ch_stat_vecopp="yes"
          fi
        else
          # Vector ops do not work
          ch_stat_vecopp="no"
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <xmmintrin.h>  // SSE]],
[[__m128 test_v4sf = _mm_add_ps(_mm_set1_ps(1.f),
                   _mm_setr_ps(1.f, 2.f, 3.f, 4.f));]])],
                            [ch_str_sse_h="xmmintrin.h"])
        fi
        if test "x$ch_str_sse_h" = "xxmmintrin.h"; then
          ch_str_veccxx_feature="sse"
          ch_int_veccxx_featurebits=1
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // General Vec]],
[[__m128 test_v4sf = _mm_add_ps(_mm_set1_ps(1.f),
                              _mm_setr_ps(1.f, 2.f, 3.f, 4.f));]])],
                            [ch_str_intelvec_h="immintrin.h"],
                            [ch_str_intelvec_h="unknown"])
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_sse_h $ch_str_intelvec_h"

#--SSE2

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-msse2" ;;
      esac
      if test "x$ch_str_have_cpusse2" = "xsse2"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <emmintrin.h>  // SSE2]],
[[__m128d test_v2df = _mm_set1_pd(1.) + _mm_setr_pd(1., 2.);]])],
                          [ch_str_sse2_h="emmintrin.h"])
        if test "x$ch_str_sse2_h" = "xemmintrin.h"; then
          if test "x$ch_stat_vecopp" != "xno"; then
            ch_stat_vecopp="yes"
          fi
        else
          # Vector ops do not work
          ch_stat_vecopp="no"
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <emmintrin.h>  // SSE2]],
[[__m128d test_v2df = _mm_add_pd(_mm_set1_pd(1.), _mm_setr_pd(1., 2.));]])],
                            [ch_str_sse2_h="emmintrin.h"])
        fi
        if test "x$ch_str_sse2_h" = "xemmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature sse2"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '2'`
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // General Vec]],
[[__m128d test_v2df = _mm_add_pd(_mm_set1_pd(1.), _mm_setr_pd(1., 2.));]])],
                            [ch_str_intelvec_h="immintrin.h"],
                            [ch_str_intelvec_h="unknown"])
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_sse2_h $ch_str_intelvec_h"

#--SSE3

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-msse3" ;;
      esac
      if test "x$ch_str_have_cpusse3" = "xsse3"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <pmmintrin.h>  // SSE3]],
[[__m128 test_v4sf = _mm_hadd_ps(_mm_set1_ps(1.f),
                               _mm_setr_ps(1.f, 2.f, 3.f, 4.f));]])],
                          [ch_str_sse3_h="pmmintrin.h"])
        if test "x$ch_str_sse3_h" = "xpmmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature sse3"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '4'`
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // General Vec]],
[[__m128 test_v4sf = _mm_hadd_ps(_mm_set1_ps(1.f),
                               _mm_setr_ps(1.f, 2.f, 3.f, 4.f));]])],
                            [ch_str_intelvec_h="immintrin.h"],
                            [ch_str_intelvec_h="unknown"])
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_sse3_h $ch_str_intelvec_h"

#--SSSE3

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-mssse3" ;;
      esac
      if test "x$ch_str_have_cpussse3" = "xssse3"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <tmmintrin.h>  // SSSE3]],
[[__m128i test_v4si = _mm_sign_epi32(_mm_set1_epi32(2),
                                   _mm_setr_epi32(1, -1, 1, -1));]])],
                          [ch_str_ssse3_h="tmmintrin.h"])
        if test "x$ch_str_ssse3_h" = "xtmmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature ssse3"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '8'`
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // General Vec]],
[[__m128i test_v4si = _mm_sign_epi32(_mm_set1_epi32(2),
                                   _mm_setr_epi32(1, -1, 1, -1));]])],
                            [ch_str_intelvec_h="immintrin.h"],
                            [ch_str_intelvec_h="unknown"])
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_ssse3_h $ch_str_intelvec_h"

#--SSE4.1

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-msse4.1" ;;
      esac
      if test "x$ch_str_have_cpusse4d1" = "xsse4.1"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <smmintrin.h>  // SSE4.1]],
[[__m128 test_v4sf = _mm_dp_ps(_mm_set1_ps(1.f),
                             _mm_setr_ps(1.f, 2.f, 3.f, 0.f),
                             113);]])],
                          [ch_str_sse4d1_h="smmintrin.h"])
        if test "x$ch_str_sse4d1_h" = "xsmmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature sse4.1"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '16'`
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // General Vec]],
[[__m128 test_v4sf = _mm_dp_ps(_mm_set1_ps(1.f),
                             _mm_setr_ps(1.f, 2.f, 3.f, 0.f),
                             113);]])],
                            [ch_str_intelvec_h="immintrin.h"],
                            [ch_str_intelvec_h="unknown"])
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_sse4d1_h $ch_str_intelvec_h"

#--SSE4.2

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-msse4.2" ;;
      esac
      if test "x$ch_str_have_cpusse4d2" = "xsse4.2"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <nmmintrin.h>  // SSE4.2]],
[[__m128i test_v2di = _mm_cmpgt_epi64(_mm_set1_epi32(1),
                                    _mm_setr_epi32(1, 2, 3, 4));]])],
                          [ch_str_sse4d2_h="nmmintrin.h"])
        if test "x$ch_str_sse4d2_h" = "xnmmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature sse4.2"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '32'`
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // General Vec]],
[[__m128i test_v2di = _mm_cmpgt_epi64(_mm_set1_epi32(1),
                                    _mm_setr_epi32(1, 2, 3, 4));]])],
                            [ch_str_intelvec_h="immintrin.h"],
                            [ch_str_intelvec_h="unknown"])
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_sse4d2_h $ch_str_intelvec_h"

#--SSE4a

      case $ch_str_cxxmake in
        # Intel likes to pretend that AMD extensions don't exist :P
        intel   ) CXXFLAGS="" ;;
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-msse4a" ;;
      esac
      if test "x$ch_str_have_cpusse4a" = "xsse4a"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <ammintrin.h>  // SSE4a]],
[[__m128i test_v2di = _mm_extracti_si64(_mm_set1_epi32(6), 1, 2);]])],
                          [ch_str_sse4a_h="ammintrin.h"])
        if test "x$ch_str_sse4a_h" = "xammintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature sse4a"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '64'`
          # The general "immintrin.h" doesn't seem to work with SSE4a
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_sse4a_h $ch_str_intelvec_h"

#--AVX

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-mavx" ;;
      esac
      if test "x$ch_str_have_cpuavx1" = "xavx1"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // AVX]],
[[__m256d test_v4df = _mm256_set1_pd(1.) + _mm256_setr_pd(1., 2., 3., 4.);]])],
                          [ch_str_avx_h="immintrin.h"])
        if test "x$ch_str_avx_h" = "ximmintrin.h"; then
          if test "x$ch_stat_vecopp" != "xno"; then
            ch_stat_vecopp="yes"
          fi
        else
          # Vector ops do not work
          ch_stat_vecopp="no"
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // AVX]],
[[__m256d test_v4df = _mm256_add_pd(_mm256_set1_pd(1.),
                                  _mm256_setr_pd(1., 2., 3., 4.));]])],
                            [ch_str_avx_h="immintrin.h"])
        fi
        if test "x$ch_str_avx_h" = "ximmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature avx"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '128'`
          # Seems we can only use immintrin.h reliably
          ch_str_intelvec_h="immintrin.h"
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_avx_h $ch_str_intelvec_h"

#--FMA

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-mfma" ;;
      esac
      if test "x$ch_str_have_cpufma" = "xfma"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // FMA]],
[[__m256d test_v4df = _mm256_fmadd_pd(_mm256_set1_pd(2.), _mm256_set1_pd(3.),
                                    _mm256_setr_pd(1., 2., 3., 4.));]])],
                          [ch_str_fma_h="immintrin.h"])
        if test "x$ch_str_fma_h" = "ximmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature fma"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '256'`
          # Seems we can only use immintrin.h reliably
          ch_str_intelvec_h="immintrin.h"
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_fma_h $ch_str_intelvec_h"

#--AVX2

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-mavx2" ;;
      esac
      if test "x$ch_str_have_cpuavx2" = "xavx2"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // AVX2]],
[[__m256i test_v8si =
  _mm256_set1_epi32(1) + _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8);]])],
                          [ch_str_avx2_h="immintrin.h"])
        if test "x$ch_str_avx2_h" = "ximmintrin.h"; then
          if test "x$ch_stat_vecopp" != "xno"; then
            ch_stat_vecopp="yes"
          fi
        else
          # Vector ops do not work
          ch_stat_vecopp="no"
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <immintrin.h>  // AVX2]],
[[__m256i test_v8si = _mm256_add_epi32(
  _mm256_set1_epi32(1), _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8));]])],
                            [ch_str_avx2_h="immintrin.h"])
        fi
        if test "x$ch_str_avx2_h" = "ximmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature avx2"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '512'`
          # Seems we can only use immintrin.h reliably
          ch_str_intelvec_h="immintrin.h"
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_avx2_h $ch_str_intelvec_h"

#--TSX HLE

      if test "x$ch_str_have_cputsxhle" = "xtsxhle"; then
        AC_MSG_WARN([The CPU supports TSX-HLE but the compiler has not been tested to see if it supports this feature])
      fi

#--TSX RTM

      if test "x$ch_str_have_cputsxrtm" = "xtsxrtm"; then
        AC_MSG_WARN([The CPU supports TSX-RTM but the compiler has not been tested to see if it supports this feature])
      fi

#--AVX512F foundation (unvalidated)

      case $ch_str_cxxmake in
        portland) CXXFLAGS="-Mvect=sse" ;;
        *       ) CXXFLAGS="-mavx512f" ;;
      esac
      if test "x$ch_str_have_cpuavx512f" = "xavx512f"; then
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <zmmintrin.h>  // AVX512F]],
[[__m512d test_v8df = _mm512_set1_pd(1.) + _mm512_setr_pd(1., 2., 3., 4., 5., 6., 7., 8.);]])],
                          [ch_str_avx512f_h="zmmintrin.h"])
        if test "x$ch_str_avx512f_h" = "xzmmintrin.h"; then
          if test "x$ch_stat_vecopp" != "xno"; then
            ch_stat_vecopp="yes"
          fi
        else
          # Vector ops do not work
          ch_stat_vecopp="no"
          AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <zmmintrin.h>  // AVX512F]],
[[__m512d test_v8df = _mm512_add_pd(_mm512_set1_pd(1.),
                                  _mm512_setr_pd(1., 2., 3., 4., 5., 6., 7., 8.));]])],
                            [ch_str_avx512f_h="zmmintrin.h"])
        fi
        if test "x$ch_str_avx512f_h" = "xzmmintrin.h"; then
          ch_str_veccxx_feature="$ch_str_veccxx_feature avx512f"
          ch_int_veccxx_featurebits=`expr "$ch_int_veccxx_featurebits" + '4096'`
          # ch_str_intelvec_h="zmmintrin.h" Is this general??
        fi
      fi
#DBG      echo "Values: $ch_str_veccxx_feature $ch_int_veccxx_featurebits $ch_str_avx512f_h $ch_str_intelvec_h"


#--AVX512CD

      if test "x$ch_str_have_cpuavx512cd" = "xavx512cd"; then
        AC_MSG_WARN([The CPU supports avx512cd but the compiler has not been tested to see if it supports this feature])
      fi

#--AVX512DQ

      if test "x$ch_str_have_cpuavx512dq" = "xavx512dq"; then
        AC_MSG_WARN([The CPU supports avx512dq but the compiler has not been tested to see if it supports this feature])
      fi

#--AVX512BW

      if test "x$ch_str_have_cpuavx512bw" = "xavx512bw"; then
        AC_MSG_WARN([The CPU supports avx512bw but the compiler has not been tested to see if it supports this feature])
      fi

#--AVX512FMA

      if test "x$ch_str_have_cpuavx512fma" = "xavx512fma"; then
        AC_MSG_WARN([The CPU supports avx512fma but the compiler has not been tested to see if it supports this feature])
      fi

#--AVX512VBMI

      if test "x$ch_str_have_cpuavx512vbmi" = "xavx512vbmi"; then
        AC_MSG_WARN([The CPU supports avx512vbmi but the compiler has not been tested to see if it supports this feature])
      fi
    fi

# Reset flags we changed
    CXXFLAGS="$ch_str_savecxxflags"

#--Results

    AC_MSG_RESULT([$ch_str_veccxx_feature])
    if test "$ch_int_veccpu_featurebits" -ne "$ch_int_veccxx_featurebits"; then
      if test "$ch_int_veccxx_featurebits" -eq "0"; then
        AC_MSG_NOTICE([The CPU supports the vector extensions: $ch_str_veccpu_feature, but compiler intrinsics do not work])
      else
        AC_MSG_NOTICE([The CPU supports the vector extensions: $ch_str_veccpu_feature, but compiler intrinsics only work for: $ch_str_veccxx_feature])
      fi
    fi
    if test "$ch_int_veccpu_featurebits" -ne "0"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_CPU_BITS],
                         [$ch_int_veccpu_featurebits])
    fi
    if test "$ch_int_veccxx_featurebits" -ne "0"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS],
                         [$ch_int_veccxx_featurebits])
    fi
    if test "x$ch_str_sse_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSE_H], [<$ch_str_sse_h>])
    fi
    if test "x$ch_str_sse2_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSE2_H], [<$ch_str_sse2_h>])
    fi
    if test "x$ch_str_sse3_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSE3_H], [<$ch_str_sse3_h>])
    fi
    if test "x$ch_str_ssse3_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSSE3_H], [<$ch_str_ssse3_h>])
    fi
    if test "x$ch_str_sse4d1_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSE4D1_H], [<$ch_str_sse4d1_h>])
    fi
    if test "x$ch_str_sse4d2_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSE4D2_H], [<$ch_str_sse4d2_h>])
    fi
    if test "x$ch_str_sse4a_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_SSE4A_H], [<$ch_str_sse4a_h>])
    fi
    if test "x$ch_str_avx_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_AVX_H], [<$ch_str_avx_h>])
    fi
    if test "x$ch_str_fma_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_FMA_H], [<$ch_str_fma_h>])
    fi
    if test "x$ch_str_avx2_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_AVX2_H], [<$ch_str_avx2_h>])
    fi
    if test "x$ch_str_avx512f_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_AVX512F_H], [<$ch_str_avx512f_h>])
    fi
    if test "x$ch_str_intelvec_h" != "xunknown"; then
      AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_X86VECEXT_INTEL_H], [<$ch_str_intelvec_h>])
    fi
    if test "x$ch_stat_vecopp" = "xyes"; then
      AC_DEFINE([CHDEF_SYSTEM_X86VECEXT_HAVE_OPERATOR])
    fi
    ;;

  *    )
    AC_MSG_RESULT([not an x86 processor])
    ;;
esac

#--------------------------------------------------------------------
# Check for a hardware popcnt instruction and intrinsic
#--------------------------------------------------------------------

AC_MSG_CHECKING([for hardware popcnt])

ch_int_popcnt_method=0
ch_str_popcnt_func=unknown
ch_str_popcnt_h=unknown

case $host_cpu in
  *86*)

#--Test hardware

    ch_stat_cpuid=no
    # Avoid using an optimization flags -- the assembly can cause problems.
    # We will also gradually add sse flags
    ch_str_savecxxflags="$CXXFLAGS"
    CXXFLAGS=""
    AC_RUN_IFELSE([AC_LANG_PROGRAM(
[[#include <fstream>]],
[[unsigned cx;
unsigned volatile res = 0;
__asm__ __volatile__ ("cpuid" : "=c"(cx) : "a"(0x1));
if (cx & (1 << 23)) { res |= (1 << 0); }
std::ofstream ofile("conftest.data");
ofile << " featurebits=" << res;
ofile.close();]])],
                  [ch_stat_cpuid="yes"])
    if test "x$ch_stat_cpuid" = "xyes" && test -r "conftest.data"; then
      ch_str_popcntcpu_feature=`cat conftest.data`
      [ch_int_popcntcpu_featurebits=`expr "x$ch_str_popcntcpu_feature" : 'x.*featurebits=\([0-9]*\)'`]

#     'popcnt' intrinscis are a bit cumbersome.  With newer versions of gcc, it
#     can be accessed with the _mm_popcnt_u32 intrinsic which is available
#     through popcntintrin.h, in earlier versions through smmintrin.h  It seems
#     easier to just use assembly.

      if test "$ch_int_popcntcpu_featurebits" -eq "1"; then

#--Try assembly

        AC_RUN_IFELSE([AC_LANG_PROGRAM(
[[#include <fstream>]],
[[unsigned x = 0x3;
int ret;
__asm__ ("popcnt %1, %0" : "=r" (ret) : "r" (x));
std::ofstream ofile("conftest.data");  // IO in case of optimization
ofile << ret;
ofile.close();]])],
                      [ch_int_popcnt_method=1])
      fi
    fi

# Reset flags we changed
    CXXFLAGS="$ch_str_savecxxflags"
    ;;
esac

#--Results

if test "$ch_int_popcnt_method" -gt "0"; then
  AC_MSG_RESULT([yes])
  AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_POPCNT_METHOD], [$ch_int_popcnt_method])
else
  AC_MSG_RESULT([no])
fi

#--------------------------------------------------------------------
# Check for CUDA
#--------------------------------------------------------------------

ch_stat_cuda=unknown
#ch_str_nvcc=unknown
ch_str_nvccpath=system
ch_str_cudatkdir=unknown
ch_str_cudahdir=unknown
ch_int_cudadevicecount=0
ch_str_cudaattributes=unknown
ch_str_cudaattributesprint=""
ch_int_cudaattr_computemajor=0
ch_int_cudaattr_computeminor=0

# nvcc should normally be in the path
AC_ARG_WITH([nvcc],
            [AC_HELP_STRING([--with-nvcc=<path>],
                            [path to Nvidia nvcc CUDA compiler])],
                            [ch_str_nvccpath=$withval])
if test "x$ch_str_nvccpath" = "xsystem"; then
  AC_CHECK_PROG([ch_str_nvcc], [nvcc], [nvcc], [unknown])
else  
  AC_PATH_PROG([ch_str_nvcc], [nvcc], [unknown], [$ch_str_nvccpath])
fi
if test "x$ch_str_nvcc" != "xunknown"; then
  ch_stat_cuda="havenvcc"
  # Cuda is normally installed in the system at /usr/lib/
  AC_CHECK_LIB([cuda], [cuInit], [ch_stat_cuda="havelib"])
  if test "x$ch_stat_cuda" = "xhavelib"; then
    # Normally installed to /usr/local/cuda
    AC_ARG_WITH([cudatk],
                [AC_HELP_STRING([--with-cuda=<path>],
                                [path to CUDA toolkit installation])],
                [ch_str_cudatkdir=$withval])
    if test "x$ch_str_cudatkdir" = "xunknown"; then
      AC_CHECK_HEADER([cuda.h],
                      [ch_str_cudahdir="system"])
      if test "x$ch_str_cudahdir" = "xunknown"; then
        AC_MSG_CHECKING([standard locations for cuda.h])
        if test -f "/usr/local/cuda/include/cuda.h"; then
          ch_str_cudahdir="/usr/local/cuda/include"
          AC_MSG_RESULT([$ch_str_cudahdir])
        else
          AC_MSG_RESULT([no])
        fi
      fi
    else
      if test -f "$ch_str_cudatkdir/include/cuda.h"; then
        ch_str_cudahdir="$ch_str_cudatkdir/include"
      fi
    fi
    if test "x$ch_str_cudahdir" = "xunknown"; then
      AC_MSG_WARN([cannot find cuda.h:  use --with-cudatk=<path>])
    else

#--Have nvcc, lib, and header file

      ch_stat_cuda="yes"
      # Avoid using an optimization flags -- the assembly can cause problems.
      # We will also gradually add sse flags
      ch_str_savecxxflags="$CXXFLAGS"
      ch_str_savecppflags="$CPPFLAGS"
      CXXFLAGS="-lcuda"
      if test "x$ch_str_cudahdir" != "xsystem"; then
         CPPFLAGS="-I$ch_str_cudahdir"
      fi

#--Just run init

      AC_RUN_IFELSE([AC_LANG_PROGRAM(
[[#include <cuda.h>]],
[[CUresult cuer = cuInit(0);
return cuer;]])],
                    [ch_stat_cuda="yes"],
                    [ch_stat_cuda="failed"])
      if test "x$ch_stat_cuda" = "xfailed"; then
        AC_MSG_WARN([The CUDA toolkit was found but a simple program failed.])
      else

#--Look for CUDA devices

        AC_MSG_CHECKING([number of CUDA devices])
        AC_RUN_IFELSE([AC_LANG_PROGRAM(
[[#include <fstream>
#include <cuda.h>]],
[[CUresult cuer = cuInit(0);
if (cuer) return 1;
int deviceCount;
cuer = cuDeviceGetCount(&deviceCount);
if (cuer) return 1;
std::ofstream ofile("conftest.data");
ofile << deviceCount;
ofile.close();
return 0;]])],
                    [ch_stat_cuda="yes"],
                    [ch_stat_cuda="failed"])
        if test "x$ch_stat_cuda" = "xyes"; then
          ch_int_cudadevicecount=`cat conftest.data`
        fi
        AC_MSG_RESULT([$ch_int_cudadevicecount])
        if test "$ch_int_cudadevicecount" -eq "0"; then
          AC_MSG_NOTICE([No CUDA devices found])
        else
          if test "$ch_int_cudadevicecount" -gt "1"; then
            AC_MSG_NOTICE([More than 1 CUDA device found.  Properties taken from first device])
          fi

#--Name and other attributes.  This test program reproduces in a similar form
#--output from Nvidia's deviceQueryDrv sample

          AC_RUN_IFELSE([AC_LANG_PROGRAM(
[[#include <fstream>
#include <iomanip>
#include <cuda.h>]],
[[CUresult cuer = cuInit(0);
if (cuer) return 1;
int dev = 0;

int driverVersion;
char deviceName[256];
int major;
int minor;
size_t totalGlobalMem;
int multiProcessorCount;
int totalConstantMem;
int sharedMemPerBlock;
int registersPerBlock;
int warpSize;
int maxThreadsPerMP;
int maxThreadsPerBlock;
int maxBlockDimX;
int maxBlockDimY;
int maxBlockDimZ;
int maxGridDimX;
int maxGridDimY;
int maxGridDimZ;

cuer = cuDeviceComputeCapability(&major, &minor, dev);
if (cuer) return 1;
cuer = cuDeviceGetName(deviceName, 256, dev);
if (cuer) return 1;
cuer = cuDriverGetVersion(&driverVersion);
if (cuer) return 1;
cuer = cuDeviceTotalMem(&totalGlobalMem, dev);
if (cuer) return 1;

cuer = cuDeviceGetAttribute(&multiProcessorCount,
                            CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&totalConstantMem,
                            CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&sharedMemPerBlock,
                            CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK,
                            dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&registersPerBlock,
                            CU_DEVICE_ATTRIBUTE_MAX_REGISTERS_PER_BLOCK, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&warpSize,
                            CU_DEVICE_ATTRIBUTE_WARP_SIZE, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxThreadsPerMP,
                            CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_MULTIPROCESSOR,
                            dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxThreadsPerBlock,
                            CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxBlockDimX,
                            CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxBlockDimY,
                            CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Y, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxBlockDimZ,
                            CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Z, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxGridDimX,
                            CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_X, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxGridDimY,
                            CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Y, dev);
if (cuer) return 1;
cuer = cuDeviceGetAttribute(&maxGridDimZ,
                            CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Z, dev);
if (cuer) return 1;

std::ofstream ofile("conftest.data");
ofile << "DRIVER_VERSION=" << driverVersion << '=';
ofile << " DEVICE_NAME=" << deviceName << '=';
ofile << " COMPUTE_MAJOR=" << major << '=';
ofile << " COMPUTE_MINOR=" << minor << '=';
ofile << " TOTAL_GLOBAL_MEMORY=" << totalGlobalMem << '=';
ofile << " MULTIPROCESSOR_COUNT=" << multiProcessorCount << '=';
ofile << " TOTAL_CONSTANT_MEMORY=" << totalConstantMem << '=';
ofile << " MAX_SHARED_MEMORY_PER_BLOCK=" << sharedMemPerBlock << '=';
ofile << " MAX_REGISTERS_PER_BLOCK=" << registersPerBlock << '=';
ofile << " WARP_SIZE=" << warpSize << '=';
ofile << " MAX_THREADS_PER_MULTIPROCESSOR=" << maxThreadsPerMP << '=';
ofile << " MAX_THREADS_PER_BLOCK=" << maxThreadsPerBlock << '=';
ofile << " MAX_BLOCK_DIM_X=" << maxBlockDimX << '=';
ofile << " MAX_BLOCK_DIM_Y=" << maxBlockDimY << '=';
ofile << " MAX_BLOCK_DIM_Z=" << maxBlockDimZ << '=';
ofile << " MAX_GRID_DIM_X=" << maxGridDimX << '=';
ofile << " MAX_GRID_DIM_Y=" << maxGridDimY << '=';
ofile << " MAX_GRID_DIM_Z=" << maxGridDimZ << '=';
ofile.close();
const char *const indent = "  ";
const int labelWidth = 47;
std::ofstream pfile("conftest.print");
pfile << indent << std::left << std::setw(labelWidth)
      << "CUDA driver version:"
      << driverVersion/1000 << '.' << (driverVersion%100)/10
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Device name:"
      << deviceName
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "CUDA compute capability:"
      << major << '.' << minor
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Total amount of global memory:"
      << totalGlobalMem << " bytes"
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Number of multiprocessors:"
      << multiProcessorCount
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Total amount of constant memory:"
      << totalConstantMem << " bytes"
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Total amount of shared memory per block:"
      << sharedMemPerBlock << " bytes"
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Total number of registers available per block:"
      << registersPerBlock
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Warp size:"
      << warpSize
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Maximum number of threads per multiprocessor:"
      << maxThreadsPerMP
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Maximum number of threads per block:"
      << maxThreadsPerBlock
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Maximum sizes of each dimension of a block:"
      << maxBlockDimX << " x " << maxBlockDimY << " x " << maxBlockDimZ
      << std::endl;
pfile << indent << std::left << std::setw(labelWidth)
      << "Maximum sizes of each dimension of a grid:"
      << maxGridDimX << " x " << maxGridDimY << " x " << maxGridDimZ
      << std::endl;
pfile.close();
return 0;]])],
                    [ch_stat_cuda="yes"],
                    [ch_stat_cuda="failed"])
          if test "x$ch_stat_cuda" = "xyes"; then
            ch_str_cudaattributes=`cat conftest.data`
            ch_str_cudaattributesprint=`cat conftest.print`
            # Compute codes
            [ch_int_cudaattr_computemajor=`expr "x$ch_str_cudaattributes" : 'x.*COMPUTE_MAJOR=\([^=]*\)'`]
            [ch_int_cudaattr_computeminor=`expr "x$ch_str_cudaattributes" : 'x.*COMPUTE_MINOR=\([^=]*\)'`]
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDA_DRIVER_VERSION],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*DRIVER_VERSION=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDA_DEVICE_NAME],
                               [["`expr "x$ch_str_cudaattributes" : 'x.*DEVICE_NAME=\([^=]*\)'`"]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDA_COMPUTE_CAPABILITY],
                               [$ch_int_cudaattr_computemajor.$ch_int_cudaattr_computeminor])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_TOTAL_GLOBAL_MEMORY],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*TOTAL_GLOBAL_MEMORY=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MULTIPROCESSOR_COUNT],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MULTIPROCESSOR_COUNT=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_TOTAL_CONSTANT_MEMORY],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*TOTAL_CONSTANT_MEMORY=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_SHARED_MEMORY_PER_BLOCK],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_SHARED_MEMORY_PER_BLOCK=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_REGISTERS_PER_BLOCK],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_REGISTERS_PER_BLOCK=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_WARP_SIZE],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*WARP_SIZE=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_THREADS_PER_MULTIPROCESSOR],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_THREADS_PER_MULTIPROCESSOR=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_THREADS_PER_BLOCK],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_THREADS_PER_BLOCK=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_BLOCK_DIM_X],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_BLOCK_DIM_X=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_BLOCK_DIM_Y],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_BLOCK_DIM_Y=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_BLOCK_DIM_Z],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_BLOCK_DIM_Z=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_GRID_DIM_X],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_GRID_DIM_X=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_GRID_DIM_Y],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_GRID_DIM_Y=\([^=]*\)'`]])
            AC_DEFINE_UNQUOTED([CHDEF_SYSTEM_CUDAATTR_MAX_GRID_DIM_Z],
                               [[`expr "x$ch_str_cudaattributes" : 'x.*MAX_GRID_DIM_Z=\([^=]*\)'`]])
          fi
# end have device
        fi
# end init passed
      fi
      CXXFLAGS="$ch_str_savecxxflags"
      CPPFLAGS="$ch_str_savecppflags"
# end have header
    fi
# end have lib
  fi
# end have nvcc
fi

#--------------------------------------------------------------------
# Finalize
#--------------------------------------------------------------------

#--Output to Make.defs.config_defaults

XTRACXXFLAGS=""
XTRACPPFLAGS=""
XTRALDFLAGS=""
NVCC=""
NVCCFLAGS=""
PTXASFLAGS=""
ch_stat_haveflags="no"

# Vector registers with gcc
if test "$ch_int_veccxx_featurebits" -ne "0"; then
  ch_str_savecxxflags="$CXXFLAGS"
  ch_stat_havegenvecflag="no"
  case "x$ch_str_cxxmake" in
    "xgcc")
      # Try compiling an sse2 example with '-march=native' to see if the option
      # is accepted by the compiler
      CXXFLAGS="-march=native"
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[#include <emmintrin.h>  // SSE2]],
[[__m128d test_v2df = _mm_add_pd(_mm_set1_pd(1.), _mm_setr_pd(1., 2.));]])],
                        [XTRACXXFLAGS="$XTRACXXFLAGS -march=native"
                         ch_stat_havegenvecflag="yes"])
    ;;
  esac
  if test "x$ch_stat_havegenvecflag" = "xyes"; then
    ch_stat_haveflags="yes"
  else
    AC_MSG_NOTICE([Unable to find generic compiler option for using vector registers])
  fi
  CXXFLAGS="$ch_str_savecxxflags"
fi

# Cuda
if test "x$ch_stat_cuda" = "xyes"; then
  NVCC=" $ch_str_nvcc"
  XTRACPPFLAGS="$XTRACPPFLAGS -I$ch_str_cudahdir"
  XTRALDFLAGS="$XTRALDFLAGS -lcuda"
  # pointer size
  if test "x$ch_int_systemptrsize" != "xunknown"; then
    NVCCFLAGS="$NVCCFLAGS -m$ch_int_systemptrsize"
  fi
  # capability 2.1 maps to code 2.1 and arch 2.0.  As of 3.5, all others map
  # directly.  Note that arch must be virtual (compute_) if code is specified
  if test "$ch_int_cudaattr_computemajor" -eq "2" && test "$ch_int_cudaattr_computeminor" -eq "1"; then
    ch_str_gpuarch="compute_20"
  else
    ch_str_gpuarch="compute_$ch_int_cudaattr_computemajor$ch_int_cudaattr_computeminor"
  fi
  ch_str_gpucode="sm_$ch_int_cudaattr_computemajor$ch_int_cudaattr_computeminor"
  NVCCFLAGS="$NVCCFLAGS -gencode=arch=$ch_str_gpuarch,code=$ch_str_gpucode"
  PTXASFLAGS="$PTXASFLAGS --ptxas-options=-v"
  ch_stat_haveflags="yes"
fi

AC_SUBST([XTRACXXFLAGS])
AC_SUBST([XTRACPPFLAGS])
AC_SUBST([XTRALDFLAGS])
AC_SUBST([NVCC])
AC_SUBST([NVCCFLAGS])
AC_SUBST([PTXASFLAGS])

AC_OUTPUT

#--------------------------------------------------------------------
# Print useful information
#--------------------------------------------------------------------

# Cuda attributes
if test "x$ch_stat_cuda" = "xyes"; then
  echo ""
  echo "CUDA Information and Attributes for Device 0"
  echo ""
  echo "$ch_str_cudaattributesprint"
fi

if test "x$ch_stat_haveflags" = "xyes"; then
  echo ""
  echo "The following defaults have been assigned to the Chombo make system."
  echo "They can be overridden using = or added to using += in Make.defs.local"
  echo ""
  if test -n "$XTRACXXFLAGS"; then
    echo "XTRACXXFLAGS   =$XTRACXXFLAGS"
  fi
  if test -n "$XTRACPPFLAGS"; then
    echo "XTRACPPFLAGS   =$XTRACPPFLAGS"
  fi
  if test -n "$XTRALDFLAGS"; then
    echo "XTRALDFLAGS    =$XTRALDFLAGS"
  fi
  if test -n "$NVCC"; then
    echo "NVCC           =$NVCC"
  fi
  if test -n "$NVCCFLAGS"; then
    echo "NVCCFLAGS      =$NVCCFLAGS"
  fi
  if test -n "$PTXASFLAGS"; then
    echo "PTXASFLAGS     =$PTXASFLAGS"
  fi
  echo ""
fi
