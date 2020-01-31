#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef CH_LANG_CC
// Constants used between C++ and Fortran to index the conserved and
// primitive variables and to know the number of each.  Also, the number
// of slopes and fluxes are defined here.
#endif

#ifndef _LGINTEGRATOR_H_
#define _LGINTEGRATOR_H_

#define GAS
#define GRAVITY

#undef  INFINITE_DOMAIN_BC

#if (CH_SPACEDIM == 2)

#define WRHO   0
#define WVELX  1
#define WVELY  2
#define WPRES  3
#define WENTR  4

#define WNUM   (5)
#define WSLOPE (5)

#define URHO  0
#define UMOMX 1
#define UMOMY 2
#define UENG  3
#define UENTR 4

#define UNUM  (5)
#define UFLUX (5)

#define MWAVE 0
#define ZWAVE 1
#define SWAVE 3
#define PWAVE 4

#endif

#if (CH_SPACEDIM == 3)

#define WRHO   0
#define WVELX  1
#define WVELY  2
#define WVELZ  3
#define WPRES  4
#define WENTR  5

#define WNUM   (6)
#define WSLOPE (6)

#define URHO  0
#define UMOMX 1
#define UMOMY 2
#define UMOMZ 3
#define UENG  4
#define UENTR 5

#define UNUM  (6)
#define UFLUX (6)

#define MWAVE 0
#define ZWAVE 1
#define SWAVE 4
#define PWAVE 5

#endif

#endif


