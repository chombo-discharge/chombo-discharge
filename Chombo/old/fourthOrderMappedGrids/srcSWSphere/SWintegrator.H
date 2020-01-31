#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SWINTEGRATOR_H_
#define _SWINTEGRATOR_H_

#if (CH_SPACEDIM == 1)

#define WHGT   0
#define WVELX  1

#define WVNORM 1

#define WNUM   2
#define WSLOPE 2

#define UHGT  0
#define UMOMX 1

#define UMOMNORM 1

#define UNUM  2
#define UFLUX 2

#define MWAVE 0
#define ZWAVE 1
#define PWAVE 2

#endif

#if (CH_SPACEDIM == 2)

#define WHGT   0
#define WVELX  1
#define WVELY  2

#define WVNORM 1
#define WVTAN 2

#define WNUM   3
#define WSLOPE 3

#define UHGT  0
#define UMOMX 1
#define UMOMY 2

#define UMOMNORM 1
#define UMOMTAN 2

#define UNUM  3
#define UFLUX 3

#define MWAVE 0
#define ZWAVE 1
#define PWAVE 3

#endif

#if (CH_SPACEDIM == 3)

#define WHGT   0
#define WVELX  1
#define WVELY  2
#define WVELZ  3

#define WVNORM 1
#define WVTAN 2
#define WVTAN2 3

#define WNUM   4
#define WSLOPE 4

#define UHGT  0
#define UMOMX 1
#define UMOMY 2
#define UMOMZ 3

#define UMOMNORM 1
#define UMOMTAN 2
#define UMOMTAN2 3

#define UNUM  4
#define UFLUX 4

#define MWAVE 0
#define ZWAVE 1
#define PWAVE 4

#endif

#endif
