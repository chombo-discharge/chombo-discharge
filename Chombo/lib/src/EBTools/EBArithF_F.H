#ifndef _EBARITHF_F_H_
#define _EBARITHF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDVOLWGTSUM 
#define GUARDVOLWGTSUM 
// Prototype for Fortran procedure VOLWGTSUM ...
//
void FORTRAN_NAME( VOLWGTSUM ,volwgtsum )(
      CHFp_CONST_FRA(src)
      ,CHFp_CONST_FRA(volfrac)
      ,CHFp_CONST_REAL(norm)
      ,CHFp_CONST_REAL(volume)
      ,CHFp_CONST_INT(comp)
      ,CHFp_CONST_INT(pval)
      ,CHFp_CONST_INT(idoreg)
      ,CHFp_CONST_INT(idoirr)
      ,CHFp_BOX(region) );

#define FORT_VOLWGTSUM FORTRAN_NAME( inlineVOLWGTSUM, inlineVOLWGTSUM)
#define FORTNT_VOLWGTSUM FORTRAN_NAME( VOLWGTSUM, volwgtsum)

inline void FORTRAN_NAME(inlineVOLWGTSUM, inlineVOLWGTSUM)(
      CHFp_CONST_FRA(src)
      ,CHFp_CONST_FRA(volfrac)
      ,CHFp_CONST_REAL(norm)
      ,CHFp_CONST_REAL(volume)
      ,CHFp_CONST_INT(comp)
      ,CHFp_CONST_INT(pval)
      ,CHFp_CONST_INT(idoreg)
      ,CHFp_CONST_INT(idoirr)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_VOLWGTSUM");
 FORTRAN_NAME( VOLWGTSUM ,volwgtsum )(
      CHFt_CONST_FRA(src)
      ,CHFt_CONST_FRA(volfrac)
      ,CHFt_CONST_REAL(norm)
      ,CHFt_CONST_REAL(volume)
      ,CHFt_CONST_INT(comp)
      ,CHFt_CONST_INT(pval)
      ,CHFt_CONST_INT(idoreg)
      ,CHFt_CONST_INT(idoirr)
      ,CHFt_BOX(region) );
}
#endif  // GUARDVOLWGTSUM 

#ifndef GUARDADDTWOFAB 
#define GUARDADDTWOFAB 
// Prototype for Fortran procedure ADDTWOFAB ...
//
void FORTRAN_NAME( ADDTWOFAB ,addtwofab )(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(region)
      ,CHFp_INT(srccomp)
      ,CHFp_INT(destcomp)
      ,CHFp_INT(numcomp) );

#define FORT_ADDTWOFAB FORTRAN_NAME( inlineADDTWOFAB, inlineADDTWOFAB)
#define FORTNT_ADDTWOFAB FORTRAN_NAME( ADDTWOFAB, addtwofab)

inline void FORTRAN_NAME(inlineADDTWOFAB, inlineADDTWOFAB)(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(region)
      ,CHFp_INT(srccomp)
      ,CHFp_INT(destcomp)
      ,CHFp_INT(numcomp) )
{
 CH_TIMELEAF("FORT_ADDTWOFAB");
 FORTRAN_NAME( ADDTWOFAB ,addtwofab )(
      CHFt_FRA(dst)
      ,CHFt_CONST_FRA(src)
      ,CHFt_BOX(region)
      ,CHFt_INT(srccomp)
      ,CHFt_INT(destcomp)
      ,CHFt_INT(numcomp) );
}
#endif  // GUARDADDTWOFAB 

#ifndef GUARDSCALEADDTWOFAB 
#define GUARDSCALEADDTWOFAB 
// Prototype for Fortran procedure SCALEADDTWOFAB ...
//
void FORTRAN_NAME( SCALEADDTWOFAB ,scaleaddtwofab )(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(src)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_BOX(region)
      ,CHFp_INT(srccomp)
      ,CHFp_INT(destcomp)
      ,CHFp_INT(numcomp) );

#define FORT_SCALEADDTWOFAB FORTRAN_NAME( inlineSCALEADDTWOFAB, inlineSCALEADDTWOFAB)
#define FORTNT_SCALEADDTWOFAB FORTRAN_NAME( SCALEADDTWOFAB, scaleaddtwofab)

inline void FORTRAN_NAME(inlineSCALEADDTWOFAB, inlineSCALEADDTWOFAB)(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(src)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_BOX(region)
      ,CHFp_INT(srccomp)
      ,CHFp_INT(destcomp)
      ,CHFp_INT(numcomp) )
{
 CH_TIMELEAF("FORT_SCALEADDTWOFAB");
 FORTRAN_NAME( SCALEADDTWOFAB ,scaleaddtwofab )(
      CHFt_FRA(dst)
      ,CHFt_CONST_FRA(src)
      ,CHFt_CONST_REAL(scale)
      ,CHFt_BOX(region)
      ,CHFt_INT(srccomp)
      ,CHFt_INT(destcomp)
      ,CHFt_INT(numcomp) );
}
#endif  // GUARDSCALEADDTWOFAB 

#ifndef GUARDAXBYFAB 
#define GUARDAXBYFAB 
// Prototype for Fortran procedure AXBYFAB ...
//
void FORTRAN_NAME( AXBYFAB ,axbyfab )(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_CONST_REAL(a)
      ,CHFp_CONST_REAL(b)
      ,CHFp_BOX(region)
      ,CHFp_INT(srccomp)
      ,CHFp_INT(destcomp)
      ,CHFp_INT(numcomp) );

#define FORT_AXBYFAB FORTRAN_NAME( inlineAXBYFAB, inlineAXBYFAB)
#define FORTNT_AXBYFAB FORTRAN_NAME( AXBYFAB, axbyfab)

inline void FORTRAN_NAME(inlineAXBYFAB, inlineAXBYFAB)(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_CONST_REAL(a)
      ,CHFp_CONST_REAL(b)
      ,CHFp_BOX(region)
      ,CHFp_INT(srccomp)
      ,CHFp_INT(destcomp)
      ,CHFp_INT(numcomp) )
{
 CH_TIMELEAF("FORT_AXBYFAB");
 FORTRAN_NAME( AXBYFAB ,axbyfab )(
      CHFt_FRA(dst)
      ,CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(y)
      ,CHFt_CONST_REAL(a)
      ,CHFt_CONST_REAL(b)
      ,CHFt_BOX(region)
      ,CHFt_INT(srccomp)
      ,CHFt_INT(destcomp)
      ,CHFt_INT(numcomp) );
}
#endif  // GUARDAXBYFAB 

#ifndef GUARDAXBYFABCOMP 
#define GUARDAXBYFABCOMP 
// Prototype for Fortran procedure AXBYFABCOMP ...
//
void FORTRAN_NAME( AXBYFABCOMP ,axbyfabcomp )(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_CONST_REAL(a)
      ,CHFp_CONST_REAL(b)
      ,CHFp_CONST_INT(destcomp)
      ,CHFp_CONST_INT(xcomp)
      ,CHFp_CONST_INT(ycomp)
      ,CHFp_BOX(region) );

#define FORT_AXBYFABCOMP FORTRAN_NAME( inlineAXBYFABCOMP, inlineAXBYFABCOMP)
#define FORTNT_AXBYFABCOMP FORTRAN_NAME( AXBYFABCOMP, axbyfabcomp)

inline void FORTRAN_NAME(inlineAXBYFABCOMP, inlineAXBYFABCOMP)(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(x)
      ,CHFp_CONST_FRA(y)
      ,CHFp_CONST_REAL(a)
      ,CHFp_CONST_REAL(b)
      ,CHFp_CONST_INT(destcomp)
      ,CHFp_CONST_INT(xcomp)
      ,CHFp_CONST_INT(ycomp)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_AXBYFABCOMP");
 FORTRAN_NAME( AXBYFABCOMP ,axbyfabcomp )(
      CHFt_FRA(dst)
      ,CHFt_CONST_FRA(x)
      ,CHFt_CONST_FRA(y)
      ,CHFt_CONST_REAL(a)
      ,CHFt_CONST_REAL(b)
      ,CHFt_CONST_INT(destcomp)
      ,CHFt_CONST_INT(xcomp)
      ,CHFt_CONST_INT(ycomp)
      ,CHFt_BOX(region) );
}
#endif  // GUARDAXBYFABCOMP 

#ifndef GUARDSUBTRACTTWOFAB 
#define GUARDSUBTRACTTWOFAB 
// Prototype for Fortran procedure SUBTRACTTWOFAB ...
//
void FORTRAN_NAME( SUBTRACTTWOFAB ,subtracttwofab )(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(region)
      ,CHFp_INT(srccomp)
      ,CHFp_INT(destcomp)
      ,CHFp_INT(numcomp) );

#define FORT_SUBTRACTTWOFAB FORTRAN_NAME( inlineSUBTRACTTWOFAB, inlineSUBTRACTTWOFAB)
#define FORTNT_SUBTRACTTWOFAB FORTRAN_NAME( SUBTRACTTWOFAB, subtracttwofab)

inline void FORTRAN_NAME(inlineSUBTRACTTWOFAB, inlineSUBTRACTTWOFAB)(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(region)
      ,CHFp_INT(srccomp)
      ,CHFp_INT(destcomp)
      ,CHFp_INT(numcomp) )
{
 CH_TIMELEAF("FORT_SUBTRACTTWOFAB");
 FORTRAN_NAME( SUBTRACTTWOFAB ,subtracttwofab )(
      CHFt_FRA(dst)
      ,CHFt_CONST_FRA(src)
      ,CHFt_BOX(region)
      ,CHFt_INT(srccomp)
      ,CHFt_INT(destcomp)
      ,CHFt_INT(numcomp) );
}
#endif  // GUARDSUBTRACTTWOFAB 

#ifndef GUARDMULTIPLYTWOFAB 
#define GUARDMULTIPLYTWOFAB 
// Prototype for Fortran procedure MULTIPLYTWOFAB ...
//
void FORTRAN_NAME( MULTIPLYTWOFAB ,multiplytwofab )(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(region)
      ,CHFp_INT(srccomp)
      ,CHFp_INT(destcomp)
      ,CHFp_INT(numcomp) );

#define FORT_MULTIPLYTWOFAB FORTRAN_NAME( inlineMULTIPLYTWOFAB, inlineMULTIPLYTWOFAB)
#define FORTNT_MULTIPLYTWOFAB FORTRAN_NAME( MULTIPLYTWOFAB, multiplytwofab)

inline void FORTRAN_NAME(inlineMULTIPLYTWOFAB, inlineMULTIPLYTWOFAB)(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(region)
      ,CHFp_INT(srccomp)
      ,CHFp_INT(destcomp)
      ,CHFp_INT(numcomp) )
{
 CH_TIMELEAF("FORT_MULTIPLYTWOFAB");
 FORTRAN_NAME( MULTIPLYTWOFAB ,multiplytwofab )(
      CHFt_FRA(dst)
      ,CHFt_CONST_FRA(src)
      ,CHFt_BOX(region)
      ,CHFt_INT(srccomp)
      ,CHFt_INT(destcomp)
      ,CHFt_INT(numcomp) );
}
#endif  // GUARDMULTIPLYTWOFAB 

#ifndef GUARDDIVIDETWOFAB 
#define GUARDDIVIDETWOFAB 
// Prototype for Fortran procedure DIVIDETWOFAB ...
//
void FORTRAN_NAME( DIVIDETWOFAB ,dividetwofab )(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(region)
      ,CHFp_INT(srccomp)
      ,CHFp_INT(destcomp)
      ,CHFp_INT(numcomp) );

#define FORT_DIVIDETWOFAB FORTRAN_NAME( inlineDIVIDETWOFAB, inlineDIVIDETWOFAB)
#define FORTNT_DIVIDETWOFAB FORTRAN_NAME( DIVIDETWOFAB, dividetwofab)

inline void FORTRAN_NAME(inlineDIVIDETWOFAB, inlineDIVIDETWOFAB)(
      CHFp_FRA(dst)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(region)
      ,CHFp_INT(srccomp)
      ,CHFp_INT(destcomp)
      ,CHFp_INT(numcomp) )
{
 CH_TIMELEAF("FORT_DIVIDETWOFAB");
 FORTRAN_NAME( DIVIDETWOFAB ,dividetwofab )(
      CHFt_FRA(dst)
      ,CHFt_CONST_FRA(src)
      ,CHFt_BOX(region)
      ,CHFt_INT(srccomp)
      ,CHFt_INT(destcomp)
      ,CHFt_INT(numcomp) );
}
#endif  // GUARDDIVIDETWOFAB 

#ifndef GUARDSUBTRACTFABR 
#define GUARDSUBTRACTFABR 
// Prototype for Fortran procedure SUBTRACTFABR ...
//
void FORTRAN_NAME( SUBTRACTFABR ,subtractfabr )(
      CHFp_FRA(dst)
      ,CHFp_CONST_REAL(src)
      ,CHFp_BOX(region) );

#define FORT_SUBTRACTFABR FORTRAN_NAME( inlineSUBTRACTFABR, inlineSUBTRACTFABR)
#define FORTNT_SUBTRACTFABR FORTRAN_NAME( SUBTRACTFABR, subtractfabr)

inline void FORTRAN_NAME(inlineSUBTRACTFABR, inlineSUBTRACTFABR)(
      CHFp_FRA(dst)
      ,CHFp_CONST_REAL(src)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_SUBTRACTFABR");
 FORTRAN_NAME( SUBTRACTFABR ,subtractfabr )(
      CHFt_FRA(dst)
      ,CHFt_CONST_REAL(src)
      ,CHFt_BOX(region) );
}
#endif  // GUARDSUBTRACTFABR 

#ifndef GUARDADDFABR 
#define GUARDADDFABR 
// Prototype for Fortran procedure ADDFABR ...
//
void FORTRAN_NAME( ADDFABR ,addfabr )(
      CHFp_FRA(dst)
      ,CHFp_CONST_REAL(src)
      ,CHFp_BOX(region) );

#define FORT_ADDFABR FORTRAN_NAME( inlineADDFABR, inlineADDFABR)
#define FORTNT_ADDFABR FORTRAN_NAME( ADDFABR, addfabr)

inline void FORTRAN_NAME(inlineADDFABR, inlineADDFABR)(
      CHFp_FRA(dst)
      ,CHFp_CONST_REAL(src)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_ADDFABR");
 FORTRAN_NAME( ADDFABR ,addfabr )(
      CHFt_FRA(dst)
      ,CHFt_CONST_REAL(src)
      ,CHFt_BOX(region) );
}
#endif  // GUARDADDFABR 

#ifndef GUARDMULTIPLYFABR 
#define GUARDMULTIPLYFABR 
// Prototype for Fortran procedure MULTIPLYFABR ...
//
void FORTRAN_NAME( MULTIPLYFABR ,multiplyfabr )(
      CHFp_FRA(dst)
      ,CHFp_CONST_REAL(src)
      ,CHFp_BOX(region) );

#define FORT_MULTIPLYFABR FORTRAN_NAME( inlineMULTIPLYFABR, inlineMULTIPLYFABR)
#define FORTNT_MULTIPLYFABR FORTRAN_NAME( MULTIPLYFABR, multiplyfabr)

inline void FORTRAN_NAME(inlineMULTIPLYFABR, inlineMULTIPLYFABR)(
      CHFp_FRA(dst)
      ,CHFp_CONST_REAL(src)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_MULTIPLYFABR");
 FORTRAN_NAME( MULTIPLYFABR ,multiplyfabr )(
      CHFt_FRA(dst)
      ,CHFt_CONST_REAL(src)
      ,CHFt_BOX(region) );
}
#endif  // GUARDMULTIPLYFABR 

#ifndef GUARDDIVIDEFABR 
#define GUARDDIVIDEFABR 
// Prototype for Fortran procedure DIVIDEFABR ...
//
void FORTRAN_NAME( DIVIDEFABR ,dividefabr )(
      CHFp_FRA(dst)
      ,CHFp_CONST_REAL(src)
      ,CHFp_BOX(region) );

#define FORT_DIVIDEFABR FORTRAN_NAME( inlineDIVIDEFABR, inlineDIVIDEFABR)
#define FORTNT_DIVIDEFABR FORTRAN_NAME( DIVIDEFABR, dividefabr)

inline void FORTRAN_NAME(inlineDIVIDEFABR, inlineDIVIDEFABR)(
      CHFp_FRA(dst)
      ,CHFp_CONST_REAL(src)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_DIVIDEFABR");
 FORTRAN_NAME( DIVIDEFABR ,dividefabr )(
      CHFt_FRA(dst)
      ,CHFt_CONST_REAL(src)
      ,CHFt_BOX(region) );
}
#endif  // GUARDDIVIDEFABR 

#ifndef GUARDEBDOTPRODUCT 
#define GUARDEBDOTPRODUCT 
// Prototype for Fortran procedure EBDOTPRODUCT ...
//
void FORTRAN_NAME( EBDOTPRODUCT ,ebdotproduct )(
      CHFp_REAL(dotprodout)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_CONST_FRA(bfab)
      ,CHFp_CONST_FRA(volfrac)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(icomp) );

#define FORT_EBDOTPRODUCT FORTRAN_NAME( inlineEBDOTPRODUCT, inlineEBDOTPRODUCT)
#define FORTNT_EBDOTPRODUCT FORTRAN_NAME( EBDOTPRODUCT, ebdotproduct)

inline void FORTRAN_NAME(inlineEBDOTPRODUCT, inlineEBDOTPRODUCT)(
      CHFp_REAL(dotprodout)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_CONST_FRA(bfab)
      ,CHFp_CONST_FRA(volfrac)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(icomp) )
{
 CH_TIMELEAF("FORT_EBDOTPRODUCT");
 FORTRAN_NAME( EBDOTPRODUCT ,ebdotproduct )(
      CHFt_REAL(dotprodout)
      ,CHFt_CONST_FRA(afab)
      ,CHFt_CONST_FRA(bfab)
      ,CHFt_CONST_FRA(volfrac)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(icomp) );
}
#endif  // GUARDEBDOTPRODUCT 

#ifndef GUARDMAXFAB 
#define GUARDMAXFAB 
// Prototype for Fortran procedure MAXFAB ...
//
void FORTRAN_NAME( MAXFAB ,maxfab )(
      CHFp_REAL(aval)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(acomp) );

#define FORT_MAXFAB FORTRAN_NAME( inlineMAXFAB, inlineMAXFAB)
#define FORTNT_MAXFAB FORTRAN_NAME( MAXFAB, maxfab)

inline void FORTRAN_NAME(inlineMAXFAB, inlineMAXFAB)(
      CHFp_REAL(aval)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(acomp) )
{
 CH_TIMELEAF("FORT_MAXFAB");
 FORTRAN_NAME( MAXFAB ,maxfab )(
      CHFt_REAL(aval)
      ,CHFt_CONST_FRA(afab)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(acomp) );
}
#endif  // GUARDMAXFAB 

#ifndef GUARDMINFAB 
#define GUARDMINFAB 
// Prototype for Fortran procedure MINFAB ...
//
void FORTRAN_NAME( MINFAB ,minfab )(
      CHFp_REAL(aval)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(acomp) );

#define FORT_MINFAB FORTRAN_NAME( inlineMINFAB, inlineMINFAB)
#define FORTNT_MINFAB FORTRAN_NAME( MINFAB, minfab)

inline void FORTRAN_NAME(inlineMINFAB, inlineMINFAB)(
      CHFp_REAL(aval)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(acomp) )
{
 CH_TIMELEAF("FORT_MINFAB");
 FORTRAN_NAME( MINFAB ,minfab )(
      CHFt_REAL(aval)
      ,CHFt_CONST_FRA(afab)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(acomp) );
}
#endif  // GUARDMINFAB 

}

#endif
