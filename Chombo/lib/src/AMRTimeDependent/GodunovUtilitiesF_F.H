#ifndef _GODUNOVUTILITIESF_F_H_
#define _GODUNOVUTILITIESF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMINFLATF 
#define GUARDMINFLATF 
// Prototype for Fortran procedure MINFLATF ...
//
void FORTRAN_NAME( MINFLATF ,minflatf )(
      CHFp_FRA1(flattening)
      ,CHFp_CONST_FRA(zetaDir)
      ,CHFp_CONST_FRA(du)
      ,CHFp_BOX(box) );

#define FORT_MINFLATF FORTRAN_NAME( inlineMINFLATF, inlineMINFLATF)
#define FORTNT_MINFLATF FORTRAN_NAME( MINFLATF, minflatf)

inline void FORTRAN_NAME(inlineMINFLATF, inlineMINFLATF)(
      CHFp_FRA1(flattening)
      ,CHFp_CONST_FRA(zetaDir)
      ,CHFp_CONST_FRA(du)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINFLATF");
 FORTRAN_NAME( MINFLATF ,minflatf )(
      CHFt_FRA1(flattening)
      ,CHFt_CONST_FRA(zetaDir)
      ,CHFt_CONST_FRA(du)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINFLATF 

#ifndef GUARDGETDPTWOF 
#define GUARDGETDPTWOF 
// Prototype for Fortran procedure GETDPTWOF ...
//
void FORTRAN_NAME( GETDPTWOF ,getdptwof )(
      CHFp_FRA1(delta2p)
      ,CHFp_CONST_FRA1(delta1p)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_GETDPTWOF FORTRAN_NAME( inlineGETDPTWOF, inlineGETDPTWOF)
#define FORTNT_GETDPTWOF FORTRAN_NAME( GETDPTWOF, getdptwof)

inline void FORTRAN_NAME(inlineGETDPTWOF, inlineGETDPTWOF)(
      CHFp_FRA1(delta2p)
      ,CHFp_CONST_FRA1(delta1p)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_GETDPTWOF");
 FORTRAN_NAME( GETDPTWOF ,getdptwof )(
      CHFt_FRA1(delta2p)
      ,CHFt_CONST_FRA1(delta1p)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDGETDPTWOF 

#ifndef GUARDGETFLATF 
#define GUARDGETFLATF 
// Prototype for Fortran procedure GETFLATF ...
//
void FORTRAN_NAME( GETFLATF ,getflatf )(
      CHFp_FRA1(zetaTwiddle)
      ,CHFp_CONST_FRA1(delta1p)
      ,CHFp_CONST_FRA1(delta2p)
      ,CHFp_CONST_REAL(smallp)
      ,CHFp_CONST_FRA1(bulkMin)
      ,CHFp_BOX(box) );

#define FORT_GETFLATF FORTRAN_NAME( inlineGETFLATF, inlineGETFLATF)
#define FORTNT_GETFLATF FORTRAN_NAME( GETFLATF, getflatf)

inline void FORTRAN_NAME(inlineGETFLATF, inlineGETFLATF)(
      CHFp_FRA1(zetaTwiddle)
      ,CHFp_CONST_FRA1(delta1p)
      ,CHFp_CONST_FRA1(delta2p)
      ,CHFp_CONST_REAL(smallp)
      ,CHFp_CONST_FRA1(bulkMin)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETFLATF");
 FORTRAN_NAME( GETFLATF ,getflatf )(
      CHFt_FRA1(zetaTwiddle)
      ,CHFt_CONST_FRA1(delta1p)
      ,CHFt_CONST_FRA1(delta2p)
      ,CHFt_CONST_REAL(smallp)
      ,CHFt_CONST_FRA1(bulkMin)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETFLATF 

#ifndef GUARDGETGRADF 
#define GUARDGETGRADF 
// Prototype for Fortran procedure GETGRADF ...
//
void FORTRAN_NAME( GETGRADF ,getgradf )(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_GETGRADF FORTRAN_NAME( inlineGETGRADF, inlineGETGRADF)
#define FORTNT_GETGRADF FORTRAN_NAME( GETGRADF, getgradf)

inline void FORTRAN_NAME(inlineGETGRADF, inlineGETGRADF)(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_GETGRADF");
 FORTRAN_NAME( GETGRADF ,getgradf )(
      CHFt_FRA1(du)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDGETGRADF 

#ifndef GUARDMIN3PTSF 
#define GUARDMIN3PTSF 
// Prototype for Fortran procedure MIN3PTSF ...
//
void FORTRAN_NAME( MIN3PTSF ,min3ptsf )(
      CHFp_FRA1(mindata)
      ,CHFp_CONST_FRA1(data)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_MIN3PTSF FORTRAN_NAME( inlineMIN3PTSF, inlineMIN3PTSF)
#define FORTNT_MIN3PTSF FORTRAN_NAME( MIN3PTSF, min3ptsf)

inline void FORTRAN_NAME(inlineMIN3PTSF, inlineMIN3PTSF)(
      CHFp_FRA1(mindata)
      ,CHFp_CONST_FRA1(data)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_MIN3PTSF");
 FORTRAN_NAME( MIN3PTSF ,min3ptsf )(
      CHFt_FRA1(mindata)
      ,CHFt_CONST_FRA1(data)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDMIN3PTSF 

#ifndef GUARDSECONDSLOPEDIFFSF 
#define GUARDSECONDSLOPEDIFFSF 
// Prototype for Fortran procedure SECONDSLOPEDIFFSF ...
//
void FORTRAN_NAME( SECONDSLOPEDIFFSF ,secondslopediffsf )(
      CHFp_FRA(deltaWC)
      ,CHFp_FRA(deltaWL)
      ,CHFp_FRA(deltaWR)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_SECONDSLOPEDIFFSF FORTRAN_NAME( inlineSECONDSLOPEDIFFSF, inlineSECONDSLOPEDIFFSF)
#define FORTNT_SECONDSLOPEDIFFSF FORTRAN_NAME( SECONDSLOPEDIFFSF, secondslopediffsf)

inline void FORTRAN_NAME(inlineSECONDSLOPEDIFFSF, inlineSECONDSLOPEDIFFSF)(
      CHFp_FRA(deltaWC)
      ,CHFp_FRA(deltaWL)
      ,CHFp_FRA(deltaWR)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_SECONDSLOPEDIFFSF");
 FORTRAN_NAME( SECONDSLOPEDIFFSF ,secondslopediffsf )(
      CHFt_FRA(deltaWC)
      ,CHFt_FRA(deltaWL)
      ,CHFt_FRA(deltaWR)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDSECONDSLOPEDIFFSF 

#ifndef GUARDFOURTHSLOPEDIFFSF 
#define GUARDFOURTHSLOPEDIFFSF 
// Prototype for Fortran procedure FOURTHSLOPEDIFFSF ...
//
void FORTRAN_NAME( FOURTHSLOPEDIFFSF ,fourthslopediffsf )(
      CHFp_FRA(delta4WC)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(delta2W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_FOURTHSLOPEDIFFSF FORTRAN_NAME( inlineFOURTHSLOPEDIFFSF, inlineFOURTHSLOPEDIFFSF)
#define FORTNT_FOURTHSLOPEDIFFSF FORTRAN_NAME( FOURTHSLOPEDIFFSF, fourthslopediffsf)

inline void FORTRAN_NAME(inlineFOURTHSLOPEDIFFSF, inlineFOURTHSLOPEDIFFSF)(
      CHFp_FRA(delta4WC)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(delta2W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_FOURTHSLOPEDIFFSF");
 FORTRAN_NAME( FOURTHSLOPEDIFFSF ,fourthslopediffsf )(
      CHFt_FRA(delta4WC)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(delta2W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDFOURTHSLOPEDIFFSF 

#ifndef GUARDAPPLYFLATF 
#define GUARDAPPLYFLATF 
// Prototype for Fortran procedure APPLYFLATF ...
//
void FORTRAN_NAME( APPLYFLATF ,applyflatf )(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA1(flattening)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_BOX(box) );

#define FORT_APPLYFLATF FORTRAN_NAME( inlineAPPLYFLATF, inlineAPPLYFLATF)
#define FORTNT_APPLYFLATF FORTRAN_NAME( APPLYFLATF, applyflatf)

inline void FORTRAN_NAME(inlineAPPLYFLATF, inlineAPPLYFLATF)(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA1(flattening)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_APPLYFLATF");
 FORTRAN_NAME( APPLYFLATF ,applyflatf )(
      CHFt_FRA(dW)
      ,CHFt_CONST_FRA1(flattening)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_BOX(box) );
}
#endif  // GUARDAPPLYFLATF 

#ifndef GUARDVANLEERLIMITERF 
#define GUARDVANLEERLIMITERF 
// Prototype for Fortran procedure VANLEERLIMITERF ...
//
void FORTRAN_NAME( VANLEERLIMITERF ,vanleerlimiterf )(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(dWleft)
      ,CHFp_CONST_FRA(dWright)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_BOX(box) );

#define FORT_VANLEERLIMITERF FORTRAN_NAME( inlineVANLEERLIMITERF, inlineVANLEERLIMITERF)
#define FORTNT_VANLEERLIMITERF FORTRAN_NAME( VANLEERLIMITERF, vanleerlimiterf)

inline void FORTRAN_NAME(inlineVANLEERLIMITERF, inlineVANLEERLIMITERF)(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(dWleft)
      ,CHFp_CONST_FRA(dWright)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VANLEERLIMITERF");
 FORTRAN_NAME( VANLEERLIMITERF ,vanleerlimiterf )(
      CHFt_FRA(dW)
      ,CHFt_CONST_FRA(dWleft)
      ,CHFt_CONST_FRA(dWright)
      ,CHFt_CONST_INT(numslopes)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVANLEERLIMITERF 

#ifndef GUARDEXTPRESERVINGVANLEERLIMITERF 
#define GUARDEXTPRESERVINGVANLEERLIMITERF 
// Prototype for Fortran procedure EXTPRESERVINGVANLEERLIMITERF ...
//
void FORTRAN_NAME( EXTPRESERVINGVANLEERLIMITERF ,extpreservingvanleerlimiterf )(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(dWleft)
      ,CHFp_CONST_FRA(dWright)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox) );

#define FORT_EXTPRESERVINGVANLEERLIMITERF FORTRAN_NAME( inlineEXTPRESERVINGVANLEERLIMITERF, inlineEXTPRESERVINGVANLEERLIMITERF)
#define FORTNT_EXTPRESERVINGVANLEERLIMITERF FORTRAN_NAME( EXTPRESERVINGVANLEERLIMITERF, extpreservingvanleerlimiterf)

inline void FORTRAN_NAME(inlineEXTPRESERVINGVANLEERLIMITERF, inlineEXTPRESERVINGVANLEERLIMITERF)(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(dWleft)
      ,CHFp_CONST_FRA(dWright)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_EXTPRESERVINGVANLEERLIMITERF");
 FORTRAN_NAME( EXTPRESERVINGVANLEERLIMITERF ,extpreservingvanleerlimiterf )(
      CHFt_FRA(dW)
      ,CHFt_CONST_FRA(dWleft)
      ,CHFt_CONST_FRA(dWright)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDEXTPRESERVINGVANLEERLIMITERF 

#ifndef GUARDPLMNORMALPREDF 
#define GUARDPLMNORMALPREDF 
// Prototype for Fortran procedure PLMNORMALPREDF ...
//
void FORTRAN_NAME( PLMNORMALPREDF ,plmnormalpredf )(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(dW)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(nSlope)
      ,CHFp_BOX(box) );

#define FORT_PLMNORMALPREDF FORTRAN_NAME( inlinePLMNORMALPREDF, inlinePLMNORMALPREDF)
#define FORTNT_PLMNORMALPREDF FORTRAN_NAME( PLMNORMALPREDF, plmnormalpredf)

inline void FORTRAN_NAME(inlinePLMNORMALPREDF, inlinePLMNORMALPREDF)(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(dW)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(nSlope)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PLMNORMALPREDF");
 FORTRAN_NAME( PLMNORMALPREDF ,plmnormalpredf )(
      CHFt_FRA(dWMinus)
      ,CHFt_FRA(dWPlus)
      ,CHFt_CONST_FRA(dW)
      ,CHFt_CONST_FRA(lambda)
      ,CHFt_CONST_REAL(dtbydx)
      ,CHFt_CONST_INT(nSlope)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPLMNORMALPREDF 

#ifndef GUARDDIVUEDGEF 
#define GUARDDIVUEDGEF 
// Prototype for Fortran procedure DIVUEDGEF ...
//
void FORTRAN_NAME( DIVUEDGEF ,divuedgef )(
      CHFp_FRA1(divu)
      ,CHFp_CONST_FRA1(uNorm)
      ,CHFp_CONST_FRA(duTan)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_DIVUEDGEF FORTRAN_NAME( inlineDIVUEDGEF, inlineDIVUEDGEF)
#define FORTNT_DIVUEDGEF FORTRAN_NAME( DIVUEDGEF, divuedgef)

inline void FORTRAN_NAME(inlineDIVUEDGEF, inlineDIVUEDGEF)(
      CHFp_FRA1(divu)
      ,CHFp_CONST_FRA1(uNorm)
      ,CHFp_CONST_FRA(duTan)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_DIVUEDGEF");
 FORTRAN_NAME( DIVUEDGEF ,divuedgef )(
      CHFt_FRA1(divu)
      ,CHFt_CONST_FRA1(uNorm)
      ,CHFt_CONST_FRA(duTan)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDDIVUEDGEF 

#ifndef GUARDARTVISCF 
#define GUARDARTVISCF 
// Prototype for Fortran procedure ARTVISCF ...
//
void FORTRAN_NAME( ARTVISCF ,artviscf )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_FRA1(divu)
      ,CHFp_CONST_REAL(coeff)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(fbox) );

#define FORT_ARTVISCF FORTRAN_NAME( inlineARTVISCF, inlineARTVISCF)
#define FORTNT_ARTVISCF FORTRAN_NAME( ARTVISCF, artviscf)

inline void FORTRAN_NAME(inlineARTVISCF, inlineARTVISCF)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_FRA1(divu)
      ,CHFp_CONST_REAL(coeff)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(fbox) )
{
 CH_TIMELEAF("FORT_ARTVISCF");
 FORTRAN_NAME( ARTVISCF ,artviscf )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_FRA1(divu)
      ,CHFt_CONST_REAL(coeff)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(fbox) );
}
#endif  // GUARDARTVISCF 

#ifndef GUARDPPMLIMITERF 
#define GUARDPPMLIMITERF 
// Prototype for Fortran procedure PPMLIMITERF ...
//
void FORTRAN_NAME( PPMLIMITERF ,ppmlimiterf )(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_BOX(box) );

#define FORT_PPMLIMITERF FORTRAN_NAME( inlinePPMLIMITERF, inlinePPMLIMITERF)
#define FORTNT_PPMLIMITERF FORTRAN_NAME( PPMLIMITERF, ppmlimiterf)

inline void FORTRAN_NAME(inlinePPMLIMITERF, inlinePPMLIMITERF)(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PPMLIMITERF");
 FORTRAN_NAME( PPMLIMITERF ,ppmlimiterf )(
      CHFt_FRA(dWMinus)
      ,CHFt_FRA(dWPlus)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPPMLIMITERF 

#ifndef GUARDCOLELLASEKORALIMITERF 
#define GUARDCOLELLASEKORALIMITERF 
// Prototype for Fortran procedure COLELLASEKORALIMITERF ...
//
void FORTRAN_NAME( COLELLASEKORALIMITERF ,colellasekoralimiterf )(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(diffW)
      ,CHFp_CONST_FRA(d2W)
      ,CHFp_CONST_FRA(dW2fcf)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox)
      ,CHFp_CONST_REAL(limitC)
      ,CHFp_CONST_REAL(eps) );

#define FORT_COLELLASEKORALIMITERF FORTRAN_NAME( inlineCOLELLASEKORALIMITERF, inlineCOLELLASEKORALIMITERF)
#define FORTNT_COLELLASEKORALIMITERF FORTRAN_NAME( COLELLASEKORALIMITERF, colellasekoralimiterf)

inline void FORTRAN_NAME(inlineCOLELLASEKORALIMITERF, inlineCOLELLASEKORALIMITERF)(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(diffW)
      ,CHFp_CONST_FRA(d2W)
      ,CHFp_CONST_FRA(dW2fcf)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox)
      ,CHFp_CONST_REAL(limitC)
      ,CHFp_CONST_REAL(eps) )
{
 CH_TIMELEAF("FORT_COLELLASEKORALIMITERF");
 FORTRAN_NAME( COLELLASEKORALIMITERF ,colellasekoralimiterf )(
      CHFt_FRA(dWMinus)
      ,CHFt_FRA(dWPlus)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(diffW)
      ,CHFt_CONST_FRA(d2W)
      ,CHFt_CONST_FRA(dW2fcf)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox)
      ,CHFt_CONST_REAL(limitC)
      ,CHFt_CONST_REAL(eps) );
}
#endif  // GUARDCOLELLASEKORALIMITERF 

#ifndef GUARDPPMFACEVALUESF 
#define GUARDPPMFACEVALUESF 
// Prototype for Fortran procedure PPMFACEVALUESF ...
//
void FORTRAN_NAME( PPMFACEVALUESF ,ppmfacevaluesf )(
      CHFp_FRA(WFace)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(dW)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_PPMFACEVALUESF FORTRAN_NAME( inlinePPMFACEVALUESF, inlinePPMFACEVALUESF)
#define FORTNT_PPMFACEVALUESF FORTRAN_NAME( PPMFACEVALUESF, ppmfacevaluesf)

inline void FORTRAN_NAME(inlinePPMFACEVALUESF, inlinePPMFACEVALUESF)(
      CHFp_FRA(WFace)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(dW)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_PPMFACEVALUESF");
 FORTRAN_NAME( PPMFACEVALUESF ,ppmfacevaluesf )(
      CHFt_FRA(WFace)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(dW)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDPPMFACEVALUESF 

#ifndef GUARDPPMNORMALPREDF 
#define GUARDPPMNORMALPREDF 
// Prototype for Fortran procedure PPMNORMALPREDF ...
//
void FORTRAN_NAME( PPMNORMALPREDF ,ppmnormalpredf )(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(nSlope)
      ,CHFp_BOX(box) );

#define FORT_PPMNORMALPREDF FORTRAN_NAME( inlinePPMNORMALPREDF, inlinePPMNORMALPREDF)
#define FORTNT_PPMNORMALPREDF FORTRAN_NAME( PPMNORMALPREDF, ppmnormalpredf)

inline void FORTRAN_NAME(inlinePPMNORMALPREDF, inlinePPMNORMALPREDF)(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(nSlope)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PPMNORMALPREDF");
 FORTRAN_NAME( PPMNORMALPREDF ,ppmnormalpredf )(
      CHFt_FRA(dWMinus)
      ,CHFt_FRA(dWPlus)
      ,CHFt_CONST_FRA(lambda)
      ,CHFt_CONST_REAL(dtbydx)
      ,CHFt_CONST_INT(nSlope)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPPMNORMALPREDF 

#ifndef GUARDGETSECONDDIFF 
#define GUARDGETSECONDDIFF 
// Prototype for Fortran procedure GETSECONDDIFF ...
//
void FORTRAN_NAME( GETSECONDDIFF ,getseconddiff )(
      CHFp_FRA(d2W)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_GETSECONDDIFF FORTRAN_NAME( inlineGETSECONDDIFF, inlineGETSECONDDIFF)
#define FORTNT_GETSECONDDIFF FORTRAN_NAME( GETSECONDDIFF, getseconddiff)

inline void FORTRAN_NAME(inlineGETSECONDDIFF, inlineGETSECONDDIFF)(
      CHFp_FRA(d2W)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_GETSECONDDIFF");
 FORTRAN_NAME( GETSECONDDIFF ,getseconddiff )(
      CHFt_FRA(d2W)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDGETSECONDDIFF 

#ifndef GUARDGETD2CELL 
#define GUARDGETD2CELL 
// Prototype for Fortran procedure GETD2CELL ...
//
void FORTRAN_NAME( GETD2CELL ,getd2cell )(
      CHFp_FRA(d2W)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_GETD2CELL FORTRAN_NAME( inlineGETD2CELL, inlineGETD2CELL)
#define FORTNT_GETD2CELL FORTRAN_NAME( GETD2CELL, getd2cell)

inline void FORTRAN_NAME(inlineGETD2CELL, inlineGETD2CELL)(
      CHFp_FRA(d2W)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_GETD2CELL");
 FORTRAN_NAME( GETD2CELL ,getd2cell )(
      CHFt_FRA(d2W)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDGETD2CELL 

#ifndef GUARDGETD2LIMFACE 
#define GUARDGETD2LIMFACE 
// Prototype for Fortran procedure GETD2LIMFACE ...
//
void FORTRAN_NAME( GETD2LIMFACE ,getd2limface )(
      CHFp_FRA(d2Wlim)
      ,CHFp_CONST_FRA(d2Wc)
      ,CHFp_CONST_FRA(d2Wcfc)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox)
      ,CHFp_CONST_REAL(limitC) );

#define FORT_GETD2LIMFACE FORTRAN_NAME( inlineGETD2LIMFACE, inlineGETD2LIMFACE)
#define FORTNT_GETD2LIMFACE FORTRAN_NAME( GETD2LIMFACE, getd2limface)

inline void FORTRAN_NAME(inlineGETD2LIMFACE, inlineGETD2LIMFACE)(
      CHFp_FRA(d2Wlim)
      ,CHFp_CONST_FRA(d2Wc)
      ,CHFp_CONST_FRA(d2Wcfc)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox)
      ,CHFp_CONST_REAL(limitC) )
{
 CH_TIMELEAF("FORT_GETD2LIMFACE");
 FORTRAN_NAME( GETD2LIMFACE ,getd2limface )(
      CHFt_FRA(d2Wlim)
      ,CHFt_CONST_FRA(d2Wc)
      ,CHFt_CONST_FRA(d2Wcfc)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(centerBox)
      ,CHFt_CONST_REAL(limitC) );
}
#endif  // GUARDGETD2LIMFACE 

#ifndef GUARDCOLELLASEKORAFACELIMITER 
#define GUARDCOLELLASEKORAFACELIMITER 
// Prototype for Fortran procedure COLELLASEKORAFACELIMITER ...
//
void FORTRAN_NAME( COLELLASEKORAFACELIMITER ,colellasekorafacelimiter )(
      CHFp_FRA(Wface)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(d2Wc)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox)
      ,CHFp_CONST_REAL(limitC) );

#define FORT_COLELLASEKORAFACELIMITER FORTRAN_NAME( inlineCOLELLASEKORAFACELIMITER, inlineCOLELLASEKORAFACELIMITER)
#define FORTNT_COLELLASEKORAFACELIMITER FORTRAN_NAME( COLELLASEKORAFACELIMITER, colellasekorafacelimiter)

inline void FORTRAN_NAME(inlineCOLELLASEKORAFACELIMITER, inlineCOLELLASEKORAFACELIMITER)(
      CHFp_FRA(Wface)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(d2Wc)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox)
      ,CHFp_CONST_REAL(limitC) )
{
 CH_TIMELEAF("FORT_COLELLASEKORAFACELIMITER");
 FORTRAN_NAME( COLELLASEKORAFACELIMITER ,colellasekorafacelimiter )(
      CHFt_FRA(Wface)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(d2Wc)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(centerBox)
      ,CHFt_CONST_REAL(limitC) );
}
#endif  // GUARDCOLELLASEKORAFACELIMITER 

#ifndef GUARDFOURTHINTERPFACES 
#define GUARDFOURTHINTERPFACES 
// Prototype for Fortran procedure FOURTHINTERPFACES ...
//
void FORTRAN_NAME( FOURTHINTERPFACES ,fourthinterpfaces )(
      CHFp_FRA(Wface)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_BOX(nextLoBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_BOX(nextHiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_FOURTHINTERPFACES FORTRAN_NAME( inlineFOURTHINTERPFACES, inlineFOURTHINTERPFACES)
#define FORTNT_FOURTHINTERPFACES FORTRAN_NAME( FOURTHINTERPFACES, fourthinterpfaces)

inline void FORTRAN_NAME(inlineFOURTHINTERPFACES, inlineFOURTHINTERPFACES)(
      CHFp_FRA(Wface)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_BOX(nextLoBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_BOX(nextHiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_FOURTHINTERPFACES");
 FORTRAN_NAME( FOURTHINTERPFACES ,fourthinterpfaces )(
      CHFt_FRA(Wface)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_BOX(nextLoBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_BOX(nextHiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDFOURTHINTERPFACES 

#ifndef GUARDSECONDINTERPFACES 
#define GUARDSECONDINTERPFACES 
// Prototype for Fortran procedure SECONDINTERPFACES ...
//
void FORTRAN_NAME( SECONDINTERPFACES ,secondinterpfaces )(
      CHFp_FRA(Wface)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_SECONDINTERPFACES FORTRAN_NAME( inlineSECONDINTERPFACES, inlineSECONDINTERPFACES)
#define FORTNT_SECONDINTERPFACES FORTRAN_NAME( SECONDINTERPFACES, secondinterpfaces)

inline void FORTRAN_NAME(inlineSECONDINTERPFACES, inlineSECONDINTERPFACES)(
      CHFp_FRA(Wface)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_SECONDINTERPFACES");
 FORTRAN_NAME( SECONDINTERPFACES ,secondinterpfaces )(
      CHFt_FRA(Wface)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDSECONDINTERPFACES 

#ifndef GUARDPPMFOURTHFACE 
#define GUARDPPMFOURTHFACE 
// Prototype for Fortran procedure PPMFOURTHFACE ...
//
void FORTRAN_NAME( PPMFOURTHFACE ,ppmfourthface )(
      CHFp_FRA(Wface)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(d2W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_PPMFOURTHFACE FORTRAN_NAME( inlinePPMFOURTHFACE, inlinePPMFOURTHFACE)
#define FORTNT_PPMFOURTHFACE FORTRAN_NAME( PPMFOURTHFACE, ppmfourthface)

inline void FORTRAN_NAME(inlinePPMFOURTHFACE, inlinePPMFOURTHFACE)(
      CHFp_FRA(Wface)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(d2W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_PPMFOURTHFACE");
 FORTRAN_NAME( PPMFOURTHFACE ,ppmfourthface )(
      CHFt_FRA(Wface)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(d2W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDPPMFOURTHFACE 

#ifndef GUARDHODIVCOEF 
#define GUARDHODIVCOEF 
// Prototype for Fortran procedure HODIVCOEF ...
//
void FORTRAN_NAME( HODIVCOEF ,hodivcoef )(
      CHFp_FRA1(divVel)
      ,CHFp_CONST_FRA1(csq)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(M0sq)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_HODIVCOEF FORTRAN_NAME( inlineHODIVCOEF, inlineHODIVCOEF)
#define FORTNT_HODIVCOEF FORTRAN_NAME( HODIVCOEF, hodivcoef)

inline void FORTRAN_NAME(inlineHODIVCOEF, inlineHODIVCOEF)(
      CHFp_FRA1(divVel)
      ,CHFp_CONST_FRA1(csq)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(M0sq)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_HODIVCOEF");
 FORTRAN_NAME( HODIVCOEF ,hodivcoef )(
      CHFt_FRA1(divVel)
      ,CHFt_CONST_FRA1(csq)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(M0sq)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDHODIVCOEF 

#ifndef GUARDHIGHORDERDIVCO 
#define GUARDHIGHORDERDIVCO 
// Prototype for Fortran procedure HIGHORDERDIVCO ...
//
void FORTRAN_NAME( HIGHORDERDIVCO ,highorderdivco )(
      CHFp_FRA1(divVel)
      ,CHFp_CONST_FRA1(csq)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_REAL(M0sq)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_HIGHORDERDIVCO FORTRAN_NAME( inlineHIGHORDERDIVCO, inlineHIGHORDERDIVCO)
#define FORTNT_HIGHORDERDIVCO FORTRAN_NAME( HIGHORDERDIVCO, highorderdivco)

inline void FORTRAN_NAME(inlineHIGHORDERDIVCO, inlineHIGHORDERDIVCO)(
      CHFp_FRA1(divVel)
      ,CHFp_CONST_FRA1(csq)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_REAL(M0sq)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_HIGHORDERDIVCO");
 FORTRAN_NAME( HIGHORDERDIVCO ,highorderdivco )(
      CHFt_FRA1(divVel)
      ,CHFt_CONST_FRA1(csq)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_REAL(M0sq)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDHIGHORDERDIVCO 

#ifndef GUARDCHECKCUBICLIMITERF 
#define GUARDCHECKCUBICLIMITERF 
// Prototype for Fortran procedure CHECKCUBICLIMITERF ...
//
void FORTRAN_NAME( CHECKCUBICLIMITERF ,checkcubiclimiterf )(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(d2W)
      ,CHFp_CONST_FRA(dW2fcf)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_BOX(nextLoBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_BOX(nextHiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(innerCenterBox)
      ,CHFp_CONST_REAL(limitC)
      ,CHFp_CONST_REAL(C3)
      ,CHFp_CONST_REAL(eps) );

#define FORT_CHECKCUBICLIMITERF FORTRAN_NAME( inlineCHECKCUBICLIMITERF, inlineCHECKCUBICLIMITERF)
#define FORTNT_CHECKCUBICLIMITERF FORTRAN_NAME( CHECKCUBICLIMITERF, checkcubiclimiterf)

inline void FORTRAN_NAME(inlineCHECKCUBICLIMITERF, inlineCHECKCUBICLIMITERF)(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(d2W)
      ,CHFp_CONST_FRA(dW2fcf)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_BOX(nextLoBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_BOX(nextHiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(innerCenterBox)
      ,CHFp_CONST_REAL(limitC)
      ,CHFp_CONST_REAL(C3)
      ,CHFp_CONST_REAL(eps) )
{
 CH_TIMELEAF("FORT_CHECKCUBICLIMITERF");
 FORTRAN_NAME( CHECKCUBICLIMITERF ,checkcubiclimiterf )(
      CHFt_FRA(dWMinus)
      ,CHFt_FRA(dWPlus)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(d2W)
      ,CHFt_CONST_FRA(dW2fcf)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_BOX(nextLoBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_BOX(nextHiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(innerCenterBox)
      ,CHFt_CONST_REAL(limitC)
      ,CHFt_CONST_REAL(C3)
      ,CHFt_CONST_REAL(eps) );
}
#endif  // GUARDCHECKCUBICLIMITERF 

}

#endif
