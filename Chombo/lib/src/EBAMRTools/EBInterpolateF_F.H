#ifndef _EBINTERPOLATEF_F_H_
#define _EBINTERPOLATEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDEBINTERPCONSTANT 
#define GUARDEBINTERPCONSTANT 
// Prototype for Fortran procedure ebinterpconstant ...
//
void FORTRAN_NAME( EBINTERPCONSTANT ,ebinterpconstant )(
      CHFp_FRA1(fine)
      ,CHFp_CONST_FRA1(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(refratio)
      ,CHFp_BOX(bref) );

#define FORT_EBINTERPCONSTANT FORTRAN_NAME( inlineEBINTERPCONSTANT, inlineEBINTERPCONSTANT)
#define FORTNT_EBINTERPCONSTANT FORTRAN_NAME( EBINTERPCONSTANT, ebinterpconstant)

inline void FORTRAN_NAME(inlineEBINTERPCONSTANT, inlineEBINTERPCONSTANT)(
      CHFp_FRA1(fine)
      ,CHFp_CONST_FRA1(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(refratio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_EBINTERPCONSTANT");
 FORTRAN_NAME( EBINTERPCONSTANT ,ebinterpconstant )(
      CHFt_FRA1(fine)
      ,CHFt_CONST_FRA1(coarse)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(refratio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDEBINTERPCONSTANT 

#ifndef GUARDEBCENTRALSLOPE 
#define GUARDEBCENTRALSLOPE 
// Prototype for Fortran procedure ebcentralslope ...
//
void FORTRAN_NAME( EBCENTRALSLOPE ,ebcentralslope )(
      CHFp_FRA1(slope)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir) );

#define FORT_EBCENTRALSLOPE FORTRAN_NAME( inlineEBCENTRALSLOPE, inlineEBCENTRALSLOPE)
#define FORTNT_EBCENTRALSLOPE FORTRAN_NAME( EBCENTRALSLOPE, ebcentralslope)

inline void FORTRAN_NAME(inlineEBCENTRALSLOPE, inlineEBCENTRALSLOPE)(
      CHFp_FRA1(slope)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_EBCENTRALSLOPE");
 FORTRAN_NAME( EBCENTRALSLOPE ,ebcentralslope )(
      CHFt_FRA1(slope)
      ,CHFt_CONST_FRA1(state)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDEBCENTRALSLOPE 

#ifndef GUARDEBHISIDESLOPE 
#define GUARDEBHISIDESLOPE 
// Prototype for Fortran procedure ebhisideslope ...
//
void FORTRAN_NAME( EBHISIDESLOPE ,ebhisideslope )(
      CHFp_FRA1(slope)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir) );

#define FORT_EBHISIDESLOPE FORTRAN_NAME( inlineEBHISIDESLOPE, inlineEBHISIDESLOPE)
#define FORTNT_EBHISIDESLOPE FORTRAN_NAME( EBHISIDESLOPE, ebhisideslope)

inline void FORTRAN_NAME(inlineEBHISIDESLOPE, inlineEBHISIDESLOPE)(
      CHFp_FRA1(slope)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_EBHISIDESLOPE");
 FORTRAN_NAME( EBHISIDESLOPE ,ebhisideslope )(
      CHFt_FRA1(slope)
      ,CHFt_CONST_FRA1(state)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDEBHISIDESLOPE 

#ifndef GUARDEBLOSIDESLOPE 
#define GUARDEBLOSIDESLOPE 
// Prototype for Fortran procedure eblosideslope ...
//
void FORTRAN_NAME( EBLOSIDESLOPE ,eblosideslope )(
      CHFp_FRA1(slope)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir) );

#define FORT_EBLOSIDESLOPE FORTRAN_NAME( inlineEBLOSIDESLOPE, inlineEBLOSIDESLOPE)
#define FORTNT_EBLOSIDESLOPE FORTRAN_NAME( EBLOSIDESLOPE, eblosideslope)

inline void FORTRAN_NAME(inlineEBLOSIDESLOPE, inlineEBLOSIDESLOPE)(
      CHFp_FRA1(slope)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_EBLOSIDESLOPE");
 FORTRAN_NAME( EBLOSIDESLOPE ,eblosideslope )(
      CHFt_FRA1(slope)
      ,CHFt_CONST_FRA1(state)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDEBLOSIDESLOPE 

#ifndef GUARDEBMAXMINMOD 
#define GUARDEBMAXMINMOD 
// Prototype for Fortran procedure ebmaxminmod ...
//
void FORTRAN_NAME( EBMAXMINMOD ,ebmaxminmod )(
      CHFp_FRA1(mmslope)
      ,CHFp_CONST_FRA1(loslope)
      ,CHFp_CONST_FRA1(hislope)
      ,CHFp_BOX(slopebox) );

#define FORT_EBMAXMINMOD FORTRAN_NAME( inlineEBMAXMINMOD, inlineEBMAXMINMOD)
#define FORTNT_EBMAXMINMOD FORTRAN_NAME( EBMAXMINMOD, ebmaxminmod)

inline void FORTRAN_NAME(inlineEBMAXMINMOD, inlineEBMAXMINMOD)(
      CHFp_FRA1(mmslope)
      ,CHFp_CONST_FRA1(loslope)
      ,CHFp_CONST_FRA1(hislope)
      ,CHFp_BOX(slopebox) )
{
 CH_TIMELEAF("FORT_EBMAXMINMOD");
 FORTRAN_NAME( EBMAXMINMOD ,ebmaxminmod )(
      CHFt_FRA1(mmslope)
      ,CHFt_CONST_FRA1(loslope)
      ,CHFt_CONST_FRA1(hislope)
      ,CHFt_BOX(slopebox) );
}
#endif  // GUARDEBMAXMINMOD 

#ifndef GUARDEBINTERPLINEAR 
#define GUARDEBINTERPLINEAR 
// Prototype for Fortran procedure ebinterplinear ...
//
void FORTRAN_NAME( EBINTERPLINEAR ,ebinterplinear )(
      CHFp_FRA1(fine)
      ,CHFp_CONST_FRA1(slope)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(refratio)
      ,CHFp_BOX(bref) );

#define FORT_EBINTERPLINEAR FORTRAN_NAME( inlineEBINTERPLINEAR, inlineEBINTERPLINEAR)
#define FORTNT_EBINTERPLINEAR FORTRAN_NAME( EBINTERPLINEAR, ebinterplinear)

inline void FORTRAN_NAME(inlineEBINTERPLINEAR, inlineEBINTERPLINEAR)(
      CHFp_FRA1(fine)
      ,CHFp_CONST_FRA1(slope)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(refratio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_EBINTERPLINEAR");
 FORTRAN_NAME( EBINTERPLINEAR ,ebinterplinear )(
      CHFt_FRA1(fine)
      ,CHFt_CONST_FRA1(slope)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(refratio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDEBINTERPLINEAR 

#ifndef GUARDEBINTERPSMOOTHERLINEAR 
#define GUARDEBINTERPSMOOTHERLINEAR 
// Prototype for Fortran procedure ebinterpsmootherlinear ...
//
void FORTRAN_NAME( EBINTERPSMOOTHERLINEAR ,ebinterpsmootherlinear )(
      CHFp_FRA1(fine)
      ,CHFp_CONST_FRA1(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(refratio)
      ,CHFp_BOX(bref) );

#define FORT_EBINTERPSMOOTHERLINEAR FORTRAN_NAME( inlineEBINTERPSMOOTHERLINEAR, inlineEBINTERPSMOOTHERLINEAR)
#define FORTNT_EBINTERPSMOOTHERLINEAR FORTRAN_NAME( EBINTERPSMOOTHERLINEAR, ebinterpsmootherlinear)

inline void FORTRAN_NAME(inlineEBINTERPSMOOTHERLINEAR, inlineEBINTERPSMOOTHERLINEAR)(
      CHFp_FRA1(fine)
      ,CHFp_CONST_FRA1(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(refratio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_EBINTERPSMOOTHERLINEAR");
 FORTRAN_NAME( EBINTERPSMOOTHERLINEAR ,ebinterpsmootherlinear )(
      CHFt_FRA1(fine)
      ,CHFt_CONST_FRA1(coarse)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(refratio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDEBINTERPSMOOTHERLINEAR 

#ifndef GUARDEBINTERPQUADRATIC 
#define GUARDEBINTERPQUADRATIC 
// Prototype for Fortran procedure ebinterpquadratic ...
//
void FORTRAN_NAME( EBINTERPQUADRATIC ,ebinterpquadratic )(
      CHFp_FRA1(fine)
      ,CHFp_CONST_FRA1(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(refratio)
      ,CHFp_BOX(bref) );

#define FORT_EBINTERPQUADRATIC FORTRAN_NAME( inlineEBINTERPQUADRATIC, inlineEBINTERPQUADRATIC)
#define FORTNT_EBINTERPQUADRATIC FORTRAN_NAME( EBINTERPQUADRATIC, ebinterpquadratic)

inline void FORTRAN_NAME(inlineEBINTERPQUADRATIC, inlineEBINTERPQUADRATIC)(
      CHFp_FRA1(fine)
      ,CHFp_CONST_FRA1(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(refratio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_EBINTERPQUADRATIC");
 FORTRAN_NAME( EBINTERPQUADRATIC ,ebinterpquadratic )(
      CHFt_FRA1(fine)
      ,CHFt_CONST_FRA1(coarse)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(refratio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDEBINTERPQUADRATIC 

#ifndef GUARDEBINTERPQUADSHIFT 
#define GUARDEBINTERPQUADSHIFT 
// Prototype for Fortran procedure ebinterpquadshift ...
//
void FORTRAN_NAME( EBINTERPQUADSHIFT ,ebinterpquadshift )(
      CHFp_FRA1(fine)
      ,CHFp_CONST_FRA1(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(refratio)
      ,CHFp_CONST_INTVECT(ishift)
      ,CHFp_BOX(bref) );

#define FORT_EBINTERPQUADSHIFT FORTRAN_NAME( inlineEBINTERPQUADSHIFT, inlineEBINTERPQUADSHIFT)
#define FORTNT_EBINTERPQUADSHIFT FORTRAN_NAME( EBINTERPQUADSHIFT, ebinterpquadshift)

inline void FORTRAN_NAME(inlineEBINTERPQUADSHIFT, inlineEBINTERPQUADSHIFT)(
      CHFp_FRA1(fine)
      ,CHFp_CONST_FRA1(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(refratio)
      ,CHFp_CONST_INTVECT(ishift)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_EBINTERPQUADSHIFT");
 FORTRAN_NAME( EBINTERPQUADSHIFT ,ebinterpquadshift )(
      CHFt_FRA1(fine)
      ,CHFt_CONST_FRA1(coarse)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(refratio)
      ,CHFt_CONST_INTVECT(ishift)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDEBINTERPQUADSHIFT 

}

#endif
