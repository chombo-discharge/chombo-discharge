#ifndef _EBAVERAGEF_F_H_
#define _EBAVERAGEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDEBAVERAGE 
#define GUARDEBAVERAGE 
// Prototype for Fortran procedure ebaverage ...
//
void FORTRAN_NAME( EBAVERAGE ,ebaverage )(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_BOX(coarbox)
      ,CHFp_CONST_INT(refrat)
      ,CHFp_BOX(bref) );

#define FORT_EBAVERAGE FORTRAN_NAME( inlineEBAVERAGE, inlineEBAVERAGE)
#define FORTNT_EBAVERAGE FORTRAN_NAME( EBAVERAGE, ebaverage)

inline void FORTRAN_NAME(inlineEBAVERAGE, inlineEBAVERAGE)(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_BOX(coarbox)
      ,CHFp_CONST_INT(refrat)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_EBAVERAGE");
 FORTRAN_NAME( EBAVERAGE ,ebaverage )(
      CHFt_FRA1(coarse)
      ,CHFt_CONST_FRA1(fine)
      ,CHFt_BOX(coarbox)
      ,CHFt_CONST_INT(refrat)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDEBAVERAGE 

#ifndef GUARDEBAVERARZ 
#define GUARDEBAVERARZ 
// Prototype for Fortran procedure ebaverarz ...
//
void FORTRAN_NAME( EBAVERARZ ,ebaverarz )(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_BOX(coarbox)
      ,CHFp_CONST_INT(refrat)
      ,CHFp_BOX(bref)
      ,CHFp_REAL(dxcoar)
      ,CHFp_REAL(dxfine) );

#define FORT_EBAVERARZ FORTRAN_NAME( inlineEBAVERARZ, inlineEBAVERARZ)
#define FORTNT_EBAVERARZ FORTRAN_NAME( EBAVERARZ, ebaverarz)

inline void FORTRAN_NAME(inlineEBAVERARZ, inlineEBAVERARZ)(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_BOX(coarbox)
      ,CHFp_CONST_INT(refrat)
      ,CHFp_BOX(bref)
      ,CHFp_REAL(dxcoar)
      ,CHFp_REAL(dxfine) )
{
 CH_TIMELEAF("FORT_EBAVERARZ");
 FORTRAN_NAME( EBAVERARZ ,ebaverarz )(
      CHFt_FRA1(coarse)
      ,CHFt_CONST_FRA1(fine)
      ,CHFt_BOX(coarbox)
      ,CHFt_CONST_INT(refrat)
      ,CHFt_BOX(bref)
      ,CHFt_REAL(dxcoar)
      ,CHFt_REAL(dxfine) );
}
#endif  // GUARDEBAVERARZ 

#ifndef GUARDEBAVERAGEFACE 
#define GUARDEBAVERAGEFACE 
// Prototype for Fortran procedure ebaverageface ...
//
void FORTRAN_NAME( EBAVERAGEFACE ,ebaverageface )(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_BOX(coarbox)
      ,CHFp_CONST_INT(refrat)
      ,CHFp_BOX(bref)
      ,CHFp_CONST_INT(idir) );

#define FORT_EBAVERAGEFACE FORTRAN_NAME( inlineEBAVERAGEFACE, inlineEBAVERAGEFACE)
#define FORTNT_EBAVERAGEFACE FORTRAN_NAME( EBAVERAGEFACE, ebaverageface)

inline void FORTRAN_NAME(inlineEBAVERAGEFACE, inlineEBAVERAGEFACE)(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_BOX(coarbox)
      ,CHFp_CONST_INT(refrat)
      ,CHFp_BOX(bref)
      ,CHFp_CONST_INT(idir) )
{
 CH_TIMELEAF("FORT_EBAVERAGEFACE");
 FORTRAN_NAME( EBAVERAGEFACE ,ebaverageface )(
      CHFt_FRA1(coarse)
      ,CHFt_CONST_FRA1(fine)
      ,CHFt_BOX(coarbox)
      ,CHFt_CONST_INT(refrat)
      ,CHFt_BOX(bref)
      ,CHFt_CONST_INT(idir) );
}
#endif  // GUARDEBAVERAGEFACE 

#ifndef GUARDEBCOARSEN 
#define GUARDEBCOARSEN 
// Prototype for Fortran procedure ebcoarsen ...
//
void FORTRAN_NAME( EBCOARSEN ,ebcoarsen )(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_CONST_FRA1(h2laplfine)
      ,CHFp_BOX(coarbox)
      ,CHFp_CONST_INT(refrat)
      ,CHFp_BOX(bref) );

#define FORT_EBCOARSEN FORTRAN_NAME( inlineEBCOARSEN, inlineEBCOARSEN)
#define FORTNT_EBCOARSEN FORTRAN_NAME( EBCOARSEN, ebcoarsen)

inline void FORTRAN_NAME(inlineEBCOARSEN, inlineEBCOARSEN)(
      CHFp_FRA1(coarse)
      ,CHFp_CONST_FRA1(fine)
      ,CHFp_CONST_FRA1(h2laplfine)
      ,CHFp_BOX(coarbox)
      ,CHFp_CONST_INT(refrat)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_EBCOARSEN");
 FORTRAN_NAME( EBCOARSEN ,ebcoarsen )(
      CHFt_FRA1(coarse)
      ,CHFt_CONST_FRA1(fine)
      ,CHFt_CONST_FRA1(h2laplfine)
      ,CHFt_BOX(coarbox)
      ,CHFt_CONST_INT(refrat)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDEBCOARSEN 

#ifndef GUARDH2LAPL1DADDITIVE 
#define GUARDH2LAPL1DADDITIVE 
// Prototype for Fortran procedure h2lapl1dadditive ...
//
void FORTRAN_NAME( H2LAPL1DADDITIVE ,h2lapl1dadditive )(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_H2LAPL1DADDITIVE FORTRAN_NAME( inlineH2LAPL1DADDITIVE, inlineH2LAPL1DADDITIVE)
#define FORTNT_H2LAPL1DADDITIVE FORTRAN_NAME( H2LAPL1DADDITIVE, h2lapl1dadditive)

inline void FORTRAN_NAME(inlineH2LAPL1DADDITIVE, inlineH2LAPL1DADDITIVE)(
      CHFp_FRA1(opphidir)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_H2LAPL1DADDITIVE");
 FORTRAN_NAME( H2LAPL1DADDITIVE ,h2lapl1dadditive )(
      CHFt_FRA1(opphidir)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDH2LAPL1DADDITIVE 

}

#endif
