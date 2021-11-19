#ifndef _CD_EBCOARSEFINEPARTICLEMESHF_F_H_
#define _CD_EBCOARSEFINEPARTICLEMESHF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDAVERAGE_GHOSTS 
#define GUARDAVERAGE_GHOSTS 
// Prototype for Fortran procedure average_ghosts ...
//
void FORTRAN_NAME( AVERAGE_GHOSTS ,average_ghosts )(
      CHFp_FRA1(dst)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_BOX(finebox)
      ,CHFp_CONST_INT(refrat)
      ,CHFp_CONST_REAL(factor) );

#define FORT_AVERAGE_GHOSTS FORTRAN_NAME( inlineAVERAGE_GHOSTS, inlineAVERAGE_GHOSTS)
#define FORTNT_AVERAGE_GHOSTS FORTRAN_NAME( AVERAGE_GHOSTS, average_ghosts)

inline void FORTRAN_NAME(inlineAVERAGE_GHOSTS, inlineAVERAGE_GHOSTS)(
      CHFp_FRA1(dst)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_BOX(finebox)
      ,CHFp_CONST_INT(refrat)
      ,CHFp_CONST_REAL(factor) )
{
 CH_TIMELEAF("FORT_AVERAGE_GHOSTS");
 FORTRAN_NAME( AVERAGE_GHOSTS ,average_ghosts )(
      CHFt_FRA1(dst)
      ,CHFt_CONST_FRA1(src)
      ,CHFt_BOX(finebox)
      ,CHFt_CONST_INT(refrat)
      ,CHFt_CONST_REAL(factor) );
}
#endif  // GUARDAVERAGE_GHOSTS 

#ifndef GUARDPIECEWISE_CONSTANT_INTERP 
#define GUARDPIECEWISE_CONSTANT_INTERP 
// Prototype for Fortran procedure piecewise_constant_interp ...
//
void FORTRAN_NAME( PIECEWISE_CONSTANT_INTERP ,piecewise_constant_interp )(
      CHFp_FRA1(refCoarData)
      ,CHFp_CONST_FRA1(coarData)
      ,CHFp_CONST_INT(refRat)
      ,CHFp_BOX(refBox)
      ,CHFp_BOX(coarBox) );

#define FORT_PIECEWISE_CONSTANT_INTERP FORTRAN_NAME( inlinePIECEWISE_CONSTANT_INTERP, inlinePIECEWISE_CONSTANT_INTERP)
#define FORTNT_PIECEWISE_CONSTANT_INTERP FORTRAN_NAME( PIECEWISE_CONSTANT_INTERP, piecewise_constant_interp)

inline void FORTRAN_NAME(inlinePIECEWISE_CONSTANT_INTERP, inlinePIECEWISE_CONSTANT_INTERP)(
      CHFp_FRA1(refCoarData)
      ,CHFp_CONST_FRA1(coarData)
      ,CHFp_CONST_INT(refRat)
      ,CHFp_BOX(refBox)
      ,CHFp_BOX(coarBox) )
{
 CH_TIMELEAF("FORT_PIECEWISE_CONSTANT_INTERP");
 FORTRAN_NAME( PIECEWISE_CONSTANT_INTERP ,piecewise_constant_interp )(
      CHFt_FRA1(refCoarData)
      ,CHFt_CONST_FRA1(coarData)
      ,CHFt_CONST_INT(refRat)
      ,CHFt_BOX(refBox)
      ,CHFt_BOX(coarBox) );
}
#endif  // GUARDPIECEWISE_CONSTANT_INTERP 

}

#endif
