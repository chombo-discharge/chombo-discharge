#ifndef _QUADCFINTERPF_F_H_
#define _QUADCFINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDQUADINTERP 
#define GUARDQUADINTERP 
// Prototype for Fortran procedure quadinterp ...
//
void FORTRAN_NAME( QUADINTERP ,quadinterp )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(phistar)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(ihilo)
      ,CHFp_CONST_REAL(h)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(scomp)
      ,CHFp_CONST_INT(ecomp)
      ,CHFp_CONST_INT(nref) );

#define FORT_QUADINTERP FORTRAN_NAME( inlineQUADINTERP, inlineQUADINTERP)
#define FORTNT_QUADINTERP FORTRAN_NAME( QUADINTERP, quadinterp)

inline void FORTRAN_NAME(inlineQUADINTERP, inlineQUADINTERP)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(phistar)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(ihilo)
      ,CHFp_CONST_REAL(h)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(scomp)
      ,CHFp_CONST_INT(ecomp)
      ,CHFp_CONST_INT(nref) )
{
 CH_TIMELEAF("FORT_QUADINTERP");
 FORTRAN_NAME( QUADINTERP ,quadinterp )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(phistar)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(ihilo)
      ,CHFt_CONST_REAL(h)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(scomp)
      ,CHFt_CONST_INT(ecomp)
      ,CHFt_CONST_INT(nref) );
}
#endif  // GUARDQUADINTERP 

#ifndef GUARDPHISTAR 
#define GUARDPHISTAR 
// Prototype for Fortran procedure PHISTAR ...
//
void FORTRAN_NAME( PHISTAR ,phistar )(
      CHFp_FRA(fPhiStar)
      ,CHFp_BOX(region)
      ,CHFp_CONST_FRA(phic)
      ,CHFp_FRA(coarslope)
      ,CHFp_FRA(coarcurva)
      ,CHFp_FRA(coarmixed)
      ,CHFp_CONST_REAL(dxf)
      ,CHFp_CONST_INT(ivar)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(sign)
      ,CHFp_CONST_INT(nRef) );

#define FORT_PHISTAR FORTRAN_NAME( inlinePHISTAR, inlinePHISTAR)
#define FORTNT_PHISTAR FORTRAN_NAME( PHISTAR, phistar)

inline void FORTRAN_NAME(inlinePHISTAR, inlinePHISTAR)(
      CHFp_FRA(fPhiStar)
      ,CHFp_BOX(region)
      ,CHFp_CONST_FRA(phic)
      ,CHFp_FRA(coarslope)
      ,CHFp_FRA(coarcurva)
      ,CHFp_FRA(coarmixed)
      ,CHFp_CONST_REAL(dxf)
      ,CHFp_CONST_INT(ivar)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(sign)
      ,CHFp_CONST_INT(nRef) )
{
 CH_TIMELEAF("FORT_PHISTAR");
 FORTRAN_NAME( PHISTAR ,phistar )(
      CHFt_FRA(fPhiStar)
      ,CHFt_BOX(region)
      ,CHFt_CONST_FRA(phic)
      ,CHFt_FRA(coarslope)
      ,CHFt_FRA(coarcurva)
      ,CHFt_FRA(coarmixed)
      ,CHFt_CONST_REAL(dxf)
      ,CHFt_CONST_INT(ivar)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(sign)
      ,CHFt_CONST_INT(nRef) );
}
#endif  // GUARDPHISTAR 

}

#endif
