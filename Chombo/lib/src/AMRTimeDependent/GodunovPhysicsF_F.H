#ifndef _GODUNOVPHYSICSF_F_H_
#define _GODUNOVPHYSICSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDFLUXDIFFF 
#define GUARDFLUXDIFFF 
// Prototype for Fortran procedure FLUXDIFFF ...
//
void FORTRAN_NAME( FLUXDIFFF ,fluxdifff )(
      CHFp_FRA(diff)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_FLUXDIFFF FORTRAN_NAME( inlineFLUXDIFFF, inlineFLUXDIFFF)
#define FORTNT_FLUXDIFFF FORTRAN_NAME( FLUXDIFFF, fluxdifff)

inline void FORTRAN_NAME(inlineFLUXDIFFF, inlineFLUXDIFFF)(
      CHFp_FRA(diff)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FLUXDIFFF");
 FORTRAN_NAME( FLUXDIFFF ,fluxdifff )(
      CHFt_FRA(diff)
      ,CHFt_CONST_FRA(F)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFLUXDIFFF 

}

#endif
