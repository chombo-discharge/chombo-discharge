#ifndef _CD_EBMULTIGRIDINTERPOLATORF_F_H_
#define _CD_EBMULTIGRIDINTERPOLATORF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMGINTERPHOMO 
#define GUARDMGINTERPHOMO 
// Prototype for Fortran procedure mginterphomo ...
//
void FORTRAN_NAME( MGINTERPHOMO ,mginterphomo )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_REAL(dxFine)
      ,CHFp_CONST_REAL(dxCoar)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(hilo)
      ,CHFp_BOX(ghostbox) );

#define FORT_MGINTERPHOMO FORTRAN_NAME( inlineMGINTERPHOMO, inlineMGINTERPHOMO)
#define FORTNT_MGINTERPHOMO FORTRAN_NAME( MGINTERPHOMO, mginterphomo)

inline void FORTRAN_NAME(inlineMGINTERPHOMO, inlineMGINTERPHOMO)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_REAL(dxFine)
      ,CHFp_CONST_REAL(dxCoar)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(hilo)
      ,CHFp_BOX(ghostbox) )
{
 CH_TIMELEAF("FORT_MGINTERPHOMO");
 FORTRAN_NAME( MGINTERPHOMO ,mginterphomo )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_REAL(dxFine)
      ,CHFt_CONST_REAL(dxCoar)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(hilo)
      ,CHFt_BOX(ghostbox) );
}
#endif  // GUARDMGINTERPHOMO 

}

#endif
