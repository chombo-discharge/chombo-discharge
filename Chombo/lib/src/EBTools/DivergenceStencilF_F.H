#ifndef _DIVERGENCESTENCILF_F_H_
#define _DIVERGENCESTENCILF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDDIVERGESTENF 
#define GUARDDIVERGESTENF 
// Prototype for Fortran procedure divergestenf ...
//
void FORTRAN_NAME( DIVERGESTENF ,divergestenf )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(divf)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dx) );

#define FORT_DIVERGESTENF FORTRAN_NAME( inlineDIVERGESTENF, inlineDIVERGESTENF)
#define FORTNT_DIVERGESTENF FORTRAN_NAME( DIVERGESTENF, divergestenf)

inline void FORTRAN_NAME(inlineDIVERGESTENF, inlineDIVERGESTENF)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(divf)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_DIVERGESTENF");
 FORTRAN_NAME( DIVERGESTENF ,divergestenf )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(divf)
      ,CHFt_CONST_FRA(flux)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(nconserved)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDDIVERGESTENF 

}

#endif
