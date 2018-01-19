#ifndef _CDR_SOLVERF_F_H_
#define _CDR_SOLVERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDADVECTIVEDERIV 
#define GUARDADVECTIVEDERIV 
// Prototype for Fortran procedure advectivederiv ...
//
void FORTRAN_NAME( ADVECTIVEDERIV ,advectivederiv )(
      CHFp_FRA1(udelrho)
      ,CHFp_CONST_FRA1(facerho)
      ,CHFp_CONST_FRA1(facevel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc) );

#define FORT_ADVECTIVEDERIV FORTRAN_NAME( inlineADVECTIVEDERIV, inlineADVECTIVEDERIV)
#define FORTNT_ADVECTIVEDERIV FORTRAN_NAME( ADVECTIVEDERIV, advectivederiv)

inline void FORTRAN_NAME(inlineADVECTIVEDERIV, inlineADVECTIVEDERIV)(
      CHFp_FRA1(udelrho)
      ,CHFp_CONST_FRA1(facerho)
      ,CHFp_CONST_FRA1(facevel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(ncomp)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_ADVECTIVEDERIV");
 FORTRAN_NAME( ADVECTIVEDERIV ,advectivederiv )(
      CHFt_FRA1(udelrho)
      ,CHFt_CONST_FRA1(facerho)
      ,CHFt_CONST_FRA1(facevel)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(ncomp)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDADVECTIVEDERIV 

}

#endif
