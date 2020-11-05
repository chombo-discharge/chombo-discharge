#ifndef _BASEFABINTPLUSF_F_H_
#define _BASEFABINTPLUSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDBASEFABINTPLUS 
#define GUARDBASEFABINTPLUS 
// Prototype for Fortran procedure basefabintplus ...
//
void FORTRAN_NAME( BASEFABINTPLUS ,basefabintplus )(
      CHFp_FIA(sum)
      ,CHFp_CONST_FIA(piece)
      ,CHFp_BOX(bx) );

#define FORT_BASEFABINTPLUS FORTRAN_NAME( inlineBASEFABINTPLUS, inlineBASEFABINTPLUS)
#define FORTNT_BASEFABINTPLUS FORTRAN_NAME( BASEFABINTPLUS, basefabintplus)

inline void FORTRAN_NAME(inlineBASEFABINTPLUS, inlineBASEFABINTPLUS)(
      CHFp_FIA(sum)
      ,CHFp_CONST_FIA(piece)
      ,CHFp_BOX(bx) )
{
 CH_TIMELEAF("FORT_BASEFABINTPLUS");
 FORTRAN_NAME( BASEFABINTPLUS ,basefabintplus )(
      CHFt_FIA(sum)
      ,CHFt_CONST_FIA(piece)
      ,CHFt_BOX(bx) );
}
#endif  // GUARDBASEFABINTPLUS 

}

#endif
