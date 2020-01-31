#ifndef _EBLEVELDATAOPSF_F_H_
#define _EBLEVELDATAOPSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDEBAVECELLTOFACE 
#define GUARDEBAVECELLTOFACE 
// Prototype for Fortran procedure ebavecelltoface ...
//
void FORTRAN_NAME( EBAVECELLTOFACE ,ebavecelltoface )(
      CHFp_FRA1(facedata)
      ,CHFp_CONST_FRA1(celldata)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox) );

#define FORT_EBAVECELLTOFACE FORTRAN_NAME( inlineEBAVECELLTOFACE, inlineEBAVECELLTOFACE)
#define FORTNT_EBAVECELLTOFACE FORTRAN_NAME( EBAVECELLTOFACE, ebavecelltoface)

inline void FORTRAN_NAME(inlineEBAVECELLTOFACE, inlineEBAVECELLTOFACE)(
      CHFp_FRA1(facedata)
      ,CHFp_CONST_FRA1(celldata)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_EBAVECELLTOFACE");
 FORTRAN_NAME( EBAVECELLTOFACE ,ebavecelltoface )(
      CHFt_FRA1(facedata)
      ,CHFt_CONST_FRA1(celldata)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDEBAVECELLTOFACE 

}

#endif
