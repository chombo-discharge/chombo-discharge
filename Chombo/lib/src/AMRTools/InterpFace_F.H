#ifndef _INTERPFACE_F_H_
#define _INTERPFACE_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDINTERPFACECONSTANT 
#define GUARDINTERPFACECONSTANT 
// Prototype for Fortran procedure interpfaceconstant ...
//
void FORTRAN_NAME( INTERPFACECONSTANT ,interpfaceconstant )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_BOX(bref)
      ,CHFp_CONST_INT(dir) );

#define FORT_INTERPFACECONSTANT FORTRAN_NAME( inlineINTERPFACECONSTANT, inlineINTERPFACECONSTANT)
#define FORTNT_INTERPFACECONSTANT FORTRAN_NAME( INTERPFACECONSTANT, interpfaceconstant)

inline void FORTRAN_NAME(inlineINTERPFACECONSTANT, inlineINTERPFACECONSTANT)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_BOX(bref)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_INTERPFACECONSTANT");
 FORTRAN_NAME( INTERPFACECONSTANT ,interpfaceconstant )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(ref_ratio)
      ,CHFt_BOX(bref)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDINTERPFACECONSTANT 

#ifndef GUARDINTERPLINEARFACE 
#define GUARDINTERPLINEARFACE 
// Prototype for Fortran procedure interplinearface ...
//
void FORTRAN_NAME( INTERPLINEARFACE ,interplinearface )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(slope)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_BOX(brefface) );

#define FORT_INTERPLINEARFACE FORTRAN_NAME( inlineINTERPLINEARFACE, inlineINTERPLINEARFACE)
#define FORTNT_INTERPLINEARFACE FORTRAN_NAME( INTERPLINEARFACE, interplinearface)

inline void FORTRAN_NAME(inlineINTERPLINEARFACE, inlineINTERPLINEARFACE)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(slope)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_BOX(brefface) )
{
 CH_TIMELEAF("FORT_INTERPLINEARFACE");
 FORTRAN_NAME( INTERPLINEARFACE ,interplinearface )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(slope)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(ref_ratio)
      ,CHFt_BOX(brefface) );
}
#endif  // GUARDINTERPLINEARFACE 

#ifndef GUARDINTERPLINEARINTERIORFACE 
#define GUARDINTERPLINEARINTERIORFACE 
// Prototype for Fortran procedure interplinearinteriorface ...
//
void FORTRAN_NAME( INTERPLINEARINTERIORFACE ,interplinearinteriorface )(
      CHFp_FRA(fine)
      ,CHFp_BOX(bcoarse)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(interiorrefbox) );

#define FORT_INTERPLINEARINTERIORFACE FORTRAN_NAME( inlineINTERPLINEARINTERIORFACE, inlineINTERPLINEARINTERIORFACE)
#define FORTNT_INTERPLINEARINTERIORFACE FORTRAN_NAME( INTERPLINEARINTERIORFACE, interplinearinteriorface)

inline void FORTRAN_NAME(inlineINTERPLINEARINTERIORFACE, inlineINTERPLINEARINTERIORFACE)(
      CHFp_FRA(fine)
      ,CHFp_BOX(bcoarse)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(interiorrefbox) )
{
 CH_TIMELEAF("FORT_INTERPLINEARINTERIORFACE");
 FORTRAN_NAME( INTERPLINEARINTERIORFACE ,interplinearinteriorface )(
      CHFt_FRA(fine)
      ,CHFt_BOX(bcoarse)
      ,CHFt_CONST_INT(ref_ratio)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(interiorrefbox) );
}
#endif  // GUARDINTERPLINEARINTERIORFACE 

#ifndef GUARDINTERPLIMITFACE 
#define GUARDINTERPLIMITFACE 
// Prototype for Fortran procedure interplimitface ...
//
void FORTRAN_NAME( INTERPLIMITFACE ,interplimitface )(
      CHFp_FRA(islope)
      ,CHFp_FRA(jslope)
      ,CHFp_FRA(kslope)
      ,CHFp_FRA(lslope)
      ,CHFp_FRA(mslope)
      ,CHFp_FRA(nslope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_BOX(bn)
      ,CHFp_BOX(validBox)
      ,CHFp_CONST_INT(normaldir) );

#define FORT_INTERPLIMITFACE FORTRAN_NAME( inlineINTERPLIMITFACE, inlineINTERPLIMITFACE)
#define FORTNT_INTERPLIMITFACE FORTRAN_NAME( INTERPLIMITFACE, interplimitface)

inline void FORTRAN_NAME(inlineINTERPLIMITFACE, inlineINTERPLIMITFACE)(
      CHFp_FRA(islope)
      ,CHFp_FRA(jslope)
      ,CHFp_FRA(kslope)
      ,CHFp_FRA(lslope)
      ,CHFp_FRA(mslope)
      ,CHFp_FRA(nslope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_BOX(bn)
      ,CHFp_BOX(validBox)
      ,CHFp_CONST_INT(normaldir) )
{
 CH_TIMELEAF("FORT_INTERPLIMITFACE");
 FORTRAN_NAME( INTERPLIMITFACE ,interplimitface )(
      CHFt_FRA(islope)
      ,CHFt_FRA(jslope)
      ,CHFt_FRA(kslope)
      ,CHFt_FRA(lslope)
      ,CHFt_FRA(mslope)
      ,CHFt_FRA(nslope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(b)
      ,CHFt_BOX(bn)
      ,CHFt_BOX(validBox)
      ,CHFt_CONST_INT(normaldir) );
}
#endif  // GUARDINTERPLIMITFACE 

}

#endif
