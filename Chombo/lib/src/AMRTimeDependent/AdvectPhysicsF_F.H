#ifndef _ADVECTPHYSICSF_F_H_
#define _ADVECTPHYSICSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDAVEFACESCALTOFACEVECT 
#define GUARDAVEFACESCALTOFACEVECT 
// Prototype for Fortran procedure avefacescaltofacevect ...
//
void FORTRAN_NAME( AVEFACESCALTOFACEVECT ,avefacescaltofacevect )(
      CHFp_FRA(facevect)
      ,CHFp_CONST_FRA1(facescal)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(vectdir)
      ,CHFp_BOX(dcalcface)
      ,CHFp_BOX(offbox) );

#define FORT_AVEFACESCALTOFACEVECT FORTRAN_NAME( inlineAVEFACESCALTOFACEVECT, inlineAVEFACESCALTOFACEVECT)
#define FORTNT_AVEFACESCALTOFACEVECT FORTRAN_NAME( AVEFACESCALTOFACEVECT, avefacescaltofacevect)

inline void FORTRAN_NAME(inlineAVEFACESCALTOFACEVECT, inlineAVEFACESCALTOFACEVECT)(
      CHFp_FRA(facevect)
      ,CHFp_CONST_FRA1(facescal)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(vectdir)
      ,CHFp_BOX(dcalcface)
      ,CHFp_BOX(offbox) )
{
 CH_TIMELEAF("FORT_AVEFACESCALTOFACEVECT");
 FORTRAN_NAME( AVEFACESCALTOFACEVECT ,avefacescaltofacevect )(
      CHFt_FRA(facevect)
      ,CHFt_CONST_FRA1(facescal)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(vectdir)
      ,CHFt_BOX(dcalcface)
      ,CHFt_BOX(offbox) );
}
#endif  // GUARDAVEFACESCALTOFACEVECT 

#ifndef GUARDAVESCALTOFACE 
#define GUARDAVESCALTOFACE 
// Prototype for Fortran procedure avescaltoface ...
//
void FORTRAN_NAME( AVESCALTOFACE ,avescaltoface )(
      CHFp_FRA1(facescal)
      ,CHFp_CONST_FRA1(cellscal)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(dcalcface) );

#define FORT_AVESCALTOFACE FORTRAN_NAME( inlineAVESCALTOFACE, inlineAVESCALTOFACE)
#define FORTNT_AVESCALTOFACE FORTRAN_NAME( AVESCALTOFACE, avescaltoface)

inline void FORTRAN_NAME(inlineAVESCALTOFACE, inlineAVESCALTOFACE)(
      CHFp_FRA1(facescal)
      ,CHFp_CONST_FRA1(cellscal)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(dcalcface) )
{
 CH_TIMELEAF("FORT_AVESCALTOFACE");
 FORTRAN_NAME( AVESCALTOFACE ,avescaltoface )(
      CHFt_FRA1(facescal)
      ,CHFt_CONST_FRA1(cellscal)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(dcalcface) );
}
#endif  // GUARDAVESCALTOFACE 

#ifndef GUARDAVECELLTOFACE 
#define GUARDAVECELLTOFACE 
// Prototype for Fortran procedure avecelltoface ...
//
void FORTRAN_NAME( AVECELLTOFACE ,avecelltoface )(
      CHFp_FRA1(facevel)
      ,CHFp_CONST_FRA(cellvel)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(dcalcface) );

#define FORT_AVECELLTOFACE FORTRAN_NAME( inlineAVECELLTOFACE, inlineAVECELLTOFACE)
#define FORTNT_AVECELLTOFACE FORTRAN_NAME( AVECELLTOFACE, avecelltoface)

inline void FORTRAN_NAME(inlineAVECELLTOFACE, inlineAVECELLTOFACE)(
      CHFp_FRA1(facevel)
      ,CHFp_CONST_FRA(cellvel)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(dcalcface) )
{
 CH_TIMELEAF("FORT_AVECELLTOFACE");
 FORTRAN_NAME( AVECELLTOFACE ,avecelltoface )(
      CHFt_FRA1(facevel)
      ,CHFt_CONST_FRA(cellvel)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(dcalcface) );
}
#endif  // GUARDAVECELLTOFACE 

#ifndef GUARDAVEFACETOCELL 
#define GUARDAVEFACETOCELL 
// Prototype for Fortran procedure avefacetocell ...
//
void FORTRAN_NAME( AVEFACETOCELL ,avefacetocell )(
      CHFp_FRA(cellvel)
      ,CHFp_CONST_FRA1(facevel)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(dcalccell) );

#define FORT_AVEFACETOCELL FORTRAN_NAME( inlineAVEFACETOCELL, inlineAVEFACETOCELL)
#define FORTNT_AVEFACETOCELL FORTRAN_NAME( AVEFACETOCELL, avefacetocell)

inline void FORTRAN_NAME(inlineAVEFACETOCELL, inlineAVEFACETOCELL)(
      CHFp_FRA(cellvel)
      ,CHFp_CONST_FRA1(facevel)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(dcalccell) )
{
 CH_TIMELEAF("FORT_AVEFACETOCELL");
 FORTRAN_NAME( AVEFACETOCELL ,avefacetocell )(
      CHFt_FRA(cellvel)
      ,CHFt_CONST_FRA1(facevel)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(dcalccell) );
}
#endif  // GUARDAVEFACETOCELL 

#ifndef GUARDMAGNITUDEF 
#define GUARDMAGNITUDEF 
// Prototype for Fortran procedure MAGNITUDEF ...
//
void FORTRAN_NAME( MAGNITUDEF ,magnitudef )(
      CHFp_FRA1(magdata)
      ,CHFp_CONST_FRA(data)
      ,CHFp_BOX(box) );

#define FORT_MAGNITUDEF FORTRAN_NAME( inlineMAGNITUDEF, inlineMAGNITUDEF)
#define FORTNT_MAGNITUDEF FORTRAN_NAME( MAGNITUDEF, magnitudef)

inline void FORTRAN_NAME(inlineMAGNITUDEF, inlineMAGNITUDEF)(
      CHFp_FRA1(magdata)
      ,CHFp_CONST_FRA(data)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAGNITUDEF");
 FORTRAN_NAME( MAGNITUDEF ,magnitudef )(
      CHFt_FRA1(magdata)
      ,CHFt_CONST_FRA(data)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAGNITUDEF 

#ifndef GUARDGETRELGRADF 
#define GUARDGETRELGRADF 
// Prototype for Fortran procedure GETRELGRADF ...
//
void FORTRAN_NAME( GETRELGRADF ,getrelgradf )(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_GETRELGRADF FORTRAN_NAME( inlineGETRELGRADF, inlineGETRELGRADF)
#define FORTNT_GETRELGRADF FORTRAN_NAME( GETRELGRADF, getrelgradf)

inline void FORTRAN_NAME(inlineGETRELGRADF, inlineGETRELGRADF)(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_GETRELGRADF");
 FORTRAN_NAME( GETRELGRADF ,getrelgradf )(
      CHFt_FRA1(du)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDGETRELGRADF 

#ifndef GUARDPOSTNORMALSOURCE 
#define GUARDPOSTNORMALSOURCE 
// Prototype for Fortran procedure POSTNORMALSOURCE ...
//
void FORTRAN_NAME( POSTNORMALSOURCE ,postnormalsource )(
      CHFp_FRA(dWminus)
      ,CHFp_FRA(dWplus)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA1(advVel)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_POSTNORMALSOURCE FORTRAN_NAME( inlinePOSTNORMALSOURCE, inlinePOSTNORMALSOURCE)
#define FORTNT_POSTNORMALSOURCE FORTRAN_NAME( POSTNORMALSOURCE, postnormalsource)

inline void FORTRAN_NAME(inlinePOSTNORMALSOURCE, inlinePOSTNORMALSOURCE)(
      CHFp_FRA(dWminus)
      ,CHFp_FRA(dWplus)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA1(advVel)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_POSTNORMALSOURCE");
 FORTRAN_NAME( POSTNORMALSOURCE ,postnormalsource )(
      CHFt_FRA(dWminus)
      ,CHFt_FRA(dWplus)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA1(advVel)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPOSTNORMALSOURCE 

#ifndef GUARDRIEMANNF 
#define GUARDRIEMANNF 
// Prototype for Fortran procedure RIEMANNF ...
//
void FORTRAN_NAME( RIEMANNF ,riemannf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_FRA1(advVel)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_RIEMANNF FORTRAN_NAME( inlineRIEMANNF, inlineRIEMANNF)
#define FORTNT_RIEMANNF FORTRAN_NAME( RIEMANNF, riemannf)

inline void FORTRAN_NAME(inlineRIEMANNF, inlineRIEMANNF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_FRA1(advVel)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RIEMANNF");
 FORTRAN_NAME( RIEMANNF ,riemannf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_FRA1(advVel)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRIEMANNF 

#ifndef GUARDQUASILINEARUPDATE 
#define GUARDQUASILINEARUPDATE 
// Prototype for Fortran procedure QUASILINEARUPDATE ...
//
void FORTRAN_NAME( QUASILINEARUPDATE ,quasilinearupdate )(
      CHFp_FRA(AdWdx)
      ,CHFp_CONST_FRA(WHalf)
      ,CHFp_CONST_FRA1(cellVel)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_QUASILINEARUPDATE FORTRAN_NAME( inlineQUASILINEARUPDATE, inlineQUASILINEARUPDATE)
#define FORTNT_QUASILINEARUPDATE FORTRAN_NAME( QUASILINEARUPDATE, quasilinearupdate)

inline void FORTRAN_NAME(inlineQUASILINEARUPDATE, inlineQUASILINEARUPDATE)(
      CHFp_FRA(AdWdx)
      ,CHFp_CONST_FRA(WHalf)
      ,CHFp_CONST_FRA1(cellVel)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_QUASILINEARUPDATE");
 FORTRAN_NAME( QUASILINEARUPDATE ,quasilinearupdate )(
      CHFt_FRA(AdWdx)
      ,CHFt_CONST_FRA(WHalf)
      ,CHFt_CONST_FRA1(cellVel)
      ,CHFt_CONST_REAL(scale)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDQUASILINEARUPDATE 

#ifndef GUARDAVECELLVECTOFACEVEC 
#define GUARDAVECELLVECTOFACEVEC 
// Prototype for Fortran procedure avecellvectofacevec ...
//
void FORTRAN_NAME( AVECELLVECTOFACEVEC ,avecellvectofacevec )(
      CHFp_FRA(facevec)
      ,CHFp_CONST_FRA(cellvec)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(dcalcface) );

#define FORT_AVECELLVECTOFACEVEC FORTRAN_NAME( inlineAVECELLVECTOFACEVEC, inlineAVECELLVECTOFACEVEC)
#define FORTNT_AVECELLVECTOFACEVEC FORTRAN_NAME( AVECELLVECTOFACEVEC, avecellvectofacevec)

inline void FORTRAN_NAME(inlineAVECELLVECTOFACEVEC, inlineAVECELLVECTOFACEVEC)(
      CHFp_FRA(facevec)
      ,CHFp_CONST_FRA(cellvec)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(dcalcface) )
{
 CH_TIMELEAF("FORT_AVECELLVECTOFACEVEC");
 FORTRAN_NAME( AVECELLVECTOFACEVEC ,avecellvectofacevec )(
      CHFt_FRA(facevec)
      ,CHFt_CONST_FRA(cellvec)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(dcalcface) );
}
#endif  // GUARDAVECELLVECTOFACEVEC 

}

#endif
