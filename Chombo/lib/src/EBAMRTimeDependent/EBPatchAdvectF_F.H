#ifndef _EBPATCHADVECTF_F_H_
#define _EBPATCHADVECTF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDADVECTSLOPEDIFFS 
#define GUARDADVECTSLOPEDIFFS 
// Prototype for Fortran procedure advectslopediffs ...
//
void FORTRAN_NAME( ADVECTSLOPEDIFFS ,advectslopediffs )(
      CHFp_FRA(deltawc)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox)
      ,CHFp_INT(iuselimiting) );

#define FORT_ADVECTSLOPEDIFFS FORTRAN_NAME( inlineADVECTSLOPEDIFFS, inlineADVECTSLOPEDIFFS)
#define FORTNT_ADVECTSLOPEDIFFS FORTRAN_NAME( ADVECTSLOPEDIFFS, advectslopediffs)

inline void FORTRAN_NAME(inlineADVECTSLOPEDIFFS, inlineADVECTSLOPEDIFFS)(
      CHFp_FRA(deltawc)
      ,CHFp_CONST_FRA(w)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox)
      ,CHFp_INT(iuselimiting) )
{
 CH_TIMELEAF("FORT_ADVECTSLOPEDIFFS");
 FORTRAN_NAME( ADVECTSLOPEDIFFS ,advectslopediffs )(
      CHFt_FRA(deltawc)
      ,CHFt_CONST_FRA(w)
      ,CHFt_CONST_INT(numslopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox)
      ,CHFt_INT(iuselimiting) );
}
#endif  // GUARDADVECTSLOPEDIFFS 

#ifndef GUARDPREDADVECT 
#define GUARDPREDADVECT 
// Prototype for Fortran procedure predadvect ...
//
void FORTRAN_NAME( PREDADVECT ,predadvect )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA1(rho)
      ,CHFp_CONST_FRA1(drho)
      ,CHFp_CONST_FRA(velcc)
      ,CHFp_FRA1(rholo)
      ,CHFp_FRA1(rhohi)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx) );

#define FORT_PREDADVECT FORTRAN_NAME( inlinePREDADVECT, inlinePREDADVECT)
#define FORTNT_PREDADVECT FORTRAN_NAME( PREDADVECT, predadvect)

inline void FORTRAN_NAME(inlinePREDADVECT, inlinePREDADVECT)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA1(rho)
      ,CHFp_CONST_FRA1(drho)
      ,CHFp_CONST_FRA(velcc)
      ,CHFp_FRA1(rholo)
      ,CHFp_FRA1(rhohi)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx) )
{
 CH_TIMELEAF("FORT_PREDADVECT");
 FORTRAN_NAME( PREDADVECT ,predadvect )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_FRA1(rho)
      ,CHFt_CONST_FRA1(drho)
      ,CHFt_CONST_FRA(velcc)
      ,CHFt_FRA1(rholo)
      ,CHFt_FRA1(rhohi)
      ,CHFt_CONST_INT(normdir)
      ,CHFt_CONST_REAL(dtbydx) );
}
#endif  // GUARDPREDADVECT 

#ifndef GUARDPOINTPREDADVECT 
#define GUARDPOINTPREDADVECT 
// Prototype for Fortran procedure pointpredadvect ...
//
void FORTRAN_NAME( POINTPREDADVECT ,pointpredadvect )(
      CHFp_REAL(dense)
      ,CHFp_REAL(denlo)
      ,CHFp_REAL(denhi)
      ,CHFp_REAL(denslope)
      ,CHFp_REALVECT(veloc)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx) );

#define FORT_POINTPREDADVECT FORTRAN_NAME( inlinePOINTPREDADVECT, inlinePOINTPREDADVECT)
#define FORTNT_POINTPREDADVECT FORTRAN_NAME( POINTPREDADVECT, pointpredadvect)

inline void FORTRAN_NAME(inlinePOINTPREDADVECT, inlinePOINTPREDADVECT)(
      CHFp_REAL(dense)
      ,CHFp_REAL(denlo)
      ,CHFp_REAL(denhi)
      ,CHFp_REAL(denslope)
      ,CHFp_REALVECT(veloc)
      ,CHFp_CONST_INT(normdir)
      ,CHFp_CONST_REAL(dtbydx) )
{
 CH_TIMELEAF("FORT_POINTPREDADVECT");
 FORTRAN_NAME( POINTPREDADVECT ,pointpredadvect )(
      CHFt_REAL(dense)
      ,CHFt_REAL(denlo)
      ,CHFt_REAL(denhi)
      ,CHFt_REAL(denslope)
      ,CHFt_REALVECT(veloc)
      ,CHFt_CONST_INT(normdir)
      ,CHFt_CONST_REAL(dtbydx) );
}
#endif  // GUARDPOINTPREDADVECT 

#ifndef GUARDPREDADVECTTRANS 
#define GUARDPREDADVECTTRANS 
// Prototype for Fortran procedure predadvecttrans ...
//
void FORTRAN_NAME( PREDADVECTTRANS ,predadvecttrans )(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA1(rho)
      ,CHFp_CONST_FRA1(drho)
      ,CHFp_CONST_FRA(velcc)
      ,CHFp_FRA1(rholo)
      ,CHFp_FRA1(rhohi)
      ,CHFp_CONST_INT(tandir)
      ,CHFp_CONST_REAL(dtbydx) );

#define FORT_PREDADVECTTRANS FORTRAN_NAME( inlinePREDADVECTTRANS, inlinePREDADVECTTRANS)
#define FORTNT_PREDADVECTTRANS FORTRAN_NAME( PREDADVECTTRANS, predadvecttrans)

inline void FORTRAN_NAME(inlinePREDADVECTTRANS, inlinePREDADVECTTRANS)(
      CHFp_BOX(dcalc)
      ,CHFp_CONST_FRA1(rho)
      ,CHFp_CONST_FRA1(drho)
      ,CHFp_CONST_FRA(velcc)
      ,CHFp_FRA1(rholo)
      ,CHFp_FRA1(rhohi)
      ,CHFp_CONST_INT(tandir)
      ,CHFp_CONST_REAL(dtbydx) )
{
 CH_TIMELEAF("FORT_PREDADVECTTRANS");
 FORTRAN_NAME( PREDADVECTTRANS ,predadvecttrans )(
      CHFt_BOX(dcalc)
      ,CHFt_CONST_FRA1(rho)
      ,CHFt_CONST_FRA1(drho)
      ,CHFt_CONST_FRA(velcc)
      ,CHFt_FRA1(rholo)
      ,CHFt_FRA1(rhohi)
      ,CHFt_CONST_INT(tandir)
      ,CHFt_CONST_REAL(dtbydx) );
}
#endif  // GUARDPREDADVECTTRANS 

#ifndef GUARDPOINTPREDADVECTTRANS 
#define GUARDPOINTPREDADVECTTRANS 
// Prototype for Fortran procedure pointpredadvecttrans ...
//
void FORTRAN_NAME( POINTPREDADVECTTRANS ,pointpredadvecttrans )(
      CHFp_REAL(dense)
      ,CHFp_REAL(denlo)
      ,CHFp_REAL(denhi)
      ,CHFp_REAL(denslope)
      ,CHFp_REALVECT(veloc)
      ,CHFp_CONST_INT(tandir)
      ,CHFp_CONST_REAL(dtbydx) );

#define FORT_POINTPREDADVECTTRANS FORTRAN_NAME( inlinePOINTPREDADVECTTRANS, inlinePOINTPREDADVECTTRANS)
#define FORTNT_POINTPREDADVECTTRANS FORTRAN_NAME( POINTPREDADVECTTRANS, pointpredadvecttrans)

inline void FORTRAN_NAME(inlinePOINTPREDADVECTTRANS, inlinePOINTPREDADVECTTRANS)(
      CHFp_REAL(dense)
      ,CHFp_REAL(denlo)
      ,CHFp_REAL(denhi)
      ,CHFp_REAL(denslope)
      ,CHFp_REALVECT(veloc)
      ,CHFp_CONST_INT(tandir)
      ,CHFp_CONST_REAL(dtbydx) )
{
 CH_TIMELEAF("FORT_POINTPREDADVECTTRANS");
 FORTRAN_NAME( POINTPREDADVECTTRANS ,pointpredadvecttrans )(
      CHFt_REAL(dense)
      ,CHFt_REAL(denlo)
      ,CHFt_REAL(denhi)
      ,CHFt_REAL(denslope)
      ,CHFt_REALVECT(veloc)
      ,CHFt_CONST_INT(tandir)
      ,CHFt_CONST_REAL(dtbydx) );
}
#endif  // GUARDPOINTPREDADVECTTRANS 

#ifndef GUARDGETDSDT 
#define GUARDGETDSDT 
// Prototype for Fortran procedure GETDSDT ...
//
void FORTRAN_NAME( GETDSDT ,getdsdt )(
      CHFp_FRA(dsdtplus)
      ,CHFp_FRA(dsdtminu)
      ,CHFp_CONST_FRA(slopeprim)
      ,CHFp_CONST_FRA(slopeupwi)
      ,CHFp_CONST_FRA(normalvel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(extrapdir)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_BOX(entirebox) );

#define FORT_GETDSDT FORTRAN_NAME( inlineGETDSDT, inlineGETDSDT)
#define FORTNT_GETDSDT FORTRAN_NAME( GETDSDT, getdsdt)

inline void FORTRAN_NAME(inlineGETDSDT, inlineGETDSDT)(
      CHFp_FRA(dsdtplus)
      ,CHFp_FRA(dsdtminu)
      ,CHFp_CONST_FRA(slopeprim)
      ,CHFp_CONST_FRA(slopeupwi)
      ,CHFp_CONST_FRA(normalvel)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(extrapdir)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_BOX(entirebox) )
{
 CH_TIMELEAF("FORT_GETDSDT");
 FORTRAN_NAME( GETDSDT ,getdsdt )(
      CHFt_FRA(dsdtplus)
      ,CHFt_FRA(dsdtminu)
      ,CHFt_CONST_FRA(slopeprim)
      ,CHFt_CONST_FRA(slopeupwi)
      ,CHFt_CONST_FRA(normalvel)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(extrapdir)
      ,CHFt_CONST_INT(numslopes)
      ,CHFt_BOX(entirebox) );
}
#endif  // GUARDGETDSDT 

#ifndef GUARDUPWINDDIFFS 
#define GUARDUPWINDDIFFS 
// Prototype for Fortran procedure upwinddiffs ...
//
void FORTRAN_NAME( UPWINDDIFFS ,upwinddiffs )(
      CHFp_FRA(slopeupwi)
      ,CHFp_CONST_FRA(primstate)
      ,CHFp_CONST_FRA(normalvel)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) );

#define FORT_UPWINDDIFFS FORTRAN_NAME( inlineUPWINDDIFFS, inlineUPWINDDIFFS)
#define FORTNT_UPWINDDIFFS FORTRAN_NAME( UPWINDDIFFS, upwinddiffs)

inline void FORTRAN_NAME(inlineUPWINDDIFFS, inlineUPWINDDIFFS)(
      CHFp_FRA(slopeupwi)
      ,CHFp_CONST_FRA(primstate)
      ,CHFp_CONST_FRA(normalvel)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(lobox)
      ,CHFp_CONST_INT(haslo)
      ,CHFp_BOX(hibox)
      ,CHFp_CONST_INT(hashi)
      ,CHFp_BOX(centerbox) )
{
 CH_TIMELEAF("FORT_UPWINDDIFFS");
 FORTRAN_NAME( UPWINDDIFFS ,upwinddiffs )(
      CHFt_FRA(slopeupwi)
      ,CHFt_CONST_FRA(primstate)
      ,CHFt_CONST_FRA(normalvel)
      ,CHFt_CONST_INT(numslopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(lobox)
      ,CHFt_CONST_INT(haslo)
      ,CHFt_BOX(hibox)
      ,CHFt_CONST_INT(hashi)
      ,CHFt_BOX(centerbox) );
}
#endif  // GUARDUPWINDDIFFS 

#ifndef GUARDADVECTIVEF 
#define GUARDADVECTIVEF 
// Prototype for Fortran procedure advectivef ...
//
void FORTRAN_NAME( ADVECTIVEF ,advectivef )(
      CHFp_FRA(udelrho)
      ,CHFp_CONST_FRA(facerho)
      ,CHFp_CONST_FRA1(facevel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc)
      ,CHFp_INT(doingvel) );

#define FORT_ADVECTIVEF FORTRAN_NAME( inlineADVECTIVEF, inlineADVECTIVEF)
#define FORTNT_ADVECTIVEF FORTRAN_NAME( ADVECTIVEF, advectivef)

inline void FORTRAN_NAME(inlineADVECTIVEF, inlineADVECTIVEF)(
      CHFp_FRA(udelrho)
      ,CHFp_CONST_FRA(facerho)
      ,CHFp_CONST_FRA1(facevel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nconserved)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(dcalc)
      ,CHFp_INT(doingvel) )
{
 CH_TIMELEAF("FORT_ADVECTIVEF");
 FORTRAN_NAME( ADVECTIVEF ,advectivef )(
      CHFt_FRA(udelrho)
      ,CHFt_CONST_FRA(facerho)
      ,CHFt_CONST_FRA1(facevel)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(nconserved)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(dcalc)
      ,CHFt_INT(doingvel) );
}
#endif  // GUARDADVECTIVEF 

#ifndef GUARDEBAVEFACETOCELL 
#define GUARDEBAVEFACETOCELL 
// Prototype for Fortran procedure ebavefacetocell ...
//
void FORTRAN_NAME( EBAVEFACETOCELL ,ebavefacetocell )(
      CHFp_FRA(cellvel)
      ,CHFp_CONST_FRA1(facevel)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(dcalc) );

#define FORT_EBAVEFACETOCELL FORTRAN_NAME( inlineEBAVEFACETOCELL, inlineEBAVEFACETOCELL)
#define FORTNT_EBAVEFACETOCELL FORTRAN_NAME( EBAVEFACETOCELL, ebavefacetocell)

inline void FORTRAN_NAME(inlineEBAVEFACETOCELL, inlineEBAVEFACETOCELL)(
      CHFp_FRA(cellvel)
      ,CHFp_CONST_FRA1(facevel)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(dcalc) )
{
 CH_TIMELEAF("FORT_EBAVEFACETOCELL");
 FORTRAN_NAME( EBAVEFACETOCELL ,ebavefacetocell )(
      CHFt_FRA(cellvel)
      ,CHFt_CONST_FRA1(facevel)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(dcalc) );
}
#endif  // GUARDEBAVEFACETOCELL 

#ifndef GUARDADVECTUPDATE 
#define GUARDADVECTUPDATE 
// Prototype for Fortran procedure advectupdate ...
//
void FORTRAN_NAME( ADVECTUPDATE ,advectupdate )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(primminu)
      ,CHFp_FRA(primplus)
      ,CHFp_CONST_FRA(primface)
      ,CHFp_CONST_FRA(normvel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nprim)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_BOX(cellbox) );

#define FORT_ADVECTUPDATE FORTRAN_NAME( inlineADVECTUPDATE, inlineADVECTUPDATE)
#define FORTNT_ADVECTUPDATE FORTRAN_NAME( ADVECTUPDATE, advectupdate)

inline void FORTRAN_NAME(inlineADVECTUPDATE, inlineADVECTUPDATE)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(primminu)
      ,CHFp_FRA(primplus)
      ,CHFp_CONST_FRA(primface)
      ,CHFp_CONST_FRA(normvel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nprim)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_BOX(cellbox) )
{
 CH_TIMELEAF("FORT_ADVECTUPDATE");
 FORTRAN_NAME( ADVECTUPDATE ,advectupdate )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(primminu)
      ,CHFt_FRA(primplus)
      ,CHFt_CONST_FRA(primface)
      ,CHFt_CONST_FRA(normvel)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(nprim)
      ,CHFt_CONST_REAL(dtbydx)
      ,CHFt_BOX(cellbox) );
}
#endif  // GUARDADVECTUPDATE 

#ifndef GUARDADVECTRIEMANN 
#define GUARDADVECTRIEMANN 
// Prototype for Fortran procedure advectriemann ...
//
void FORTRAN_NAME( ADVECTRIEMANN ,advectriemann )(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(primgdnv)
      ,CHFp_CONST_FRA(primleft)
      ,CHFp_CONST_FRA(primrigh)
      ,CHFp_CONST_FRA1(advectvel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nprim)
      ,CHFp_CONST_INT(curcomp)
      ,CHFp_CONST_INT(doingvel) );

#define FORT_ADVECTRIEMANN FORTRAN_NAME( inlineADVECTRIEMANN, inlineADVECTRIEMANN)
#define FORTNT_ADVECTRIEMANN FORTRAN_NAME( ADVECTRIEMANN, advectriemann)

inline void FORTRAN_NAME(inlineADVECTRIEMANN, inlineADVECTRIEMANN)(
      CHFp_BOX(dcalc)
      ,CHFp_FRA(primgdnv)
      ,CHFp_CONST_FRA(primleft)
      ,CHFp_CONST_FRA(primrigh)
      ,CHFp_CONST_FRA1(advectvel)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_CONST_INT(nprim)
      ,CHFp_CONST_INT(curcomp)
      ,CHFp_CONST_INT(doingvel) )
{
 CH_TIMELEAF("FORT_ADVECTRIEMANN");
 FORTRAN_NAME( ADVECTRIEMANN ,advectriemann )(
      CHFt_BOX(dcalc)
      ,CHFt_FRA(primgdnv)
      ,CHFt_CONST_FRA(primleft)
      ,CHFt_CONST_FRA(primrigh)
      ,CHFt_CONST_FRA1(advectvel)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_CONST_INT(nprim)
      ,CHFt_CONST_INT(curcomp)
      ,CHFt_CONST_INT(doingvel) );
}
#endif  // GUARDADVECTRIEMANN 

}

#endif
