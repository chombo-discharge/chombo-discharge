#ifndef _DIRICHLETPOISSONDOMAINBCF_F_H_
#define _DIRICHLETPOISSONDOMAINBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETDIRICHLETFACEFLUX 
#define GUARDSETDIRICHLETFACEFLUX 
// Prototype for Fortran procedure SETDIRICHLETFACEFLUX ...
//
void FORTRAN_NAME( SETDIRICHLETFACEFLUX ,setdirichletfaceflux )(
      CHFp_FRA1(faceFlux)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(value)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside)
      ,CHFp_BOX(box) );

#define FORT_SETDIRICHLETFACEFLUX FORTRAN_NAME( inlineSETDIRICHLETFACEFLUX, inlineSETDIRICHLETFACEFLUX)
#define FORTNT_SETDIRICHLETFACEFLUX FORTRAN_NAME( SETDIRICHLETFACEFLUX, setdirichletfaceflux)

inline void FORTRAN_NAME(inlineSETDIRICHLETFACEFLUX, inlineSETDIRICHLETFACEFLUX)(
      CHFp_FRA1(faceFlux)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(value)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SETDIRICHLETFACEFLUX");
 FORTRAN_NAME( SETDIRICHLETFACEFLUX ,setdirichletfaceflux )(
      CHFt_FRA1(faceFlux)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(value)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iside)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSETDIRICHLETFACEFLUX 

#ifndef GUARDSETHODIRICHLETFACEFLUX 
#define GUARDSETHODIRICHLETFACEFLUX 
// Prototype for Fortran procedure SETHODIRICHLETFACEFLUX ...
//
void FORTRAN_NAME( SETHODIRICHLETFACEFLUX ,sethodirichletfaceflux )(
      CHFp_FRA1(faceFlux)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(value)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside)
      ,CHFp_BOX(box) );

#define FORT_SETHODIRICHLETFACEFLUX FORTRAN_NAME( inlineSETHODIRICHLETFACEFLUX, inlineSETHODIRICHLETFACEFLUX)
#define FORTNT_SETHODIRICHLETFACEFLUX FORTRAN_NAME( SETHODIRICHLETFACEFLUX, sethodirichletfaceflux)

inline void FORTRAN_NAME(inlineSETHODIRICHLETFACEFLUX, inlineSETHODIRICHLETFACEFLUX)(
      CHFp_FRA1(faceFlux)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(value)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iside)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SETHODIRICHLETFACEFLUX");
 FORTRAN_NAME( SETHODIRICHLETFACEFLUX ,sethodirichletfaceflux )(
      CHFt_FRA1(faceFlux)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(value)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iside)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSETHODIRICHLETFACEFLUX 

}

#endif
