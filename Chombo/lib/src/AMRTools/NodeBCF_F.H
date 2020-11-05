#ifndef _NODEBCF_F_H_
#define _NODEBCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDFACENODEBC 
#define GUARDFACENODEBC 
// Prototype for Fortran procedure FACENODEBC ...
//
void FORTRAN_NAME( FACENODEBC ,facenodebc )(
      CHFp_FRA(state)
      ,CHFp_CONST_FRA(neumfac)
      ,CHFp_CONST_FRA(dircfac)
      ,CHFp_CONST_FRA(inhmval)
      ,CHFp_BOX(facebox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) );

#define FORT_FACENODEBC FORTRAN_NAME( inlineFACENODEBC, inlineFACENODEBC)
#define FORTNT_FACENODEBC FORTRAN_NAME( FACENODEBC, facenodebc)

inline void FORTRAN_NAME(inlineFACENODEBC, inlineFACENODEBC)(
      CHFp_FRA(state)
      ,CHFp_CONST_FRA(neumfac)
      ,CHFp_CONST_FRA(dircfac)
      ,CHFp_CONST_FRA(inhmval)
      ,CHFp_BOX(facebox)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) )
{
 CH_TIMELEAF("FORT_FACENODEBC");
 FORTRAN_NAME( FACENODEBC ,facenodebc )(
      CHFt_FRA(state)
      ,CHFt_CONST_FRA(neumfac)
      ,CHFt_CONST_FRA(dircfac)
      ,CHFt_CONST_FRA(inhmval)
      ,CHFt_BOX(facebox)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(startcomp)
      ,CHFt_CONST_INT(endcomp) );
}
#endif  // GUARDFACENODEBC 

}

#endif
