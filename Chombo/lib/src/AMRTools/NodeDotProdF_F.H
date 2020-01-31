#ifndef _NODEDOTPRODF_F_H_
#define _NODEDOTPRODF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDTRAPWEIGHTS 
#define GUARDTRAPWEIGHTS 
// Prototype for Fortran procedure TRAPWEIGHTS ...
//
void FORTRAN_NAME( TRAPWEIGHTS ,trapweights )(
      CHFp_FRA1(wt)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dv) );

#define FORT_TRAPWEIGHTS FORTRAN_NAME( inlineTRAPWEIGHTS, inlineTRAPWEIGHTS)
#define FORTNT_TRAPWEIGHTS FORTRAN_NAME( TRAPWEIGHTS, trapweights)

inline void FORTRAN_NAME(inlineTRAPWEIGHTS, inlineTRAPWEIGHTS)(
      CHFp_FRA1(wt)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dv) )
{
 CH_TIMELEAF("FORT_TRAPWEIGHTS");
 FORTRAN_NAME( TRAPWEIGHTS ,trapweights )(
      CHFt_FRA1(wt)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dv) );
}
#endif  // GUARDTRAPWEIGHTS 

#ifndef GUARDNODEDOTPRODUCT 
#define GUARDNODEDOTPRODUCT 
// Prototype for Fortran procedure NODEDOTPRODUCT ...
//
void FORTRAN_NAME( NODEDOTPRODUCT ,nodedotproduct )(
      CHFp_REAL(dotprodout)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_CONST_FRA(bfab)
      ,CHFp_CONST_FRA1(wfab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) );

#define FORT_NODEDOTPRODUCT FORTRAN_NAME( inlineNODEDOTPRODUCT, inlineNODEDOTPRODUCT)
#define FORTNT_NODEDOTPRODUCT FORTRAN_NAME( NODEDOTPRODUCT, nodedotproduct)

inline void FORTRAN_NAME(inlineNODEDOTPRODUCT, inlineNODEDOTPRODUCT)(
      CHFp_REAL(dotprodout)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_CONST_FRA(bfab)
      ,CHFp_CONST_FRA1(wfab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) )
{
 CH_TIMELEAF("FORT_NODEDOTPRODUCT");
 FORTRAN_NAME( NODEDOTPRODUCT ,nodedotproduct )(
      CHFt_REAL(dotprodout)
      ,CHFt_CONST_FRA(afab)
      ,CHFt_CONST_FRA(bfab)
      ,CHFt_CONST_FRA1(wfab)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(startcomp)
      ,CHFt_CONST_INT(endcomp) );
}
#endif  // GUARDNODEDOTPRODUCT 

#ifndef GUARDNODENORM 
#define GUARDNODENORM 
// Prototype for Fortran procedure NODENORM ...
//
void FORTRAN_NAME( NODENORM ,nodenorm )(
      CHFp_REAL(normout)
      ,CHFp_CONST_INT(p)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_CONST_FRA1(wfab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) );

#define FORT_NODENORM FORTRAN_NAME( inlineNODENORM, inlineNODENORM)
#define FORTNT_NODENORM FORTRAN_NAME( NODENORM, nodenorm)

inline void FORTRAN_NAME(inlineNODENORM, inlineNODENORM)(
      CHFp_REAL(normout)
      ,CHFp_CONST_INT(p)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_CONST_FRA1(wfab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) )
{
 CH_TIMELEAF("FORT_NODENORM");
 FORTRAN_NAME( NODENORM ,nodenorm )(
      CHFt_REAL(normout)
      ,CHFt_CONST_INT(p)
      ,CHFt_CONST_FRA(afab)
      ,CHFt_CONST_FRA1(wfab)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(startcomp)
      ,CHFt_CONST_INT(endcomp) );
}
#endif  // GUARDNODENORM 

#ifndef GUARDNODEINTEGRAL 
#define GUARDNODEINTEGRAL 
// Prototype for Fortran procedure NODEINTEGRAL ...
//
void FORTRAN_NAME( NODEINTEGRAL ,nodeintegral )(
      CHFp_REAL(ans)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_CONST_FRA1(wfab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) );

#define FORT_NODEINTEGRAL FORTRAN_NAME( inlineNODEINTEGRAL, inlineNODEINTEGRAL)
#define FORTNT_NODEINTEGRAL FORTRAN_NAME( NODEINTEGRAL, nodeintegral)

inline void FORTRAN_NAME(inlineNODEINTEGRAL, inlineNODEINTEGRAL)(
      CHFp_REAL(ans)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_CONST_FRA1(wfab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) )
{
 CH_TIMELEAF("FORT_NODEINTEGRAL");
 FORTRAN_NAME( NODEINTEGRAL ,nodeintegral )(
      CHFt_REAL(ans)
      ,CHFt_CONST_FRA(afab)
      ,CHFt_CONST_FRA1(wfab)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(startcomp)
      ,CHFt_CONST_INT(endcomp) );
}
#endif  // GUARDNODEINTEGRAL 

#ifndef GUARDNODEMAXNORM 
#define GUARDNODEMAXNORM 
// Prototype for Fortran procedure NODEMAXNORM ...
//
void FORTRAN_NAME( NODEMAXNORM ,nodemaxnorm )(
      CHFp_REAL(normout)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) );

#define FORT_NODEMAXNORM FORTRAN_NAME( inlineNODEMAXNORM, inlineNODEMAXNORM)
#define FORTNT_NODEMAXNORM FORTRAN_NAME( NODEMAXNORM, nodemaxnorm)

inline void FORTRAN_NAME(inlineNODEMAXNORM, inlineNODEMAXNORM)(
      CHFp_REAL(normout)
      ,CHFp_CONST_FRA(afab)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(startcomp)
      ,CHFp_CONST_INT(endcomp) )
{
 CH_TIMELEAF("FORT_NODEMAXNORM");
 FORTRAN_NAME( NODEMAXNORM ,nodemaxnorm )(
      CHFt_REAL(normout)
      ,CHFt_CONST_FRA(afab)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(startcomp)
      ,CHFt_CONST_INT(endcomp) );
}
#endif  // GUARDNODEMAXNORM 

}

#endif
