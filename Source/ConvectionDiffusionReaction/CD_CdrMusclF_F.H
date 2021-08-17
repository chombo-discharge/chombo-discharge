#ifndef _CD_CDRMUSCLF_F_H_
#define _CD_CDRMUSCLF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMUSCL_SLOPES 
#define GUARDMUSCL_SLOPES 
// Prototype for Fortran procedure muscl_slopes ...
//
void FORTRAN_NAME( MUSCL_SLOPES ,muscl_slopes )(
      CHFp_FRA1(deltaC)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(region) );

#define FORT_MUSCL_SLOPES FORTRAN_NAME( inlineMUSCL_SLOPES, inlineMUSCL_SLOPES)
#define FORTNT_MUSCL_SLOPES FORTRAN_NAME( MUSCL_SLOPES, muscl_slopes)

inline void FORTRAN_NAME(inlineMUSCL_SLOPES, inlineMUSCL_SLOPES)(
      CHFp_FRA1(deltaC)
      ,CHFp_CONST_FRA1(state)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_MUSCL_SLOPES");
 FORTRAN_NAME( MUSCL_SLOPES ,muscl_slopes )(
      CHFt_FRA1(deltaC)
      ,CHFt_CONST_FRA1(state)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(region) );
}
#endif  // GUARDMUSCL_SLOPES 

#ifndef GUARDMUSCL_UPWIND 
#define GUARDMUSCL_UPWIND 
// Prototype for Fortran procedure muscl_upwind ...
//
void FORTRAN_NAME( MUSCL_UPWIND ,muscl_upwind )(
      CHFp_FRA1(face_states)
      ,CHFp_CONST_FRA1(cell_slopes)
      ,CHFp_CONST_FRA1(cell_states)
      ,CHFp_CONST_FRA1(face_vel)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(face_box) );

#define FORT_MUSCL_UPWIND FORTRAN_NAME( inlineMUSCL_UPWIND, inlineMUSCL_UPWIND)
#define FORTNT_MUSCL_UPWIND FORTRAN_NAME( MUSCL_UPWIND, muscl_upwind)

inline void FORTRAN_NAME(inlineMUSCL_UPWIND, inlineMUSCL_UPWIND)(
      CHFp_FRA1(face_states)
      ,CHFp_CONST_FRA1(cell_slopes)
      ,CHFp_CONST_FRA1(cell_states)
      ,CHFp_CONST_FRA1(face_vel)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(face_box) )
{
 CH_TIMELEAF("FORT_MUSCL_UPWIND");
 FORTRAN_NAME( MUSCL_UPWIND ,muscl_upwind )(
      CHFt_FRA1(face_states)
      ,CHFt_CONST_FRA1(cell_slopes)
      ,CHFt_CONST_FRA1(cell_states)
      ,CHFt_CONST_FRA1(face_vel)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(face_box) );
}
#endif  // GUARDMUSCL_UPWIND 

}

#endif
