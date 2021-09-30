#ifndef _CD_CDRPLASMASTEPPERF_F_H_
#define _CD_CDRPLASMASTEPPERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTE_RELAXATION_TIME 
#define GUARDCOMPUTE_RELAXATION_TIME 
// Prototype for Fortran procedure compute_relaxation_time ...
//
void FORTRAN_NAME( COMPUTE_RELAXATION_TIME ,compute_relaxation_time )(
      CHFp_FRA1(relax_time)
      ,CHFp_CONST_FRA1(electric_field)
      ,CHFp_CONST_FRA1(current_density)
      ,CHFp_CONST_REAL(eps0)
      ,CHFp_CONST_REAL(safety)
      ,CHFp_BOX(cell_box) );

#define FORT_COMPUTE_RELAXATION_TIME FORTRAN_NAME( inlineCOMPUTE_RELAXATION_TIME, inlineCOMPUTE_RELAXATION_TIME)
#define FORTNT_COMPUTE_RELAXATION_TIME FORTRAN_NAME( COMPUTE_RELAXATION_TIME, compute_relaxation_time)

inline void FORTRAN_NAME(inlineCOMPUTE_RELAXATION_TIME, inlineCOMPUTE_RELAXATION_TIME)(
      CHFp_FRA1(relax_time)
      ,CHFp_CONST_FRA1(electric_field)
      ,CHFp_CONST_FRA1(current_density)
      ,CHFp_CONST_REAL(eps0)
      ,CHFp_CONST_REAL(safety)
      ,CHFp_BOX(cell_box) )
{
 CH_TIMELEAF("FORT_COMPUTE_RELAXATION_TIME");
 FORTRAN_NAME( COMPUTE_RELAXATION_TIME ,compute_relaxation_time )(
      CHFt_FRA1(relax_time)
      ,CHFt_CONST_FRA1(electric_field)
      ,CHFt_CONST_FRA1(current_density)
      ,CHFt_CONST_REAL(eps0)
      ,CHFt_CONST_REAL(safety)
      ,CHFt_BOX(cell_box) );
}
#endif  // GUARDCOMPUTE_RELAXATION_TIME 

}

#endif
