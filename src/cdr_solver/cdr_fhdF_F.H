#ifndef _CDR_FHDF_F_H_
#define _CDR_FHDF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDHEAVISIDE_MEAN 
#define GUARDHEAVISIDE_MEAN 
// Prototype for Fortran procedure heaviside_mean ...
//
void FORTRAN_NAME( HEAVISIDE_MEAN ,heaviside_mean )(
      CHFp_FRA1(facedata)
      ,CHFp_CONST_FRA1(celldata)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(facebox) );

#define FORT_HEAVISIDE_MEAN FORTRAN_NAME( inlineHEAVISIDE_MEAN, inlineHEAVISIDE_MEAN)
#define FORTNT_HEAVISIDE_MEAN FORTRAN_NAME( HEAVISIDE_MEAN, heaviside_mean)

inline void FORTRAN_NAME(inlineHEAVISIDE_MEAN, inlineHEAVISIDE_MEAN)(
      CHFp_FRA1(facedata)
      ,CHFp_CONST_FRA1(celldata)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(facebox) )
{
 CH_TIMELEAF("FORT_HEAVISIDE_MEAN");
 FORTRAN_NAME( HEAVISIDE_MEAN ,heaviside_mean )(
      CHFt_FRA1(facedata)
      ,CHFt_CONST_FRA1(celldata)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(facebox) );
}
#endif  // GUARDHEAVISIDE_MEAN 

}

#endif
