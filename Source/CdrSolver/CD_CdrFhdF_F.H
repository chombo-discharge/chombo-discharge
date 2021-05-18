/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_CdrFhdF_F.H
  @brief  Fortran headers for scalar fluctuating hydrodynamics
  @author Robert Marskar
*/

#ifndef CD_CdrFhdF_F_H
#define CD_CdrFhdF_F_H

// Chombo includes
#include <FORT_PROTO.H>
#include <CH_Timer.H>
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>

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

#include <CD_NamespaceFooter.H>

#endif
