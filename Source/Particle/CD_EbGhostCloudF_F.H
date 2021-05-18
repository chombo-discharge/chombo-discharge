/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_EbGhostCloudF_F.H
  @brief  Fortran headers for EbGhostCloud
  @author Robert Marskar
*/

#ifndef CD_EbGhostCloudF_F_H
#define CD_EbGhostCloudF_F_H

// Chombo includes
#include <FORT_PROTO.H>
#include <CH_Timer.H>
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>

extern "C"
{

#ifndef GUARDAVERAGE_GHOSTS 
#define GUARDAVERAGE_GHOSTS 
  // Prototype for Fortran procedure average_ghosts ...
  //
  void FORTRAN_NAME( AVERAGE_GHOSTS ,average_ghosts )(
						      CHFp_FRA1(dst)
						      ,CHFp_CONST_FRA1(src)
						      ,CHFp_BOX(finebox)
						      ,CHFp_CONST_INT(refrat)
						      ,CHFp_CONST_REAL(factor) );

#define FORT_AVERAGE_GHOSTS FORTRAN_NAME( inlineAVERAGE_GHOSTS, inlineAVERAGE_GHOSTS)
#define FORTNT_AVERAGE_GHOSTS FORTRAN_NAME( AVERAGE_GHOSTS, average_ghosts)

  inline void FORTRAN_NAME(inlineAVERAGE_GHOSTS, inlineAVERAGE_GHOSTS)(
								       CHFp_FRA1(dst)
								       ,CHFp_CONST_FRA1(src)
								       ,CHFp_BOX(finebox)
								       ,CHFp_CONST_INT(refrat)
								       ,CHFp_CONST_REAL(factor) )
  {
    CH_TIMELEAF("FORT_AVERAGE_GHOSTS");
    FORTRAN_NAME( AVERAGE_GHOSTS ,average_ghosts )(
						   CHFt_FRA1(dst)
						   ,CHFt_CONST_FRA1(src)
						   ,CHFt_BOX(finebox)
						   ,CHFt_CONST_INT(refrat)
						   ,CHFt_CONST_REAL(factor) );
  }
#endif  // GUARDAVERAGE_GHOSTS
  
}

#include <CD_NamespaceFooter.H>

#endif
