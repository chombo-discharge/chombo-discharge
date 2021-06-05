/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzOpF_F.H
  @brief  Fortran headers for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzOpF_F_H
#define CD_EBHelmholtzOpF_F_H

// Chombo includes
#include <FORT_PROTO.H>
#include <REAL.H>

extern "C"
{

#ifndef GUARDADDBCOTERMTOINVRELCOEF 
#define GUARDADDBCOTERMTOINVRELCOEF
  
// Function which adds the diagonal term proportional to bco to the relaxation coefficient
void FORTRAN_NAME( ADDBCOTERMTOINVRELCOEF ,addbcotermtoinvrelcoef )(
      CHFp_FRA1(relcoef)
      ,CHFp_CONST_FRA1(bcoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(cellbox) );

#define FORT_ADDBCOTERMTOINVRELCOEF FORTRAN_NAME( inlineADDBCOTERMTOINVRELCOEF, inlineADDBCOTERMTOINVRELCOEF)
#define FORTNT_ADDBCOTERMTOINVRELCOEF FORTRAN_NAME( ADDBCOTERMTOINVRELCOEF, addbcotermtoinvrelcoef)

inline void FORTRAN_NAME(inlineADDBCOTERMTOINVRELCOEF, inlineADDBCOTERMTOINVRELCOEF)(
      CHFp_FRA1(relcoef)
      ,CHFp_CONST_FRA1(bcoef)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(cellbox))
{
 CH_TIMELEAF("FORT_ADDBCOTERMTOINVRELCOEF");
 FORTRAN_NAME( ADDBCOTERMTOINVRELCOEF ,addbcotermtoinvrelcoef )(
      CHFt_FRA1(relcoef)
      ,CHFt_CONST_FRA1(bcoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(cellbox));
}
#endif  // GUARDADDBCOTERMTOINVRELCOEF 

}

#endif
