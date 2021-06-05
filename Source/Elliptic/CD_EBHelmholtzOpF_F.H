#ifndef _CD_EBHELMHOLTZOPF_F_H_
#define _CD_EBHELMHOLTZOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDADDBCOTERMTOINVRELCOEF 
#define GUARDADDBCOTERMTOINVRELCOEF 
// Prototype for Fortran procedure addbcotermtoinvrelcoef ...
//
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
      ,CHFp_BOX(cellbox) )
{
 CH_TIMELEAF("FORT_ADDBCOTERMTOINVRELCOEF");
 FORTRAN_NAME( ADDBCOTERMTOINVRELCOEF ,addbcotermtoinvrelcoef )(
      CHFt_FRA1(relcoef)
      ,CHFt_CONST_FRA1(bcoef)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(cellbox) );
}
#endif  // GUARDADDBCOTERMTOINVRELCOEF 

#ifndef GUARDINVERTRELAXATIONCOEFFICIENT 
#define GUARDINVERTRELAXATIONCOEFFICIENT 
// Prototype for Fortran procedure invertrelaxationcoefficient ...
//
void FORTRAN_NAME( INVERTRELAXATIONCOEFFICIENT ,invertrelaxationcoefficient )(
      CHFp_FRA1(relcoef)
      ,CHFp_BOX(cellbox) );

#define FORT_INVERTRELAXATIONCOEFFICIENT FORTRAN_NAME( inlineINVERTRELAXATIONCOEFFICIENT, inlineINVERTRELAXATIONCOEFFICIENT)
#define FORTNT_INVERTRELAXATIONCOEFFICIENT FORTRAN_NAME( INVERTRELAXATIONCOEFFICIENT, invertrelaxationcoefficient)

inline void FORTRAN_NAME(inlineINVERTRELAXATIONCOEFFICIENT, inlineINVERTRELAXATIONCOEFFICIENT)(
      CHFp_FRA1(relcoef)
      ,CHFp_BOX(cellbox) )
{
 CH_TIMELEAF("FORT_INVERTRELAXATIONCOEFFICIENT");
 FORTRAN_NAME( INVERTRELAXATIONCOEFFICIENT ,invertrelaxationcoefficient )(
      CHFt_FRA1(relcoef)
      ,CHFt_BOX(cellbox) );
}
#endif  // GUARDINVERTRELAXATIONCOEFFICIENT 

}

#endif
