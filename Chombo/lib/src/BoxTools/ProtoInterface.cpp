#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ProtoInterface.H"
#include "SPACE.H"
#include "UsingNamespace.H"

#ifdef USE_PROTO

///get point from intvect
Point
ProtoCh::getPoint( const IntVect& a_iv)
{
  Point retval;
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    retval[idir] = a_iv[idir];
  }
  return retval;

}

/// gets proto box from chombo box
Proto::Box  
ProtoCh::getProtoBox( const Box& a_box)
{
  Point ptlo = getPoint(a_box.smallEnd());
  Point pthi = getPoint(a_box.bigEnd());
  return Proto::Box(ptlo, pthi);
}

///get intvect from point
IntVect 
ProtoCh::getIntVect(const  Point  & a_pt)
{
  IntVect retval;
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    retval[idir] = a_pt[idir];
  }
  return retval;
}

///get chombo box from proto box
Box 
ProtoCh::getBox(const Proto::Box & a_bx)
{
  IntVect ivlo = getIntVect(a_bx.low());
  IntVect ivhi = getIntVect(a_bx.high());
  return Box(ivlo, ivhi);
}



#endif


