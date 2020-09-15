#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _REDUCTIONOPS_H_
#define _REDUCTIONOPS_H_

#include "BoxLayoutData.H"
#include "FArrayBox.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"

void computeFaceReductionWeights(LevelData<FluxBox>& a_weights);
/// Class to do summing operation in the context of CopyTo
/**
   This class will perform a summing operation of the data in src
   in the summingDir direction(s), multiplying by the scale, and placing
   the sum in the corresponding location in dest
*/

class SumOp : public LDOperator<FArrayBox>
{
public:
  Real scale;
  Vector<int> m_summingDir;
  SumOp();

  // single-direction
  SumOp(int a_summingDir);

  // allows for a multidimensional reduction
  SumOp(const Vector<int>& a_summingDir);

  virtual void linearOut(const FArrayBox& arg, void* buf, const Box& R,
                         const Interval& comps) const;

  virtual void linearIn(FArrayBox& arg,  void* buf, const Box& R,
                        const Interval& comps) const;
  void op(FArrayBox& dest,
          const Box& RegionFrom,
          const Interval& Cdest,
          const Box& RegionTo,
          const FArrayBox& src,
          const Interval& Csrc) const;

};

// (DFM 11/13/08) as currently implemented, FaceSumOp doesn't
// do the right thing for multiple adjoining grids -- it will
// double-count overlying faces where boxes adjoin. Since we don't
// actually need the face-centered summing operator at the moment,
// take the cowardly path of just commenting it out to revisit if it
// becomes a needed member of the Chombo family


class FaceSumOp: public LDOperator<FluxBox>
{
public:

  FaceSumOp();
  FaceSumOp(const int& a_summingDir);
  
  virtual ~FaceSumOp()
  {
  }
  
  void setDir  (const int& a_sumDir);
  void setScale(const Real& a_scale );
  
  void op(FluxBox       & dataTo     ,
          const Box     & RegionFrom ,
          const Interval& Cdest      ,
          const Box     & RegionTo   ,
          const FluxBox & src        ,
          const Interval& Csrc       ) const;



  virtual void linearIn(FluxBox       & a_arg  ,  
                        void*           a_buf  , 
                        const Box     & a_R    ,
                        const Interval& a_comps) const;


  virtual void linearOut(const FluxBox & a_arg  ,         
                         void*           a_buf  , 
                         const Box     & a_R    ,
                         const Interval& a_comps) const;
  
  int size(const FluxBox & a_fluxBox,
           const Box     & a_bx     , 
           const Interval& a_comps  ) const;


  Real m_scale     ;
  int  m_summingDir;
  
};

/// Class to do spreading operation in the context of CopyTo
/**
   This class will perform a spreading operation of the data in src
   along the summingDir direction(s), multiplying by the scale, and placing
   the src values in the corresponding locations in dest.
*/

class SpreadingOp : public LDOperator<FArrayBox>
{
public:
  Real scale;
  Vector<int> m_spreadingDir;
  SpreadingOp();

  // single spreading direction
  SpreadingOp(int a_spreadingDir);

  SpreadingOp(const Vector<int>& a_spreadingDir);

  virtual void linearIn(FArrayBox& arg,  void* buf, const Box& R,
                        const Interval& comps) const;
  void op(FArrayBox& dest,
          const Box& RegionFrom,
          const Interval& Cdest,
          const Box& RegionTo,
          const FArrayBox& src,
          const Interval& Csrc) const;

  void applyOp(FArrayBox& dest,
               const Box& RegionFrom,
               const Interval& Cdest,
               const Box& RegionTo,
               const FArrayBox& src,
               const Interval& Csrc,
               Real a_scale) const;

};

/// Class to do spreading operation in the context of CopyTo
/**
   This class will perform a spreading operation of the data in src
   along the summingDir direction(s), multiplying by the scale, and placing
   the src values in the corresponding locations in dest.
*/

class FaceSpreadingOp : public LDOperator<FluxBox>
{
public:
  Real scale;
  Vector<int> m_spreadingDir;
  FaceSpreadingOp();

  // single spreading direction
  FaceSpreadingOp(int a_spreadingDir);

  // allows for multiple directions
  FaceSpreadingOp(const Vector<int>& a_spreadingDir);

  virtual void linearIn(FluxBox& arg,  void* buf, const Box& R,
                        const Interval& comps) const;
  void op(FluxBox& dest,
          const Box& RegionFrom,
          const Interval& Cdest,
          const Box& RegionTo,
          const FluxBox& src,
          const Interval& Csrc) const;

  void applyOp(FluxBox& dest,
               const Box& RegionFrom,
               const Interval& Cdest,
               const Box& RegionTo,
               const FluxBox& src,
               const Interval& Csrc,
               Real a_scale) const;

};

#include "NamespaceFooter.H"

#endif
