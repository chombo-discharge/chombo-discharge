#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// EdgeToCell.cpp
// Dan Martin, Fri, Jan 14, 2000

#include <cassert>

#include "DataIterator.H"
#include "EdgeToCell.H"
#include "EdgeToCellF_F.H"
#include "ProtoInterface.H"
#ifdef USE_PROTO
#include "Proto_Stencil.H"
#endif
#include "NamespaceHeader.H"

#ifdef USE_PROTO

using Proto::Point;
using Proto::BoxData;
using Proto::Var;
using Proto::Stencil;
using Proto::Shift;
using CH_XD::IntVect;
using CH_XD::Box;
using CH_XD::BaseFab;

void
ProtoEdgeToCellPatch(BaseFab<Real>      & a_cellData, 
                     const int          & a_cellComp,
                     const BaseFab<Real>& a_edgeData, 
                     const int          & a_edgeComp, 
                     const Box          & a_cellBox, 
                     const int          & a_idir)
{
  BoxData<double, 1> bdcell, bdedge;
  ProtoCh::aliasBoxData<double, 1>(bdedge, a_edgeData, a_edgeComp);
  ProtoCh::aliasBoxData<double, 1>(bdcell, a_cellData, a_cellComp);
  Proto::Box cellbx = ProtoCh::getProtoBox(a_cellBox);
  Stencil<double> edgeToCellSten = (0.5)*Shift::Basis(a_idir, 1) + (0.5)*Shift::Zeros();
  bdcell |= edgeToCellSten(bdedge, cellbx);
}

void 
maxValLoHi(Var<double,  1>   & a_max,
           Var<double,  1>   & a_loval,
           Var<double,  1>   & a_hival)
{
  double retval =  std::max(a_loval(0), a_hival(0));;
  a_max(0) = retval;
}

/*********************/

void
ProtoEdgeToCellPatchMax(BaseFab<Real>      & a_cellData, 
                        const int          & a_cellComp,
                        const BaseFab<Real>& a_edgeData, 
                        const int          & a_edgeComp, 
                        const Box          & a_cellBox, 
                        const int          & a_idir)
{
  BoxData<double, 1> bdcell, bdedge, bdlo, bdhi;
  ProtoCh::aliasBoxData<double, 1>(bdedge, a_edgeData, a_edgeComp);
  ProtoCh::aliasBoxData<double, 1>(bdcell, a_cellData, a_cellComp);
  Proto::Box cellbx = ProtoCh::getProtoBox(a_cellBox);
  Stencil<double> loSten = (1.0)*Shift::Basis(a_idir, -1);
  Stencil<double> hiSten = (1.0)*Shift::Zeros();
  BoxData<double> lovals(cellbx);
  BoxData<double> hivals(cellbx);
  bdcell = Proto::forall<double>(maxValLoHi, cellbx, lovals, hivals);
}

#endif
// ----------------------------------------------------------------
void EdgeToCell(const LevelData<FluxBox>& a_edgeData,
                LevelData<FArrayBox>& a_cellData)
{
  // this is just a wrapper around the single-grid version
  DataIterator dit = a_edgeData.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      EdgeToCell(a_edgeData[dit()], a_cellData[dit()]);
    }
}

// ----------------------------------------------------------------
void EdgeToCell(const FluxBox& a_edgeData,
                FArrayBox& a_cellData)
{
  // loop over components -- assumption is that in cell-centered
  // data, direction changes faster than component.
  for (int comp = 0; comp<a_edgeData.nComp(); comp++)
    {
      // loop over directions
      for (int dir = 0; dir < SpaceDim; dir++)
        {
          Box cellBox = a_cellData.box();
          cellBox &= a_edgeData.box();
          int cellcomp = SpaceDim*comp + dir;

#ifdef USE_PROTO
          ProtoEdgeToCellPatch(a_cellData, 
                               cellcomp,
                               a_edgeData[dir], 
                               comp, 
                               cellBox, 
                               dir);
#else
          FORT_EDGETOCELL(CHF_CONST_FRA1(a_edgeData[dir],comp),
                          CHF_FRA1(a_cellData, cellcomp),
                          CHF_BOX(cellBox),
                          CHF_CONST_INT(dir));
#endif
        } // end loop over directions
    } // end loop over components
}

void EdgeToCell(const FluxBox& a_edgeData, const int a_edgeComp,
                FArrayBox& a_cellData, const int a_cellComp,
                const int a_dir)
{
   Box cellBox = a_cellData.box();
   cellBox &= a_edgeData.box();

#ifdef USE_PROTO
   ProtoEdgeToCellPatch(a_cellData, 
                        a_cellComp,
                        a_edgeData[a_dir], 
                        a_edgeComp, 
                        cellBox, 
                        a_dir);
#else
   FORT_EDGETOCELL(CHF_CONST_FRA1(a_edgeData[a_dir],a_edgeComp),
                   CHF_FRA1(a_cellData, a_cellComp),
                   CHF_BOX(cellBox),
                   CHF_CONST_INT(a_dir));
#endif
}

void EdgeToCell(const FluxBox& a_edgeData, const int a_edgeComp,
                FArrayBox& a_cellData, const int a_cellComp,
                const Box& a_cellBox, const int a_dir)
{


#ifdef USE_PROTO
   ProtoEdgeToCellPatch(a_cellData, 
                        a_cellComp,
                        a_edgeData[a_dir], 
                        a_edgeComp, 
                        a_cellBox, 
                        a_dir);
#else
  FORT_EDGETOCELL(CHF_CONST_FRA1(a_edgeData[a_dir],a_edgeComp),
                  CHF_FRA1(a_cellData, a_cellComp),
                  CHF_BOX(a_cellBox),
                  CHF_CONST_INT(a_dir));
#endif
}

// max functions

// ----------------------------------------------------------------
void EdgeToCellMax(const LevelData<FluxBox>& a_edgeData,
                   LevelData<FArrayBox>& a_cellData)

{
  // this is just a wrapper around the single-grid version
  DataIterator dit = a_edgeData.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
    {
      EdgeToCellMax(a_edgeData[dit()], a_cellData[dit()]);
    }
}

// ----------------------------------------------------------------
void EdgeToCellMax(const FluxBox& a_edgeData,
                   FArrayBox& a_cellData)
{

  // loop over components -- assumption is that in cell-centered
  // data, direction changes faster than component.
  for (int comp = 0; comp<a_edgeData.nComp(); comp++)
    {
      // loop over directions
      for (int dir = 0; dir < SpaceDim; dir++)
        {
          Box cellBox = a_cellData.box();
          cellBox &= a_edgeData.box();
          int cellcomp = SpaceDim*comp + dir;
#ifdef USE_PROTO
          ProtoEdgeToCellPatchMax(a_cellData, 
                                  cellcomp,
                                  a_edgeData[dir], 
                                  comp,
                                  cellBox, 
                                  dir);
#else
          FORT_EDGETOCELLMAX(CHF_CONST_FRA1(a_edgeData[dir],comp),
                             CHF_FRA1(a_cellData, cellcomp),
                             CHF_BOX(cellBox),
                             CHF_CONST_INT(dir));
#endif
        } // end loop over directions
    } // end loop over components
}

void EdgeToCellMax(const FluxBox& a_edgeData, const int a_edgeComp,
                   FArrayBox& a_cellData, const int a_cellComp,
                   const int a_dir)
{

  const Box& cellBox = a_cellData.box();
#ifdef USE_PROTO
  ProtoEdgeToCellPatchMax(a_cellData, 
                          a_cellComp,
                          a_edgeData[a_dir], 
                          a_edgeComp,
                          cellBox, 
                          a_dir);
#else
  FORT_EDGETOCELLMAX(CHF_CONST_FRA1(a_edgeData[a_dir],a_edgeComp),
                     CHF_FRA1(a_cellData, a_cellComp),
                     CHF_BOX(cellBox),
                     CHF_CONST_INT(a_dir));
#endif
}

void EdgeToCellMax(const FluxBox& a_edgeData, const int a_edgeComp,
                   FArrayBox& a_cellData, const int a_cellComp,
                   const Box& a_cellBox, const int a_dir)
{

#ifdef USE_PROTO
  ProtoEdgeToCellPatchMax(a_cellData, 
                          a_cellComp,
                          a_edgeData[a_dir], 
                          a_edgeComp,
                          a_cellBox, 
                          a_dir);
#else
  FORT_EDGETOCELLMAX(CHF_CONST_FRA1(a_edgeData[a_dir],a_edgeComp),
                     CHF_FRA1(a_cellData, a_cellComp),
                     CHF_BOX(a_cellBox),
                     CHF_CONST_INT(a_dir));
#endif
}

#include "NamespaceFooter.H"
