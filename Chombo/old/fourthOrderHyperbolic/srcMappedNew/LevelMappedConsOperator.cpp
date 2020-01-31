#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>

#include "LevelMappedConsOperator.H"
#include "PatchMappedConsOperator.H"
#include "FourthOrderUtil.H" // functions with names beginning "fourthOrder"
#include "CellToEdge.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
LevelMappedConsOperator::LevelMappedConsOperator()
{
  m_defined = false;
  m_dx = 0.0;
  m_refineCoarse = 0;
  m_patchConsOperatorPtr = new PatchMappedConsOperator();
  m_coordSysPtr = NULL;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
LevelMappedConsOperator::~LevelMappedConsOperator()
{
}


//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::setCoordSys(NewFourthOrderCoordSys* a_coordSysPtr)
{
  m_coordSysPtr = a_coordSysPtr;
  PatchMappedConsOperator& patchConsOperator =
    (PatchMappedConsOperator&) *m_patchConsOperatorPtr;
  patchConsOperator.setCoordSys(a_coordSysPtr);
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::computeFaceAverages(
                                       LevelData<FluxBox>& a_face_data,
                                       const LevelData<FArrayBox>& a_cell_data ) const
{
  if (m_spaceOrder == 4)
    {
      fourthOrderCellToFace( a_face_data, a_cell_data );
    }
  else if (m_spaceOrder == 2)
    {
      CellToEdge( a_cell_data, a_face_data );
    }
  else
    {
      MayDay::Error("Bad Space Order in LevelMappedConsOperator");
    }
}

//////////////////////////////////////////////////////////////////////////////

inline void
LevelMappedConsOperator::getPhysicalCellVolumes(
                                                LevelData<FArrayBox>& a_volumes ) const
{
  // Do I need this?
  // a_volumes.define( m_coordSysPtr->getCellVolumes() );
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::computeCompFaceFluxes( LevelData<FluxBox>& a_uTimesV,
                                          const LevelData<FluxBox>& a_u,
                                          const LevelData<FluxBox>& a_v) const
{
   // Compute the SpaceDim-by-nComp face-averaged fluxes in computational
   // space, where a_v is the SpaceDim-dimensional velocity vector and
   // a_u is the nComp-dim state vector
   int ncomp = a_u.nComp();
   CH_assert(a_v.nComp() == SpaceDim);
   CH_assert(a_uTimesV.nComp() == SpaceDim * ncomp);
   PatchMappedConsOperator& patchConsOperator =
     (PatchMappedConsOperator&) *m_patchConsOperatorPtr;
   // loop over boxes
   DataIterator dit = a_u.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
     {
       FluxBox& thisUV = a_uTimesV[dit];
       const FluxBox& thisU = a_u[dit];
       const FluxBox& thisV = a_v[dit];
       patchConsOperator.computeCompFaceFluxes(thisUV, thisU, thisV);
     } // end loop over boxes
}

//////////////////////////////////////////////////////////////////////////////
void LevelMappedConsOperator::cellUJToCellU(LevelData<FArrayBox>& a_Uavg,
                                            const LevelData<FArrayBox>& a_UJavg) const
{
  PatchMappedConsOperator& patchConsOperator =
    (PatchMappedConsOperator&) *m_patchConsOperatorPtr;
  // loop over boxes
  DataIterator dit = a_UJavg.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      setPatchIndex(dit());

      FArrayBox& thisU = a_Uavg[dit];
      const FArrayBox& thisUJ = a_UJavg[dit];
      patchConsOperator.cellUJToCellU(thisU, thisUJ);
    } // end loop over boxes
}

#include "NamespaceFooter.H"
