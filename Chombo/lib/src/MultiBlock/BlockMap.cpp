#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BlockMap.H"
#include "RealVect.H"
#include "FArrayBox.H"
#include "BoxIterator.H"
#include "UGIO.H"
#include "LoadBalance.H"
#include "CH_HDF5.H"
#include "CellCentersF_F.H"

#include "NamespaceHeader.H"

BlockMap BlockMap::Identity;

// ---------------------------------------------------------
void BlockMap::getCornerCoordinates(FArrayBox&     a_cornerCoords,
                                    const Box&     a_box) const
{
  Box bxNodes = surroundingNodes(a_box);
  FArrayBox mapFab(bxNodes, SpaceDim);
  for (BoxIterator bit(bxNodes); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      RealVect cartesianCoords(iv);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          mapFab(iv, idir) = cartesianCoords[idir];
        }
    }
  setPhysicalFromMap(a_cornerCoords, bxNodes, mapFab);
}


// ---------------------------------------------------------
void BlockMap::getCellCenterCoordinates(FArrayBox&     a_cellCenterCoords,
                                        const Box&     a_box) const
{
  FArrayBox cartesianCenters(a_box, SpaceDim);
  // is unchanged by refinement:  use box indices
  FORT_SETCELLCENTERS(CHF_FRA(cartesianCenters),
                      CHF_BOX(a_box));
  setPhysicalFromMap(a_cellCenterCoords, a_box, cartesianCenters);
}


// ---------------------------------------------------------
void BlockMap::setPhysicalFromMap(/// physical coordinates, SpaceDim components
                                  FArrayBox&         a_physFab,
                                  /// box on which to set physical coordinates
                                  const Box&         a_bx,
                                  /// mapped coordinates, SpaceDim components
                                  const FArrayBox&   a_mapFab) const
{ // Default identity map is unchanged by refinement.
  a_physFab.copy(a_mapFab, a_bx);
}


// ---------------------------------------------------------
void BlockMap::setMapFromPhysical(/// mapped coordinates, SpaceDim components
                                  FArrayBox&         a_mapFab,
                                  /// box on which to set physical coordinates
                                  const Box&         a_bx,
                                  /// physical coordinates, SpaceDim components
                                  const FArrayBox&   a_physFab) const
{ // Default identity map is unchanged by refinement.
  a_mapFab.copy(a_physFab, a_bx);
}


// ---------------------------------------------------------
void BlockMap::setJacobian(/// jacobian on a_bx
                           FArrayBox&         a_jacobianFab,
                           /// box on which to set jacobian
                           const Box&         a_bx,
                           /// mapped coordinates, SpaceDim components
                           const FArrayBox&   a_mapFab) const
{ // Default identity map is unchanged by refinement.
  a_jacobianFab.setVal(0., a_bx, 0);
}


// ---------------------------------------------------------
// This should be overridden by derived class.
bool BlockMap::inThisBlock(const RealVect&   a_physPoint) const
{
  return true;
}


// ---------------------------------------------------------
void BlockMap::cartesianToReal(RealVect& a_point) const
{
  // Going through FArrayBox because we're calling
  // setPhysicalFromMap in a derived class.
  Box box0(IntVect::Zero, IntVect::Zero);
  FArrayBox cartesianFab(box0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      cartesianFab(IntVect::Zero, idir) = a_point[idir];
    }
  FArrayBox realFab(box0, SpaceDim);
  setPhysicalFromMap(realFab, box0, cartesianFab);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_point[idir] = realFab(IntVect::Zero, idir);
    }
}

// ---------------------------------------------------------
void BlockMap::realToCartesian(RealVect& a_point) const
{
  // Going through FArrayBox because we're calling
  // setMapFromPhysical in a derived class.
  Box box0(IntVect::Zero, IntVect::Zero);
  FArrayBox realFab(box0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      realFab(IntVect::Zero, idir) = a_point[idir];
    }
  FArrayBox cartesianFab(box0, SpaceDim);
  setMapFromPhysical(cartesianFab, box0, realFab);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_point[idir] = cartesianFab(IntVect::Zero, idir);
    }
}

// ---------------------------------------------------------
std::string BlockMap::name() const
{
  return "default";
}


// ---------------------------------------------------------
BlockMap* BlockMap::read(std::istream& is) const
{
  //reads nothing
  return new BlockMap();

}

// ---------------------------------------------------------
void BlockMap::write(std::ostream& os) const
{
  // writes nothing
}

// ---------------------------------------------------------
BlockMap* BlockMap::factory(std::istream& is)
{
  char buffer[80];
  memset(buffer, 0, 80);
  is.getline(buffer, 79, ' ');
  std::string n(buffer);
  BlockMap* fac = NULL; // MapConstructors[n];
  if (fac == NULL)
    {
      MayDay::Abort("unrecognized BlockMap type");
    }
  return fac->read(is);
}

#include "NamespaceFooter.H"
