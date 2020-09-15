#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  ANAG, LBNL

#include <cmath>

#include "EBData.H"
#include "VoFIterator.H"
#include "FaceIterator.H"
#include "EBISBox.H"
#include "PolyGeom.H"
#include "EBArith.H"
#include "NamespaceHeader.H"
VolIndex f_debugVoF2(IntVect(D_DECL(1,1,1)), 0);
bool EBDataImplem::s_verbose = false;

/************************/
void calculateEBNormalMoments(IndMomSpaceDim       &   a_EBNormalMomentsDir,
                              const IndMomSpaceDim &   a_partialDerOfNormal,
                              const IndMomSpaceDim &   a_EBMoments)
{
  CH_TIME("calculateEBNormalMoments");
  a_EBNormalMomentsDir.setToZero();
  //here we are computing sum(m_p) = integral(ni(x-xbar)^p1/p! (x-xbar)^p2)
  // = sum_p2 ((n^p2)/p2!) m_p+p12
  for(MomItSpaceDim momitOuter; momitOuter.ok(); ++momitOuter)
    {
      //this is where the answer goes 
      IndexTM<int,SpaceDim>  p = momitOuter();

      for(MomItSpaceDim momitInner; momitInner.ok(); ++momitInner)
        {
          IndexTM<int,SpaceDim> p2 = momitInner();
          IndexTM<int,SpaceDim> pplusp2 = p + p2;
          if(pplusp2.sum() <= CH_EBIS_ORDER)
            {
              Real p2fact = pfactorial(p2);
              Real incr = (a_partialDerOfNormal[p2]/p2fact)*a_EBMoments[pplusp2];
              a_EBNormalMomentsDir[p] += incr;
            }
        }
    }
}
/*******************************/

bool 
EBDataImplem::
irregFace(const FaceIndex& a_face) const
{
  bool retval = (m_volData.getIVS().contains(a_face.gridIndex(Side::Lo)) || m_volData.getIVS().contains(a_face.gridIndex(Side::Hi)));
  return retval;
}
/*******************************/
bool 
EBDataImplem::
irregVoF(const VolIndex& a_vof) const
{
  bool retval = (m_volData.getIVS().contains(a_vof.gridIndex()));
  return retval;
}

/*******************************/
IndMomSDMinOne 
EBDataImplem::
getFaceMoments(const FaceIndex& a_face) const
{
  //if this fails, this was built using a geometry generation object 
  //that does not support moments (try using WrappedGShop instead of GeometryShop)
  if(!m_hasMoments)
  {
    MayDay::Error("geometry generation used did not support higher order moments");
  }

  IndMomSDMinOne retval;
  if(irregFace(a_face))
  {
    retval = m_faceData[a_face.direction()](a_face, 0).m_faceMoments;
  }
  else
  {
    retval= m_regularAreaMoments;
  }
  return retval;
}
/*******************************/
IndMomSpaceDim 
EBDataImplem::
getVolumeMoments(const VolIndex& a_vof) const
{
  //if this fails, this was built using a geometry generation object 
  //that does not support moments (try using WrappedGShop instead of GeometryShop)
  if(!m_hasMoments)
  {
    MayDay::Error("geometry generation used did not support higher order moments");
  }
  IndMomSpaceDim retval;
  bool isIrreg = irregVoF(a_vof);
  if(isIrreg)
  {
    retval = m_volData(a_vof, 0).m_volumeMoments;
  }
  else
  {
    retval = m_regularVolumeMoments;
  }
  return retval;
}

/*******************************/
IndMomSpaceDim 
EBDataImplem::
getEBMoments(const VolIndex& a_vof) const
{
  //if this fails, this was built using a geometry generation object 
  //that does not support moments (try using WrappedGShop instead of GeometryShop)
  if(!m_hasMoments)
  {
    MayDay::Error("the geometry generation tool you chose did not support higher order moments");
  }
  IndMomSpaceDim retval;
  if(irregVoF(a_vof))
  {
    retval = m_volData(a_vof, 0).m_averageFace.m_EBMoments;
  }
  else
  {
    retval.setToZero();
  }
  return retval;
}
/*******************************/
IndMomSpaceDim 
EBDataImplem::
getEBNormalPartialDerivs(const VolIndex& a_vof, int a_normalComponent) const
{
  //if this fails, this was built using a geometry generation object 
  //that does not support moments (try using WrappedGShop instead of GeometryShop)
  if(!m_hasMoments)
  {
    MayDay::Error("the geometry generation tool you chose did not support higher order moments");
  }
  CH_assert(a_normalComponent >= 0);
  CH_assert(a_normalComponent < SpaceDim);
  
  IndMomSpaceDim retval;
  if(irregVoF(a_vof))
  {
    retval = m_volData(a_vof, 0).m_averageFace.m_normalPartialDeriv[a_normalComponent];
  }
  else
  {
    retval.setToZero();
  }
  return retval;
}
/*******************************/
IndMomSpaceDim 
EBDataImplem::
getEBNormalMoments(const VolIndex& a_vof, int a_normalComponent) const
{
  CH_assert(a_normalComponent >= 0);
  CH_assert(a_normalComponent < SpaceDim);
  //if this fails, this was built using a geometry generation object 
  //that does not support moments (try using WrappedGShop instead of GeometryShop)
  if(!m_hasMoments)
  {
    MayDay::Error("the geometry generation tool you chose did not support higher order moments");
  }
  
  IndMomSpaceDim retval;
  if(irregVoF(a_vof))
    {
      retval = m_volData(a_vof, 0).m_averageFace.m_EBNormalMoments[a_normalComponent];
    }
  else
    {
      retval.setToZero();
    }
  return retval;
}
/*******************************/
BoundaryData::BoundaryData()
  :
  m_bndryArea(0.0),
  m_normal(RealVect::Zero),
  m_bndryCentroid(RealVect::Zero),
  m_bndryPhase(-1)
{
  m_volIndex = VolIndex(IntVect::Zero, 0);
}

/*******************************/
VolData::VolData()
  : m_volFrac(1.0),
  m_volCentroid(RealVect::Zero)
{
}

/*******************************/
//BoundaryData specialization of linearSize
template < >
int linearSize(const BoundaryData& a_data)
{

  int retval = 0;
  retval += sizeof(Real); //bndryarea
  retval += linearSize(a_data.m_normal);
  retval += linearSize(a_data.m_bndryCentroid);
  retval += sizeof(int) ;     //bndryPhase
  retval += a_data.m_volIndex.linearSize(); //a_data.m_volIndex
  retval += IndMomSpaceDim::linearSize(); //eb moments
  retval += SpaceDim*IndMomSpaceDim::linearSize(); //normalderivs
  retval += SpaceDim*IndMomSpaceDim::linearSize(); //ebnormalmoments
  return retval;
}


/*******************************/
template < >
void linearIn(BoundaryData& a_data, const void* const a_voidBuf)
{
  unsigned char* buf = (unsigned char*) a_voidBuf;
  size_t incr = 0;
  size_t sze  = 0;
  
  Real* rbuf = (Real*) buf;
  a_data.m_bndryArea = *rbuf;
  incr = sizeof(Real);
  buf += incr;
  sze += incr;

  linearIn(a_data.m_normal, buf);
  incr = linearSize(a_data.m_normal);
  buf += incr;
  sze += incr;

  linearIn(a_data.m_bndryCentroid, buf);
  incr = linearSize(a_data.m_bndryCentroid);
  buf += incr;
  sze += incr;

  int* intbuf = (int*) buf;
  a_data.m_bndryPhase = *intbuf;
  incr = sizeof(int);
  buf += incr;
  sze += incr;

  linearIn(a_data.m_volIndex, buf);
  incr = linearSize(a_data.m_volIndex);
  buf += incr;
  sze += incr;

  linearIn(a_data.m_EBMoments, buf);
  incr = linearSize(a_data.m_EBMoments);
  buf += incr;
  sze += incr;

  for(int idir = 0; idir < SpaceDim; idir++)
  {
    linearIn(a_data.m_normalPartialDeriv[idir], buf);
    incr = linearSize(a_data.m_normalPartialDeriv[idir]);
    buf += incr;
    sze += incr;

    linearIn(a_data.m_EBNormalMoments[idir], buf);
    incr = linearSize(a_data.m_EBNormalMoments[idir]);
    buf += incr;
    sze += incr;
  }
}
/*******************************/
//BoundaryData specialization of linearOut
template < >
void linearOut(void* const a_voidBuf, const BoundaryData& a_data)
{
  unsigned char* buf = (unsigned char*) a_voidBuf;
  size_t incr = 0;
  size_t sze  = 0;
  
  Real* rbuf = (Real*) buf;
  *rbuf = a_data.m_bndryArea;
  incr = sizeof(Real);
  buf += incr;
  sze += incr;

  linearOut(buf, a_data.m_normal);
  incr = linearSize(a_data.m_normal);
  buf += incr;
  sze += incr;

  linearOut(buf, a_data.m_bndryCentroid);
  incr = linearSize(a_data.m_bndryCentroid);
  buf += incr;
  sze += incr;

  int* intbuf = (int*) buf;
  *intbuf = a_data.m_bndryPhase;
  incr = sizeof(int);
  buf += incr;
  sze += incr;

  if(!a_data.m_volIndex.isDefined())
  {
    MayDay::Error("about to linearout undefined volindex");
  }
  linearOut(buf, a_data.m_volIndex);
  incr = linearSize(a_data.m_volIndex);
  buf += incr;
  sze += incr;

  linearOut(buf, a_data.m_EBMoments);
  incr = linearSize(a_data.m_EBMoments);
  buf += incr;
  sze += incr;

  for(int idir = 0; idir < SpaceDim; idir++)
  {
    linearOut(buf, a_data.m_normalPartialDeriv[idir]);
    incr = linearSize(a_data.m_normalPartialDeriv[idir]);
    buf += incr;
    sze += incr;

    linearOut(buf, a_data.m_EBNormalMoments[idir]);
    incr = linearSize(a_data.m_EBNormalMoments[idir]);
    buf += incr;
    sze += incr;
  }
}
/*******************************/
//VolData specialization of linearSize
template < >
int linearSize(const VolData& a_data)
{
  int retval = 0;
  retval += sizeof(Real); //volfrac
  retval += linearSize(    a_data.m_volCentroid);
  retval += linearSize(    a_data.m_averageFace);
  retval += linearSize(    a_data.m_volumeMoments);
  retval += linearListSize(a_data.m_phaseFaces);
  return retval;
}

//VolData specialization of linearIn
template < >
void linearIn(VolData& a_data, const void* const a_voidBuf)
{

  unsigned char* buf = (unsigned char*) a_voidBuf;
  size_t incr = 0;
  size_t sze  = 0;

  //get the vector size
  Real* rbuf =  (Real*) buf;
  a_data.m_volFrac = *rbuf;
  incr = sizeof(Real);
  buf += incr;
  sze += incr;
  
  linearIn(a_data.m_volCentroid, buf);
  incr = linearSize(a_data.m_volCentroid);
  buf += incr;
  sze += incr;

  linearIn(a_data.m_averageFace, buf);
  incr = linearSize(a_data.m_averageFace);
  buf += incr;
  sze += incr;

  linearIn(a_data.m_volumeMoments, buf);
  incr = linearSize(a_data.m_volumeMoments);
  buf += incr;
  sze += incr;

  linearListIn(a_data.m_phaseFaces, buf);
  incr = linearListSize(a_data.m_phaseFaces);
  buf += incr;
  sze += incr;
}
//VolData specialization of linearOut
template < >
void linearOut(void* const a_voidBuf, const VolData& a_data)
{
  unsigned char* buf = (unsigned char*) a_voidBuf;
  size_t incr = 0;
  size_t sze  = 0;

  //get the vector size
  Real* rbuf =  (Real*) buf;
  *rbuf = a_data.m_volFrac;
  incr = sizeof(Real);
  buf += incr;
  sze += incr;
  
  linearOut(buf, a_data.m_volCentroid);
  incr = linearSize(a_data.m_volCentroid);
  buf += incr;
  sze += incr;

  linearOut(buf, a_data.m_averageFace);
  incr = linearSize(a_data.m_averageFace);
  buf += incr;
  sze += incr;

  linearOut(buf, a_data.m_volumeMoments);
  incr = linearSize(a_data.m_volumeMoments);
  buf += incr;
  sze += incr;

  linearListOut(buf, a_data.m_phaseFaces);
  incr = linearListSize(a_data.m_phaseFaces);
  buf += incr;
  sze += incr;
}
/*****/
template < >
int linearSize(const FaceData& a_data)
{
  int retval = 0;
  retval += sizeof(Real); //areaFrac
  retval += linearSize(a_data.m_faceCentroid);
  retval += linearSize(a_data.m_faceMoments);
  return retval;
}
/*****/
template < >
void linearIn(FaceData& a_data, const void* const a_voidBuf)
{
  unsigned char* buf = (unsigned char*) a_voidBuf;
  size_t incr = 0;
  size_t sze  = 0;

  Real* rbuf = (Real*) buf;
  a_data.m_areaFrac = *rbuf;
  incr = sizeof(Real);
  buf += incr;
  sze += incr;

  linearIn(a_data.m_faceCentroid, buf);
  incr = linearSize(a_data.m_faceCentroid);
  buf += incr;
  sze += incr;

  linearIn(a_data.m_faceMoments, buf);
  incr = linearSize(a_data.m_faceMoments);
  buf += incr;
  sze += incr;
}
/*****/
//VolData specialization of linearOut
template < >
void linearOut(void* const a_voidBuf, const FaceData& a_data)
{
  unsigned char* buf = (unsigned char*) a_voidBuf;
  size_t incr = 0;
  size_t sze  = 0;

  Real* rbuf = (Real*) buf;
  *rbuf = a_data.m_areaFrac;
  incr = sizeof(Real);
  buf += incr;
  sze += incr;

  linearOut(buf, a_data.m_faceCentroid);
  incr = linearSize(a_data.m_faceCentroid);
  buf += incr;
  sze += incr;

  linearOut(buf, a_data.m_faceMoments);
  incr = linearSize(a_data.m_faceMoments);
  buf += incr;
  sze += incr;
}

/************************/
void EBDataImplem::setVerbose(bool a_verbose)
{
  s_verbose = a_verbose;
}

bool EBDataImplem::s_verboseDebug = false;
/************************/
void EBDataImplem::setVerboseDebug(bool a_verboseDebug)
{
  s_verboseDebug = a_verboseDebug;
}
/************************/
void EBDataImplem::
addFullIrregularVoFs(const IntVectSet& a_vofsToChange,
                     const EBGraph&    a_newGraph,
                     const BaseIVFAB<VolData>&     a_grownData,
                     const EBGraph&    a_oldGraph)
{
  if (!a_vofsToChange.isEmpty())
    {
      //calculate set by adding in new intvects
      const IntVectSet& ivsOld = m_volData.getIVS();
      IntVectSet ivsNew = ivsOld | a_vofsToChange;

      //for copying purposes
      Box minBox = ivsOld.minBox();
      Interval interv(0,0);

      //save old data into temporarys
      BaseIVFAB<VolData>  oldVolData(ivsOld, a_oldGraph, 1);
      oldVolData.copy(minBox, interv, minBox, m_volData, interv);
      BaseIFFAB<FaceData> oldFaceData[SpaceDim];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          oldFaceData[idir].define(ivsOld, a_oldGraph, idir, 1);
          oldFaceData[idir].copy(minBox, interv, minBox, m_faceData[idir], interv);
        }

      //redefine member data with new IVS and copy old data back
      //into it.
      m_volData.define(ivsNew, a_newGraph, 1);
      m_volData.copy(minBox, interv, minBox, oldVolData, interv);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_faceData[idir].define(ivsNew, a_newGraph, idir, 1);
          m_faceData[idir].copy(minBox, interv, minBox, oldFaceData[idir], interv);
        }

      //now put correct data into new vofs.  the volume fraction of a formally
      //regular cell is always unity
      VolData  fullVol;
      fullVol.setToRegular(m_dx);

      for (VoFIterator vofit(a_vofsToChange, a_newGraph); vofit.ok(); ++vofit)
        {
          m_volData(vofit(), 0) = fullVol;
        }

      //there are rare cases that arise from coarsening where the area fractions will not be one.
      //we can check this by seeing the number of faces on the side
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          FaceIterator faceit(a_vofsToChange, a_newGraph, idir, FaceStop::SurroundingWithBoundary);
          for (faceit.reset(); faceit.ok(); ++faceit)
            {

              FaceData fullFace;
              fullFace.setToRegular(m_dx);
              if (!faceit().isBoundary())
                {
                  IntVect ivhi = faceit().gridIndex(Side::Hi);
                  IntVect ivlo = faceit().gridIndex(Side::Lo);
                  Vector<FaceIndex> allFaces = a_newGraph.getAllFaces(ivhi, idir, Side::Lo);
                  if (allFaces.size() > 1)
                    {
                      //now we have the wacky case where two full, single-valued faces were coarsened
                      //to make a multi-valued face on the coarse level.
                      //also know as the bjorn geometry.  We need to infer the centroids from the
                      //volume centroids of the vofs involved.   The value of the centroid will be
                      //0.5 so i just need to get the sign.
                      //doing the best i can here.   there might be cases where this gets the
                      //centroid wrong but i cannot think of them
                      fullFace.m_areaFrac = 1.0/Real(allFaces.size());
                      VolIndex otherVoF;
                      if (a_vofsToChange.contains(ivlo) && (!a_vofsToChange.contains(ivhi)))
                        {
                          otherVoF = faceit().getVoF(Side::Hi);
                        }
                      else if (a_vofsToChange.contains(ivhi) && (!a_vofsToChange.contains(ivlo)))
                        {
                          otherVoF = faceit().getVoF(Side::Lo);
                        }
                      else
                        {
                          //you really should only be changing one of the vofs if there is a multivalued
                          //face between them
                          MayDay::Error("vofsToChange contains internal mutlivalued face");
                        }
                      const VolData& otherVolData = a_grownData(otherVoF, 0);
                      for (int tandir = 0; tandir < SpaceDim; tandir++)
                        {
                          fullFace.m_faceCentroid[tandir] = 1.0/Real(allFaces.size());
                          if (otherVolData.m_volCentroid[tandir] < 0)
                            {
                              fullFace.m_faceCentroid[tandir] *= -1;
                            }
                        }
                    }
                }
              m_faceData[idir](faceit(), 0) = fullFace;
            }
        }
    }
}
/************************/
/************************/
void EBDataImplem::
addEmptyIrregularVoFs(const IntVectSet& a_vofsToChange,
                      const EBGraph&    a_newGraph)
{
  if (!a_vofsToChange.isEmpty())
    {
      //calculate set by adding in new intvects
      const IntVectSet& ivsOld = m_volData.getIVS();
      IntVectSet ivsNew = ivsOld | a_vofsToChange;

      //for copying purposes
      Box minBox = ivsOld.minBox();
      Interval interv(0,0);

      //save old data into temporarys
      BaseIVFAB<VolData>  oldVolData(ivsOld, a_newGraph, 1);
      oldVolData.copy(minBox, interv, minBox, m_volData, interv);
      BaseIFFAB<FaceData> oldFaceData[SpaceDim];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          oldFaceData[idir].define(ivsOld, a_newGraph, idir, 1);
          oldFaceData[idir].copy(minBox, interv, minBox, m_faceData[idir], interv);
        }

      //redefine member data with new IVS and copy old data back
      //into it.
      m_volData.define(ivsNew, a_newGraph, 1);
      m_volData.copy(minBox, interv, minBox, oldVolData, interv);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_faceData[idir].define(ivsNew, a_newGraph, idir, 1);
          m_faceData[idir].copy(minBox, interv, minBox, oldFaceData[idir], interv);
        }

      //now put correct data into new vofs.  the volume fraction of a formally
      //regular cell is always unity
      VolData  emptyVol;
      emptyVol.setToCovered();
      for (VoFIterator vofit(a_vofsToChange, a_newGraph); vofit.ok(); ++vofit)
        {
          m_volData(vofit(), 0) = emptyVol;
        }
    }
}
/************************/
EBDataImplem::
EBDataImplem()
{
  m_isVoFDataDefined = false;
  m_isFaceDataDefined = false;
  m_hasMoments = false;
}
/************************/
EBDataImplem::
~EBDataImplem()
{
}
/************************/
void EBDataImplem::
copy(const Box&      a_regionFrom,
     const Interval& a_Cd,
     const Box&      a_regionTo,
     const EBDataImplem&   a_source,
     const Interval& a_Cs)
{
  CH_assert(m_isVoFDataDefined);
  CH_assert(m_isFaceDataDefined);
  CH_assert(a_source.m_isVoFDataDefined);
  CH_assert(a_source.m_isFaceDataDefined);
  Interval ivsca(0,0);
  m_volData.copy(a_regionFrom, ivsca, a_regionTo, a_source.m_volData ,ivsca);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_faceData[idir].copy(a_regionFrom, ivsca, a_regionTo, a_source.m_faceData[idir], ivsca);
    }

}
/************************/
void 
EBDataImplem::
define (const EBGraph&           a_graph,
        const Box&               a_region,
        const Real& a_dx,
        bool a_hasMoments)
{
  m_dx = a_dx;
  m_hasMoments = a_hasMoments;
  defineVoFData(a_graph, a_region);
  defineFaceData(a_graph, a_region);
  m_regularAreaMoments.setRegular(m_dx);
  m_regularVolumeMoments.setRegular(m_dx);
}
/************************/
void EBDataImplem::
defineVoFData(const EBGraph& a_graph,
              const Box& a_validBox)
{
  CH_TIME("EBDataImpem::defineVoFData");
  m_isVoFDataDefined = true;

  IntVectSet ivsIrreg = a_graph.getIrregCells(a_validBox);
  m_volData.define(ivsIrreg, a_graph, 1);
}
/************************/
void EBDataImplem::
defineFaceData(const EBGraph& a_graph,
               const Box& a_validBox)
{
  CH_TIME("EBDataImpem::defineFaceData");
  m_isFaceDataDefined = true;
  IntVectSet ivsIrreg = a_graph.getIrregCells(a_validBox);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_faceData[idir].define(ivsIrreg, a_graph, idir, 1);
    }
}
/************************/
void
EBDataImplem::define(const EBGraph&           a_graph,
                     const Vector<IrregNode>& a_irregGraph,
                     const Box&               a_validBox,
                     const Real &             a_dx,
                     bool                     a_hasMoments)

{
  CH_TIME("EBDataImpem::define");
  m_dx = a_dx;
  m_hasMoments = a_hasMoments;
  m_regularAreaMoments.setRegular(m_dx);
  m_regularVolumeMoments.setRegular(m_dx);

  defineVoFData( a_graph, a_validBox);
  defineFaceData(a_graph, a_validBox);

  if (a_graph.hasIrregular())
  {
    for (int inode = 0; inode < a_irregGraph.size(); inode++)
    {
      const IrregNode& node = a_irregGraph[inode];
      const IntVect& iv = node.m_cell;
      if (a_validBox.contains(iv))
      {
        const int&  cellInd = node.m_cellIndex;
        VolIndex vof(iv, cellInd);

        VolData& vol = m_volData(vof, 0);
        vol.m_volFrac       = node.m_volFrac;
        vol.m_volCentroid   = node.m_volCentroid;
        vol.m_averageFace.m_volIndex      = vof;
        vol.m_averageFace.m_bndryCentroid = node.m_bndryCentroid;
        vol.m_averageFace.m_bndryPhase    = -1;  // EB

        for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          for (SideIterator sit; sit.ok(); ++sit)
          {
            Vector<FaceIndex> faces = a_graph.getFaces(vof, faceDir, sit());
            int nodeind = node.index(faceDir, sit());
            Vector<Real> areaFracs         = node.m_areaFrac[nodeind];
            Vector<RealVect> faceCentroids = node.m_faceCentroid[nodeind];
            for (int iface = 0; iface < faces.size(); iface++)
            {
              const Real&     areaFracNode     = areaFracs[iface];
              const RealVect& faceCentroidNode = faceCentroids[iface];
              const FaceIndex& face = faces[iface];

              m_faceData[faceDir](face,0).m_areaFrac     = areaFracNode;
              m_faceData[faceDir](face,0).m_faceCentroid = faceCentroidNode;
            }
          }
        }
        //put in higher order moment stuff if it is there
        if(a_hasMoments)
        {
          VolData& vol = m_volData(vof, 0);
          vol.m_volumeMoments = node.m_volumeMoments;
          vol.m_averageFace.m_EBMoments = node.m_EBMoments;
          for (int idir = 0; idir < SpaceDim; idir++)
          {
            vol.m_averageFace.m_normalPartialDeriv[idir] = node.m_normalPartialDeriv[idir];
          }
          //use the eb moments and the partial derivatives of the normal 
          //to calculate the eb normal moments
          for (int normDir = 0; normDir < SpaceDim; normDir++)
          {
            calculateEBNormalMoments(vol.m_averageFace.m_EBNormalMoments[normDir], 
                                     vol.m_averageFace.m_normalPartialDeriv[normDir],
                                     vol.m_averageFace.m_EBMoments);
          }
          for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
          {
            for (SideIterator sit; sit.ok(); ++sit)
            {
              Vector<FaceIndex> faces = a_graph.getFaces(vof, faceDir, sit());
              int iindex = node.index(faceDir, sit());
              for (int iface = 0; iface < faces.size(); iface++)
              {
                const FaceIndex& face = faces[iface];
                m_faceData[faceDir](face,0).m_faceMoments = node.m_faceMoments[iindex];
              }//end list of faces on side
            }//end side iterator
          }//end face direction iterator

        }// end if has moments
      } //end if inside box
    } //end loop over nodes
  } //end if(hasIrregular)
}
/*******************************/
const Real& EBDataImplem::volFrac(const VolIndex& a_vof) const
{
  return m_volData(a_vof, 0).m_volFrac;
}

/*******************************/
const Real& EBDataImplem::bndryArea(const VolIndex& a_vof, int face) const
{
  const VolData& v =  m_volData(a_vof, 0);
  if (v.m_phaseFaces.size()>0)
    return v.m_phaseFaces[face].m_bndryArea;
  CH_assert(face == 0);
  return v.m_averageFace.m_bndryArea;
}

const Real& EBDataImplem::bndryArea(const VolIndex& a_vof) const
{
  static Real zero = 0;
  if (m_volData.getIVS().contains(a_vof.gridIndex()))
    return m_volData(a_vof, 0).m_averageFace.m_bndryArea;

  return zero;
}

/*******************************/
const RealVect& EBDataImplem::normal(const VolIndex& a_vof) const
{
  return m_volData(a_vof, 0).m_averageFace.m_normal;

}
const RealVect& EBDataImplem::normal(const VolIndex& a_vof, int face) const
{
  const VolData& v =  m_volData(a_vof, 0);
  if (v.m_phaseFaces.size()>0)
    return v.m_phaseFaces[face].m_normal;
  CH_assert(face == 0);
  return v.m_averageFace.m_normal;
}

/*******************************/
const RealVect& EBDataImplem::centroid(const VolIndex& a_vof) const
{
  return m_volData(a_vof, 0).m_volCentroid;
}
/*******************************/
const RealVect& EBDataImplem::bndryCentroid(const VolIndex& a_vof) const
{
  return m_volData(a_vof, 0).m_averageFace.m_bndryCentroid;
}
const RealVect& EBDataImplem::bndryCentroid(const VolIndex& a_vof, int face) const
{
  const VolData& v =  m_volData(a_vof, 0);
  if (v.m_phaseFaces.size()>0)
    return v.m_phaseFaces[face].m_bndryCentroid;
  CH_assert(face == 0);
  return v.m_averageFace.m_bndryCentroid;

}

  /// used by multi-fluid code
int EBDataImplem::facePhase(const VolIndex& a_vof, int face) const
{
  return m_volData(a_vof, 0).m_phaseFaces[face].m_bndryPhase;
}

  /// used by multi-fluid code
const VolIndex& EBDataImplem::faceIndex(const VolIndex& a_vof, int face) const
{
  return m_volData(a_vof, 0).m_phaseFaces[face].m_volIndex;
}

  /// used by multi-fluid code
void EBDataImplem::setFacePhase(const VolIndex& a_vof, int face, int phase)
{
  VolData& voldat = m_volData(a_vof, 0);
  voldat.m_phaseFaces[face].m_bndryPhase=phase;
}

  /// used by multi-fluid code
void EBDataImplem::setFaceIndex(const VolIndex& a_vof, int face, const VolIndex& index)
{
  m_volData(a_vof, 0).m_phaseFaces[face].m_volIndex = index;
}

///
int EBDataImplem::numFacePhase(const VolIndex& a_vof) const
{
  return  m_volData(a_vof, 0).m_phaseFaces.size();
}

/*******************************/
const RealVect& EBDataImplem::centroid(const FaceIndex& a_face) const
{
  int faceDir = a_face.direction();
  return m_faceData[faceDir](a_face, 0).m_faceCentroid;
}
/*******************************/
const Real& EBDataImplem::areaFrac(const FaceIndex& a_face) const
{
  int faceDir = a_face.direction();
  return m_faceData[faceDir](a_face, 0).m_areaFrac;
}
/*******************************/
void
EBDataImplem::
coarsenBoundaryAreaAndNormal(Real&                    a_bndryAreaCoar,
                             RealVect&                a_normalCoar,
                             const Vector<Real>&      a_bndryAreaFine,
                             const Vector<RealVect>&  a_normalFine)
{
  CH_TIME("EBDataImplem::coarsenBoundaryAreaAndNormal");

   Real faceCoarsenFactor = D_TERM(1.0, * 0.5, *0.5);
   // Real voluCoarsenFactor = D_TERM(0.5, * 0.5, *0.5);

  //A^B_xi, C = sum(A^B_xi, F/2^D-1)
  RealVect bndryAreaVec= RealVect::Zero;
  for (int ifine = 0; ifine < a_normalFine.size(); ifine++)
    {
      RealVect normalFine = a_normalFine[ifine];
      Real bndryAreaFine =  a_bndryAreaFine[ifine];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          bndryAreaVec[idir] += normalFine[idir]*bndryAreaFine;
        }
    }
  bndryAreaVec *= faceCoarsenFactor;
  //A^B,C = ||A^B_xi,C||
  a_bndryAreaCoar = 0.;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_bndryAreaCoar += bndryAreaVec[idir]*bndryAreaVec[idir];
    }
  a_bndryAreaCoar = sqrt(a_bndryAreaCoar);
//  for (int ifine = 0; ifine < a_normalFine.size(); ifine++)
//    {
//      a_bndryAreaCoar += a_bndryAreaFine[ifine];
//    }
//  a_bndryAreaCoar *= faceCoarsenFactor;


  //n_xi, C = A^B_xi,C/AB,c
  if (a_bndryAreaCoar > 0.)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          a_normalCoar[idir] = bndryAreaVec[idir]/a_bndryAreaCoar;
        }
    }
  else
    {
      a_normalCoar = RealVect::Zero;
    }
}


void   EBDataImplem::setBoundaryPhase(int phase)
{
  CH_assert(m_volData.nComp() == 1);
  VolData*  p = m_volData.dataPtr(0);
  for (int i=0; i<m_volData.numVoFs(); ++i)
    {
      p[i].m_averageFace.m_bndryPhase = phase;
    }

}

void   EBDataImplem::clearMultiBoundaries()
{
  CH_assert(m_volData.nComp() == 1);
  if (m_volData.numVoFs() > 0)
    {
      VolData*  p = m_volData.dataPtr(0);
      for (int i=0; i<m_volData.numVoFs(); ++i)
        {
          p[i].m_phaseFaces.clear();
        }
    }

}
/*******************************/
void 
EBDataImplem::
setAreaMomentsToZero(const FaceIndex& a_face)
{
  m_faceData[a_face.direction()](a_face, 0).m_faceMoments.setToZero();
}
/*******************************/
///shift input by shift and increment output
void
EBDataImplem::
shiftAndIncrement(IndMomSpaceDim& a_output, const IndMomSpaceDim& a_input,const RealVect& a_shiftRV)
{
  IndexTM<Real, SpaceDim> shiftVec;
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    shiftVec[idir] =  a_shiftRV[idir];
  }
  IndMomSpaceDim increment = a_input;
  increment.shift(shiftVec);
  
  a_output += increment;
}

/*******************************/
///shift input by shift and increment output
void
EBDataImplem::
shiftAndIncrement(IndMomSDMinOne& a_output, const IndMomSDMinOne& a_input,const RealVect& a_shiftRV, int faceDir)
{
  IndexTM<Real, SpaceDim-1> shiftVec;
  int iindex = 0;
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    if(idir != faceDir)
    {
      shiftVec[iindex] =  a_shiftRV[idir];
      iindex++;
    }
  }
  IndMomSDMinOne increment = a_input;
  increment.shift(shiftVec);
  
  a_output += increment;
}
/*******************************/
void 
EBDataImplem::
setVolumeMomentsToZero(const VolIndex& a_vof)
{
  m_volData(a_vof, 0).m_volumeMoments.setToZero();
  m_volData(a_vof, 0).m_averageFace.m_EBMoments.setToZero();
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    m_volData(a_vof, 0).m_averageFace.m_normalPartialDeriv[idir].setToZero();
    m_volData(a_vof, 0).m_averageFace.m_EBNormalMoments[idir].setToZero();
  }
}
/*******************************/

void EBDataImplem::
coarsenVoFs(const EBDataImplem&  a_fineEBDataImplem,
            const EBGraph&       a_fineGraph,
            const EBGraph&       a_coarGraph,
            const Box&           a_validRegion)
{
  CH_TIME("EBDataImplem::coarsenVoFs");

  defineVoFData(a_coarGraph, a_validRegion);

  if (a_coarGraph.hasIrregular())
    {
      IntVectSet ivsIrreg = a_coarGraph.getIrregCells(a_validRegion);
      std::list<BoundaryData>  boundary;

      for (VoFIterator vofit(ivsIrreg, a_coarGraph); vofit.ok(); ++vofit)
        {
          CH_TIME("EBDataImplem::coarsenVoFs_VoFIterator");
          const VolIndex& vofCoar = vofit();
          Vector<VolIndex> vofsFine = a_coarGraph.refine(vofCoar);
          int nFine = vofsFine.size();
          Vector<Real> bndryAreaFine(nFine);
          Vector<Real> volFracFine(nFine);
          Vector<int>  phase(nFine);
          Vector<RealVect> bndryCentroidFine(nFine);
          Vector<RealVect> volCentroidFine(nFine);
          Vector<RealVect> normalFine(nFine);

          for (int ifine = 0; ifine < nFine; ifine++)
            {
              CH_TIME("EBDataImplem::coarsenVoFs_fine");
              const VolIndex& vofFine =vofsFine[ifine];

              if (a_fineGraph.isIrregular(vofFine.gridIndex()))
                {
                  const VolData& vol = a_fineEBDataImplem.m_volData(vofFine,0);
                  bndryAreaFine[ifine] = vol.m_averageFace.m_bndryArea;
                  //a_fineEBDataImplem.bndryArea(vofFine);

                  volFracFine[ifine] = vol.m_volFrac;
                  //  a_fineEBDataImplem.volFrac(vofFine);
                  bndryCentroidFine[ifine] = vol.m_averageFace.m_bndryCentroid;
                  // a_fineEBDataImplem.bndryCentroid(vofFine);
                  volCentroidFine[ifine] = vol.m_volCentroid;
                  // a_fineEBDataImplem.centroid(vofFine);
                  normalFine[ifine] = vol.m_averageFace.m_normal;
                  //  a_fineEBDataImplem.normal(vofFine);
                  if (vol.m_phaseFaces.size()>0)
                  {
                    for (int i=0; i<vol.m_phaseFaces.size(); i++)
                    {
                      boundary.push_back(vol.m_phaseFaces[i]);
                    }
                  }
                  else
                  {
                    boundary.push_back(vol.m_averageFace);
                  }

                }
              else
                {
                  CH_assert(a_fineGraph.isRegular(vofFine.gridIndex()));
                  bndryAreaFine[ifine] = 0.0;
                  volFracFine[ifine] = 1.0;
                  bndryCentroidFine[ifine] = RealVect::Zero;
                  volCentroidFine[ifine] = RealVect::Zero;
                  normalFine[ifine] = RealVect::Zero;
                }
            }

          Real volFracCoar, bndryAreaCoar;
          RealVect volCentroidCoar, bndryCentroidCoar, normalCoar;

          coarsenVolFracAndCentroid(volFracCoar,
                                    volCentroidCoar,
                                    volFracFine,
                                    volCentroidFine,
                                    vofsFine,
                                    vofCoar);

          coarsenBoundaryAreaAndNormal(bndryAreaCoar,
                                       normalCoar,
                                       bndryAreaFine,
                                       normalFine);

          coarsenBndryCentroid(bndryCentroidCoar,
                               bndryCentroidFine,
                               bndryAreaFine,
                               vofsFine,
                               vofCoar);

          VolData& vol = m_volData(vofCoar,0);

          vol.m_volFrac       = volFracCoar;
          vol.m_averageFace.m_bndryArea     = bndryAreaCoar;
          vol.m_averageFace.m_volIndex      = vofCoar;
          vol.m_averageFace.m_bndryCentroid = bndryCentroidCoar;
          vol.m_averageFace.m_normal        = normalCoar;
          vol.m_volCentroid   = volCentroidCoar;

          vol.m_phaseFaces.resize(boundary.size());
          for (int i=0; i<vol.m_phaseFaces.size(); ++i)
          {
            vol.m_phaseFaces[i] = boundary.front();
            boundary.pop_front();
            //vol.m_bndryPhase = phase[0];  // two-phase flow assumption
          }

          if(m_hasMoments)
          {
            setVolumeMomentsToZero(vofCoar);
            //grab the fine volume moments, shift them and add them up.
            IndMomSpaceDim& coarVoMom = m_volData(vofCoar, 0).m_volumeMoments;
            IndMomSpaceDim& coarEBMom = m_volData(vofCoar, 0).m_averageFace.m_EBMoments;

            for(int ivoffine = 0; ivoffine < vofsFine.size(); ivoffine++)
            {
              Real dxFine = a_fineEBDataImplem.m_dx;
              const VolIndex& vofFine = vofsFine[ivoffine];
              RealVect fineLoc = EBArith::getVoFLocation(vofFine, dxFine, RealVect::Zero);
              RealVect coarLoc = EBArith::getVoFLocation(vofCoar,   m_dx, RealVect::Zero);
              //              RealVect shiftAmt= coarLoc - fineLoc;
              RealVect shiftAmt= fineLoc - coarLoc;
              IndMomSpaceDim fineVoMom = a_fineEBDataImplem.getVolumeMoments(vofFine);
              IndMomSpaceDim fineEBMom = a_fineEBDataImplem.getEBMoments(vofFine);
              shiftAndIncrement(coarVoMom, fineVoMom, shiftAmt);
              shiftAndIncrement(coarEBMom, fineEBMom, shiftAmt); 
              for(int idir = 0; idir < SpaceDim; idir++)
              {
                IndMomSpaceDim fineEBNormMom = a_fineEBDataImplem.getEBNormalMoments(vofFine, idir);
                shiftAndIncrement(m_volData(vofCoar, 0).m_averageFace.m_EBNormalMoments[idir], fineEBNormMom, shiftAmt);
                //no idea what else to do here
                IndMomSpaceDim fineEBDeriv = a_fineEBDataImplem.getEBNormalPartialDerivs(vofFine, idir);
                m_volData(vofCoar, 0).m_averageFace.m_normalPartialDeriv[idir] = fineEBDeriv;
              }

            }
          }

        }
    }
}
/*******************************/
void EBDataImplem::
coarsenFaces(const EBDataImplem& a_fineEBDataImplem,
             const EBGraph&      a_fineGraph,
             const EBGraph&      a_coarGraph,
             const Box&          a_validRegion)
{
  CH_TIME("EBDataImplem::coarsenFaces");

  defineFaceData(a_coarGraph, a_validRegion);
  IntVectSet ivsIrreg = a_coarGraph.getIrregCells(a_validRegion);
  Box fineRegion = a_fineGraph.getRegion();
  if (a_coarGraph.hasIrregular())
    {
      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          CH_TIME("EBDataImplem::coarsenFaces_faceDir");

          FaceIterator faceit(ivsIrreg, a_coarGraph, faceDir,
                              FaceStop::SurroundingWithBoundary);

          for (faceit.reset(); faceit.ok(); ++faceit)
            {
              CH_TIME("EBDataImplem::coarsenFaces_FaceIterator");

              const FaceIndex&  faceCoar  = faceit();
              Vector<FaceIndex> facesFine = a_coarGraph.refine(faceCoar, a_fineGraph);

              Vector<Real>     areaFracsFine(facesFine.size());
              Vector<RealVect> centroidsFine(facesFine.size());
              for (int ifine = 0; ifine < facesFine.size(); ifine++)
                {
                  CH_TIME("EBDataImplem::coarsenFaces_fine");

                  const FaceIndex& faceFine = facesFine[ifine];
                  IntVect loiv = faceFine.gridIndex(Side::Lo);
                  IntVect hiiv = faceFine.gridIndex(Side::Hi);
                  if ((fineRegion.contains(loiv) && a_fineGraph.isIrregular(loiv)) ||
                     (fineRegion.contains(hiiv) && a_fineGraph.isIrregular(hiiv)))
                    {
                      areaFracsFine[ifine] = a_fineEBDataImplem.areaFrac(faceFine);
                      centroidsFine[ifine] = a_fineEBDataImplem.centroid(faceFine);
                    }
                  else
                    {
                      areaFracsFine[ifine] = 1.0;
                      centroidsFine[ifine] = RealVect::Zero;
                    }
                }
              Real areaFracCoar;
              RealVect centroidCoar;
              coarsenAreaFrac(areaFracCoar, areaFracsFine);

              coarsenFaceCentroid(centroidCoar, centroidsFine,
                                  areaFracsFine, facesFine,
                                  faceCoar);

              m_faceData[faceDir](faceCoar, 0).m_areaFrac     = areaFracCoar;
              m_faceData[faceDir](faceCoar, 0).m_faceCentroid = centroidCoar;
              if(m_hasMoments)
              {
                BaseIFFAB<FaceData>& coarFaceData = m_faceData[faceDir];
                setAreaMomentsToZero(faceCoar);
                RealVect coarLoc = EBArith::getFaceLocation(faceCoar, m_dx, RealVect::Zero);

                IndMomSDMinOne& coarFaceMom = coarFaceData(faceCoar, 0).m_faceMoments;
                for(int iface = 0; iface < facesFine.size(); iface++)
                {
                  const FaceIndex& faceFine = facesFine[iface];
                  Real dxFine = a_fineEBDataImplem.m_dx;
                  RealVect fineLoc = EBArith::getFaceLocation(faceFine, dxFine, RealVect::Zero);
                  RealVect shiftAmt= fineLoc - coarLoc;
                  shiftAmt[faceCoar.direction()] = 0;
                  IndMomSDMinOne fineFaceMom = a_fineEBDataImplem.getFaceMoments(faceFine);
                  shiftAndIncrement(coarFaceMom, fineFaceMom, shiftAmt, faceDir);
                }
              }

            } //end loop over faces
        } //end loop over face directions
    }
}
/*******************************/
void EBDataImplem::
coarsenFaceCentroid(RealVect&                a_centroidCoar,
                    const Vector<RealVect>&  a_centroidsFine,
                    const Vector<Real>&      a_areaFracFine,
                    const Vector<FaceIndex>& a_facesFine,
                    const FaceIndex&         a_faceCoar)
{
  CH_TIME("EBDataImplem::coarsenFaceCentroid");

  Real totalFaceArea = 0.0;
  for (int ifineFace = 0; ifineFace < a_facesFine.size(); ifineFace++)
    {
      totalFaceArea += a_areaFracFine[ifineFace];
    }

  a_centroidCoar = RealVect::Zero;
  if (totalFaceArea > 0.0)
    {
      const IntVect& coarIV = a_faceCoar.gridIndex(Side::Lo);
      for (int ifineFace = 0; ifineFace < a_facesFine.size(); ifineFace++)
        {
          const FaceIndex& fineFace        = a_facesFine[ifineFace];
          IntVect fineIV = fineFace.gridIndex(Side::Lo);

          RealVect centroidFine = a_centroidsFine[ifineFace];
          centroidFine = fineToCoarseTransform(centroidFine,
                                               coarIV, fineIV);

          const Real&   areaFracFine= a_areaFracFine[ifineFace];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              a_centroidCoar[idir] +=
                centroidFine[idir]*(areaFracFine/totalFaceArea);
            }
        }
    }
}
/*******************************/
void EBDataImplem::
coarsenAreaFrac(Real& a_areaFracCoar,
                const Vector<Real>& a_areaFracFine)
{
  CH_TIME("EBDataImplem::coarsenAreaFrac");
  //this is the factor by which the area of a fine
  //face is smaller than the area of a coarse face.
  Real faceCoarsenFactor = 1.0;
  for (int dir = 0; dir < SpaceDim-1; ++dir)
    {
      faceCoarsenFactor *= 0.5;
    }
  a_areaFracCoar = 0.0;
  for (int ifine = 0; ifine < a_areaFracFine.size(); ifine++)
    {
      a_areaFracCoar += a_areaFracFine[ifine];
    }
  a_areaFracCoar *= faceCoarsenFactor;
}
/*******************************/
void
EBDataImplem::
coarsenVolFracAndCentroid(Real&                   a_volFracCoar,
                          RealVect&               a_volCentroidCoar,
                          const Vector<Real>&     a_volFracFine,
                          const Vector<RealVect>& a_volCentroidFine,
                          const Vector<VolIndex>& a_fineVoFs,
                          const VolIndex&         a_coarVoF)
{
  CH_TIME("EBDataImplem::coarsenVolFracAndCentroid");

  Real volCoarsenFactor = 1.0;
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      volCoarsenFactor *= 0.5;
    }

  Real totalVol = 0;
  for (int ifine = 0; ifine < a_fineVoFs.size(); ifine++)
    {
      totalVol  += a_volFracFine[ifine];
    }
  a_volFracCoar = volCoarsenFactor*totalVol;

  a_volCentroidCoar = RealVect::Zero;
  for (int ifine = 0; ifine < a_fineVoFs.size(); ifine++)
    {
      const VolIndex& fineVoF= a_fineVoFs[ifine];
      const IntVect& coarIV = a_coarVoF.gridIndex();
      const IntVect& fineIV = fineVoF.gridIndex();
      RealVect volCentroidFine = a_volCentroidFine[ifine];
      Real volFracFine=  a_volFracFine[ifine];
      volCentroidFine = fineToCoarseTransform(volCentroidFine,
                                              coarIV, fineIV);

      if (totalVol > 0.0)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              a_volCentroidCoar[idir] +=
                volCentroidFine[idir]*(volFracFine/totalVol);
            }
        }
    }
}
/*******************************/
RealVect
EBDataImplem::
fineToCoarseTransform(const RealVect& a_finePoint,
                      const IntVect&  a_coarCell,
                      const IntVect&  a_fineCell)
{
  RealVect retval;
  //assuming nref = 2. make dxf = 1
  Real dxc = 2.0;
  RealVect fineCellLoc, coarCellLoc;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      fineCellLoc[idir] = Real(a_fineCell[idir]) + 0.5;
      coarCellLoc[idir] = dxc*(Real(a_coarCell[idir]) + 0.5);
    }
  retval = a_finePoint+ fineCellLoc;
  retval -= coarCellLoc;
  //to put it into a space where dxc = 1, the normalized answer
  retval /= dxc;
  return retval;
}
/*******************************/
void
EBDataImplem::
coarsenBndryCentroid(RealVect&               a_bndryCentroidCoar,
                     const Vector<RealVect>& a_bndryCentroidFine,
                     const Vector<Real>&     a_bndryAreaFine,
                     const Vector<VolIndex>& a_fineVoFs,
                     const VolIndex&         a_coarVoF)
{
  CH_TIME("EBDataImplem::coarsenBndryCentroid");

  Real totalArea = 0;
  for (int ifine = 0; ifine < a_fineVoFs.size(); ifine++)
    {
      totalArea += a_bndryAreaFine[ifine];
    }

  a_bndryCentroidCoar = RealVect::Zero;
  for (int ifine = 0; ifine < a_fineVoFs.size(); ifine++)
    {
      const VolIndex& fineVoF= a_fineVoFs[ifine];
      const IntVect& fineIV=  fineVoF.gridIndex();
      const IntVect& coarIV=a_coarVoF.gridIndex();

      Real bndryAreaFine = a_bndryAreaFine[ifine];

      RealVect bndryCentroidFine = a_bndryCentroidFine[ifine];

      bndryCentroidFine =
        fineToCoarseTransform(bndryCentroidFine,coarIV, fineIV);
      if (totalArea > 0.0)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              a_bndryCentroidCoar[idir] +=
                bndryCentroidFine[idir]*(bndryAreaFine/totalArea);
            }
        }
    }
}
/*******************************/
int
EBDataImplem::size(const Box& a_region, const Interval& a_comps) const
{
  CH_TIME("ebdataimplem::size");
  CH_assert(m_isFaceDataDefined);
  CH_assert(m_isVoFDataDefined);
  int linearSize = m_volData.size(a_region, a_comps);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      linearSize += m_faceData[idir].size(a_region, a_comps);
    }
  if (s_verbose)
    {
      pout() << " EBDataImplem::size returning " << linearSize << " for box " << a_region << endl;
    }
  return linearSize;
}
/*******************************/
void
EBDataImplem::linearOut(void* a_buf,
                        const Box& a_region,
                        const Interval& a_comps) const
{
  CH_TIME("ebdataimplem::linearout");
  CH_assert(m_isFaceDataDefined);
  CH_assert(m_isVoFDataDefined);

  if (s_verboseDebug)
    {
      BaseIVFAB<VolData>::setVerboseDebug(true);
    }
//   int linearSize = 0;
  unsigned char* buffer = (unsigned char* ) a_buf;
  m_volData.linearOut(buffer, a_region, a_comps);
  buffer     +=   m_volData.size(a_region, a_comps);
  if (s_verboseDebug)
    {
      pout() << " EBDataImplem::linearOut volData size " << m_volData.size(a_region, a_comps) << " for box " << a_region << endl;
    }
  //linearSize += m_volData.size(a_region, a_comps);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_faceData[idir].linearOut(buffer, a_region, a_comps);
      buffer +=       m_faceData[idir].size(a_region, a_comps);
      if (s_verboseDebug)
        {
          pout() << " EBDataImplem::linearOut faceData size " << m_faceData[idir].size(a_region, a_comps) << " for box " << a_region << endl;
        }
      //linearSize += m_faceData[idir].size(a_region, a_comps);
    }
  if (s_verboseDebug)
    {
      BaseIVFAB<VolData>::setVerboseDebug(false);
    }
}
/*******************************/
void
EBDataImplem::linearIn(void*           a_buf,
                       const Box&      a_region,
                       const Interval& a_comps)
{
  CH_TIME("ebdataimplem::linearin");
  CH_assert(m_isFaceDataDefined);
  CH_assert(m_isVoFDataDefined);

  if (s_verboseDebug)
    {
      BaseIVFAB<VolData>::setVerboseDebug(true);
    }
//   int linearSize = 0;
  unsigned char* buffer = (unsigned char* ) a_buf;
  m_volData.linearIn(buffer, a_region, a_comps);

  buffer     +=   m_volData.size(a_region, a_comps);
  if (s_verboseDebug)
    {
      pout() << " EBDataImplem::linearIn volData size " << m_volData.size(a_region, a_comps) << " for box " << a_region << endl;
    }
//   linearSize += m_volData.size(a_region, a_comps);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
//       if (s_verboseDebug)
//         {
//           pout() << " EBDataImplem::linearIn faceData size before " << m_faceData[idir].size(a_region, a_comps) << " for box " << a_region << endl;
//         }
      m_faceData[idir].linearIn(buffer, a_region, a_comps);
      buffer +=       m_faceData[idir].size(a_region, a_comps);
      if (s_verboseDebug)
        {
          pout() << " EBDataImplem::linearIn faceData size after " << m_faceData[idir].size(a_region, a_comps) << " for box " << a_region << endl;
        }
//       linearSize += m_faceData[idir].size(a_region, a_comps);
    }

  if (s_verboseDebug)
    {
      BaseIVFAB<VolData>::setVerboseDebug(false);
    }
}
/*******************************/
/*******************************/
/*******************************/
void EBData::
computeNormalsAndBoundaryAreas(const EBGraph& a_graph,
                               const Box&     a_validRegion)
{
  EBISBox ebisBox;
  DataIndex dummydit;
  ebisBox.define(a_graph, *this, dummydit);
  BaseIVFAB<VolData>& volData = m_implem->getVolData();
  if (a_graph.hasIrregular())
    {
      IntVectSet ivsIrreg = a_graph.getIrregCells(a_validRegion);
      for (VoFIterator vofit(ivsIrreg, a_graph); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Real bndryArea  =  PolyGeom::bndryArea(vof, ebisBox);
          RealVect normal =  PolyGeom::normal(   vof, ebisBox, bndryArea);

          //copy results to volData
          volData(vof,0).m_averageFace.m_bndryArea = bndryArea;
          volData(vof,0).m_averageFace.m_normal    = normal;
        }
    }

}

#include "NamespaceFooter.H"
