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

#ifndef _EBDATA_H_
#define _EBDATA_H_

#include "REAL.H"
#include "RealVect.H"

#include "EBGraph.H"
#include "IrregNode.H"
#include "BaseIVFAB.H"
#include "BaseIFFAB.H"

#include "NamespaceHeader.H"

class EBIndexSpace;

template <class T>
int linearSize(const T& inputT);

template <class T>
void linearIn(T& a_outputT, const void* const inBuf);

template <class T>
void linearOut(void* const a_outBuf, const T& inputT);

class BoundaryData
{
public:
  BoundaryData();

  BoundaryData(const Real& bndryArea, const RealVect& normal,
               const RealVect& bndryCentroid, int phase,const VolIndex& index)
    :
    m_bndryArea(bndryArea),
    m_normal(normal),
    m_bndryCentroid(bndryCentroid),
    m_bndryPhase(phase),
    m_volIndex(index)
    {
    }

  Real     m_bndryArea;
  RealVect m_normal;
  RealVect m_bndryCentroid;
  int      m_bndryPhase;
  VolIndex m_volIndex;


  /// the moment at  the irregular face associated with the  monomial with the input exponents 
  /**
     Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
     returns integral_over_irregular_area((x^p y^q z^r) dA)
  **/
  IndMomSpaceDim m_EBMoments;


  /// the derivatives of the normal
  /**
     (d_p)(n_d)
  **/
  IndMomSpaceDim m_normalPartialDeriv[SpaceDim];

  /// ebmoments with the the normal integrated in
  /**
     Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
     returns integral_over_irregular_area((x^p y^q z^r)n_i dA)
  **/
  IndMomSpaceDim m_EBNormalMoments[SpaceDim];

  void
  setToZero()
    {
      m_bndryArea = 0.;
      m_bndryCentroid = RealVect::Zero;
      m_normal = BASISV(0); //have to give it something non-zero
      m_EBMoments.setToZero();
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        m_normalPartialDeriv[idir].setToZero();
        m_EBNormalMoments   [idir].setToZero();
        
      }
    }
  
  BoundaryData& operator=(const BoundaryData& a_in)
    {
      if(&a_in != this)
      {
        m_bndryArea         = a_in.m_bndryArea;
        m_bndryCentroid     = a_in.m_bndryCentroid;
        m_bndryPhase        = a_in.m_bndryPhase;
        m_normal            = a_in.m_normal;
        m_volIndex          = a_in.m_volIndex;
        m_EBMoments         = a_in.m_EBMoments;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          m_normalPartialDeriv[idir] = a_in.m_normalPartialDeriv[idir];
          m_EBNormalMoments   [idir] = a_in.m_EBNormalMoments   [idir];
        }
      }
      return *this;
    }
};

template < >
int linearSize(const BoundaryData& vdata);
template < >
void linearIn(BoundaryData& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const BoundaryData& a_inputT);
///
/**
 */
class VolData
{
public:
  VolData();
  Real                 m_volFrac;
  RealVect             m_volCentroid;
  BoundaryData         m_averageFace;
  IndMomSpaceDim       m_volumeMoments;
  Vector<BoundaryData> m_phaseFaces;
  void setToRegular(Real a_dx)
    {
      m_volFrac = 1.0;
      m_volCentroid = RealVect::Zero;
      m_volumeMoments.setRegular(a_dx);
      m_averageFace.setToZero();
      m_phaseFaces.resize(0);
    }

  void setToCovered()
    {
      m_volFrac = 0.0;
      m_volCentroid = RealVect::Zero;
      m_averageFace.setToZero();
      m_volumeMoments.setToZero();
      m_phaseFaces.resize(0);
    }
};

/// VolData specializations of linear* functions used by LevelData
//[NOTE: see also SPMD.H. <dbs>]
template < >
int linearSize(const VolData& vdata);
template < >
void linearIn(VolData& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const VolData& a_inputT);

///
/**
 */
class FaceData
{
public:
  Real     m_areaFrac;
  RealVect m_faceCentroid;
  IndMomSDMinOne  m_faceMoments;
  void setToRegular(Real a_dx)
    {
      m_areaFrac = 0.;
      m_faceCentroid = RealVect::Zero;
      m_faceMoments.setRegular(a_dx);
    }
  void setToCovered()
    {
      m_areaFrac = 0.;
      m_faceCentroid = RealVect::Zero;
      m_faceMoments.setToZero();
    }
};


template < >
int linearSize(const FaceData& vdata);
template < >
void linearIn(FaceData& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const FaceData& a_inputT);

///
/**
   This class contains all the geometric information
   for an ebisbox.
*/
class EBDataImplem
{

public:

  ///
  EBDataImplem();

  ///
  ~EBDataImplem();

  ///
  /**
     noop---real defines are more complicated
  */
  void define(const Box& box, int comps)
    {
    }


  ///
  /**
     noop --- real constructor more complicated
  */
  EBDataImplem(const Box& a_box, int a_comps)
    {
      define(a_box, a_comps);
    }

  ///
  /**
     Copy the information from a_source to the over the intersection
     of the box a_region, the box of the current EBDataImplem and
     the box of a_source.  The Interval arguments are ignored.
     This function is required by LevelData.
  */
  void copy(const Box&      a_regionFrom,
            const Interval& a_Cd,
            const Box&      a_regionto,
            const EBDataImplem&   a_source,
            const Interval& a_Cs);

  ///define the whole thing
  void
  define(const EBGraph&           a_graph,
         const Vector<IrregNode>& a_irregData,
         const Box&               a_validBox,
         const Real             & a_dx,
         bool                     a_hasMoments );

  //this is the define to use when you are just going to copy the data into it.  region should include ghost cells
  void define(const EBGraph&           a_graph,
              const Box&               a_region,
              const Real&              a_dx,
              bool a_hasMoments) ;

  ///
  void
  coarsenVoFs(const EBDataImplem&   a_fineEBDataImplem,
              const EBGraph&        a_fineGraph,
              const EBGraph&        a_coarGraph,
              const Box&           a_validRegion);

  ///
  void
  coarsenFaces(const EBDataImplem& a_fineEBDataImplem,
               const EBGraph&      a_fineGraph,
               const EBGraph&      a_coarGraph,
               const Box&          a_validRegion);
  ///
  const Real& volFrac(const VolIndex& a_vof) const;

  ///
  const Real& areaFrac(const FaceIndex& a_face1) const;

  ///
  const RealVect& centroid(const FaceIndex& facein) const;

  ///
  const RealVect& centroid(const VolIndex& a_vof) const;

  ///
  const RealVect& bndryCentroid(const VolIndex& a_vof) const;
  const RealVect& bndryCentroid(const VolIndex& a_vof, int face) const;

  ///
  const Real& bndryArea(const VolIndex& a_vof) const;
  const Real& bndryArea(const VolIndex& a_vof, int face) const;

  ///
  const RealVect& normal(const VolIndex& a_vof) const;
  const RealVect& normal(const VolIndex& a_vof, int face) const;

  //  New set of boundary inquiry functions to support multi-fluid applications
  /// used by multi-fluid applications

  /// used by multi-fluid code
  int facePhase(const VolIndex& a_vof, int aface) const ;

  /// used by multi-fluid code
  const VolIndex& faceIndex(const VolIndex& a_vof, int face) const ;

  /// used by multi-fluid code
  void setFacePhase(const VolIndex& a_vof, int face, int phase);
  ///
  /// used by multi-fluid code
  void setFaceIndex(const VolIndex& a_vof, int face, const VolIndex& index);

  ///
  int numFacePhase(const VolIndex& a_vof) const ;

  void clearMultiBoundaries(); 

  void setBoundaryPhase(int phase);
  /**
     This stuff required by LevelData in parallel:
  */
  int size(const Box& R, const Interval& comps) const;

  ///
  void linearOut(void* buf, const Box& R, const Interval& comps) const;

  ///
  void linearIn(void* buf, const Box& R, const Interval& comps);

  ///
  static int preAllocatable()
    {
      return 2; // dyanmic allocatable.
    }

  BaseIVFAB<VolData>& getVolData()
    {
      return m_volData;
    }
  const BaseIVFAB<VolData>& getVolData() const
    {
      return m_volData;
    }


  ///return true if higher order moments are available
  bool hasMoments() const
    {
      return m_hasMoments;
    }
  ///
  void addFullIrregularVoFs(const IntVectSet& a_vofsToChange,
                            const EBGraph&    a_newGhostGraph,
                            const BaseIVFAB<VolData>&     a_newGhostData,
                            const EBGraph&    a_oldGhostGraph);

  ///multifluid angels dancing on the heads of pins.
  void addEmptyIrregularVoFs(const IntVectSet& a_vofsToChange,
                             const EBGraph&    a_newGraph);
  ///
  static void setVerbose(bool a_verbose);
  ///
  static void setVerboseDebug(bool a_verboseDebug);


  /// get the moment at  the VoF associated with the  monomial with the input exponents
  /**
     Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
     returns integral_over_VoF(x^p y^q z^r dV) for all p q
  **/
  IndMomSpaceDim getVolumeMoments(const VolIndex& a_vof) const;

  /// get the normal the irregular face associated with the  monomial with the input exponents 
  /**
     Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
     returns integral_over_VoF((x^p y^q z^r) dV) for p q r
  **/
  IndMomSpaceDim getEBMoments(const VolIndex& a_vof) const;

  /// get the moment at  the face associated with the  monomial with the input exponents
  /**
     Given face variables x, y,  p = mono(0), q = mono(1)),
     returns integral_over_face(x^p y^q dA) for all p q
  **/
  IndMomSDMinOne getFaceMoments(const FaceIndex& a_face) const;

  /// I am sick of guessing this name wrong
  IndMomSDMinOne getAreaMoments(const FaceIndex& a_face) const
    {
      return getFaceMoments(a_face);
    }

  /// 
  /**
     gets the partial derivs of the normal component.   
  **/
  IndMomSpaceDim getEBNormalPartialDerivs(const VolIndex& a_vof, int normalComponent) const;

  /// get the normal*moment at  the irregular face associated with the  monomial with the input exponents 
  /**
     Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
     returns integral_over_VoF(((x-x0)^p (y-y0)^q (z-z0)^r)*n_i dV) for p q r,
     and (x0,y0,z0) is the center of the cell 
  **/
  IndMomSpaceDim getEBNormalMoments(const VolIndex& a_vof, int normalComponent) const;

private:


  bool 
  irregFace(const FaceIndex& a_face) const;

  bool 
  irregVoF(const VolIndex& a_vof) const;

  /// each data holder is defined over the irregular cells of the graph
  void
  defineVoFData(const EBGraph& a_graph, const Box& a_region);

  /// each data holder is defined over the irregular cells of the graph
  void
  defineFaceData(const EBGraph& a_graph, const Box& a_region);

  void setVolumeMomentsToZero(const VolIndex& a_vof);

  void setAreaMomentsToZero(const FaceIndex& a_face);

  void setCoveredAndRegular();

  void shiftAndIncrement(IndMomSpaceDim& a_output, const IndMomSpaceDim& a_input,const RealVect& a_shiftRV);
  void shiftAndIncrement(IndMomSDMinOne& a_output, const IndMomSDMinOne& a_input,const RealVect& a_shiftRV, int faceDir);

  static   bool s_verbose;
  static   bool s_verboseDebug;
  ///
  BaseIVFAB<VolData>         m_volData;

  ///
  IndMomSDMinOne m_regularAreaMoments;
 
  ///
  IndMomSpaceDim m_regularVolumeMoments;

  ///
  BaseIFFAB<FaceData>        m_faceData[SpaceDim];

  ///
  bool m_isFaceDataDefined;

  ///
  bool m_isVoFDataDefined;

  ///
  bool m_hasMoments;

  ///
  Real m_dx;
  void operator=(const EBDataImplem& ebiin)
    {;}

  EBDataImplem(const EBDataImplem& ebiin)
    {;}

  void
  coarsenFaceCentroid(RealVect&                a_centroidCoar,
                      const Vector<RealVect>&  a_centroidsFine,
                      const Vector<Real>&      a_areaFracFine,
                      const Vector<FaceIndex>& a_facesFine,
                      const FaceIndex&         a_faceCoar);
  void
  coarsenAreaFrac(Real& a_areaFracCoar,
                  const Vector<Real>& a_areaFracFine);

  void
  coarsenVolFracAndCentroid(Real&                   a_volFracCoar,
                            RealVect&               a_volCentroidCoar,
                            const Vector<Real>&     a_volFracFine,
                            const Vector<RealVect>& a_volCentroidFine,
                            const Vector<VolIndex>& a_fineVoFs,
                            const VolIndex&         a_coarVoF);

  void
  coarsenBoundaryAreaAndNormal(Real&                    a_bndryAreaCoar,
                               RealVect&                a_normalCoar,
                               const Vector<Real>&      a_bndryAreaFine,
                               const Vector<RealVect>&  a_normalFine);

  RealVect
  fineToCoarseTransform(const RealVect& a_finePoint,
                        const IntVect&  a_coarCell,
                        const IntVect&  a_fineCell);

  void
  coarsenBndryCentroid(RealVect&               a_bndryCentroidCoar,
                       const Vector<RealVect>& a_bndryCentroidFine,
                       const Vector<Real>&     a_bndryAreaFine,
                       const Vector<VolIndex>& a_fineVoFs,
                       const VolIndex&         a_coarVoF);

  void fetch(std::list<const VolData*>& fineVols, const Vector<VolIndex>& vofsFine) const;

  friend class EBIndexSpace;
  friend class EBISLevel;
};

///
/**
   Ref-counted version of EBDataImplem.
*/
class EBData
{
public:

  ///
  EBData(): m_implem( RefCountedPtr<EBDataImplem>( new EBDataImplem() ) )
    {
    }

  ///
  ~EBData()
    {
    }

  ///
  /**
   */
  void define(const Box& a_box, int a_comps)
    {
      m_implem->define(a_box, a_comps);
    }


  ///
  /**
   */
  EBData(const Box& a_box, int a_comps)
  :m_implem(new EBDataImplem(a_box, a_comps))
    {
    }

  ///
  /**
     Copy the information from a_source to the over the intersection
     of the box a_region, the box of the current EBData and
     the box of a_source.  The Interval arguments are ignored.
     This function is required by LevelData.
  */
  void copy(const Box&      a_regionFrom,
            const Interval& a_Cd,
            const Box&      a_regionto,
            const EBData&   a_source,
            const Interval& a_Cs)
    {
      m_implem->copy(a_regionFrom, a_Cd, a_regionto, *a_source.m_implem, a_Cs);
    }


  ///define the whole thing
  void
  define(const EBGraph&           a_graph,
         const Vector<IrregNode>& a_irregData,
         const Box&               a_validBox,
         const Real &             a_dx,
         bool                     a_hasMoments)
    {
      m_implem->define(a_graph, a_irregData, a_validBox, a_dx, a_hasMoments);
      computeNormalsAndBoundaryAreas(a_graph, a_validBox);
    }

  //this is the define to use when you are just going to copy the data into it.  region should include ghost cells
  void define(const EBGraph&           a_graph,
              const Box&               a_region,
              const Real&              a_dx,
              bool a_hasMoments) 
    {
      m_implem->define(a_graph, a_region, a_dx, a_hasMoments);
    }
  ///
  void
  coarsenVoFs(const EBData&   a_fineEBData,
              const EBGraph&  a_fineGraph,
              const EBGraph&  a_coarGraph,
              const Box&      a_validRegion)
    {
      m_implem->coarsenVoFs(*a_fineEBData.m_implem, a_fineGraph, a_coarGraph, a_validRegion);
    }

  void 
  coarsenFaces(const EBData& a_fineEBData,
               const EBGraph& a_fineGraph,
               const EBGraph& a_coarGraph,
               const Box&     a_validRegion)
    {
      m_implem->coarsenFaces(*a_fineEBData.m_implem, a_fineGraph, a_coarGraph, a_validRegion);
    }

  /// get the moment at  the VoF associated with the  monomial with the input exponents
  /**
     Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
     returns integral_over_VoF(x^p y^q z^r dV) for all p q
  **/
  IndMomSpaceDim getVolumeMoments(const VolIndex& a_vof) const
    {
      return m_implem->getVolumeMoments(a_vof);
    }

  /// get the normal the irregular face associated with the  monomial with the input exponents 
  /**
     Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
     returns integral_over_VoF((x^p y^q z^r) dV) for p q r
  **/
  IndMomSpaceDim getEBMoments(const VolIndex& a_vof) const
    {
      return m_implem->getEBMoments(a_vof);
    }

  /// get the moment at  the face associated with the  monomial with the input exponents
  /**
     Given face variables x, y,  p = mono(0), q = mono(1)),
     returns integral_over_face(x^p y^q dA) for all p q
  **/
  IndMomSDMinOne getFaceMoments(const FaceIndex& a_face) const
    {
      return m_implem->getFaceMoments(a_face);
    }

  /// I am sick of guessing this name wrong
  IndMomSDMinOne getAreaMoments(const FaceIndex& a_face) const
    {
      return getFaceMoments(a_face);
    }

  /// 
  /**
     gets the partial derivs of the normal component.   
  **/
  IndMomSpaceDim getEBNormalPartialDerivs(const VolIndex& a_vof, int a_normalComponent) const
    {
      return m_implem->getEBNormalPartialDerivs(a_vof, a_normalComponent);
    }

  /// get the normal*moment at  the irregular face associated with the  monomial with the input exponents 
  /**
     Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
     returns integral_over_VoF((x^p y^q z^r)*n_i dV) for p q r
     returns the OUTWARD (away from the center of the vof) normal comps
  **/
  IndMomSpaceDim getEBNormalMoments(const VolIndex& a_vof, int a_normalComponent) const
    {
      return m_implem->getEBNormalMoments(a_vof, a_normalComponent);
    }

  ///
  const Real& volFrac(const VolIndex& a_vof) const;

  ///
  const Real& areaFrac(const FaceIndex& a_face1) const;

  ///
  const RealVect& centroid(const FaceIndex& facein) const;

  ///
  const RealVect& centroid(const VolIndex& a_vof) const;

  ///return true if higher order moments are available
  bool hasMoments() const
    {
      return m_implem->hasMoments();
    }
  ///
  const RealVect& bndryCentroid(const VolIndex& a_vof) const;
  const RealVect& bndryCentroid(const VolIndex& a_vof, int face) const;

  ///
  const Real& bndryArea(const VolIndex& a_vof) const;
  const Real& bndryArea(const VolIndex& a_vof, int face) const;

  ///
  const RealVect& normal(const VolIndex& a_vof) const;
  const RealVect& normal(const VolIndex& a_vof, int face) const;

  /// used by multi-fluid code
  int facePhase(const VolIndex& a_vof, int a_face) const 
    {
      return m_implem->facePhase(a_vof, a_face);
    }

  /// used by multi-fluid code
  const VolIndex& faceIndex(const VolIndex& a_vof, int face) const 
    {
      return m_implem->faceIndex(a_vof, face);
    }

  /// used by multi-fluid code
  void setFacePhase(const VolIndex& a_vof, int face, int phase)
    {
      m_implem->setFacePhase(a_vof, face, phase);
    }

  /// used by multi-fluid code
  void setFaceIndex(const VolIndex& a_vof, int face, const VolIndex& index)
    {
      m_implem->setFaceIndex(a_vof, face, index);
    }
  ///
  int numFacePhase(const VolIndex& a_vof) const 
    {
      return m_implem->numFacePhase(a_vof);
    }

  void clearMultiBoundaries()
    {
      m_implem->clearMultiBoundaries();
    }
  void setBoundaryPhase(int phase)
    {
      m_implem->setBoundaryPhase(phase);
    }
  ///
  EBData(const EBData& a_ebiin)
    {
      m_implem = a_ebiin.m_implem;
    }

  ///
  /**
     This is a pointer comparison.
  */
  bool operator==(const EBData& a_ebiin)
    {
      return ((&(*m_implem)) == (&(*a_ebiin.m_implem)));
    }

  ///
  EBData& operator=(const EBData& a_ebiin)
    {
      m_implem = a_ebiin.m_implem;
      return *this;
    }

  ///
  /**
     This stuff required by LevelData in parallel:
  */
  int size(const Box& R, const Interval& comps) const
    {
      return  (m_implem->size(R, comps));
    }

  ///
  void linearOut(void* buf, const Box& R, const Interval& comps) const
    {
      m_implem->linearOut(buf, R, comps);
    }

  ///
  void linearIn(void* buf, const Box& R, const Interval& comps)
    {
      m_implem->linearIn(buf, R, comps);
    }

  ///
  void
  addFullIrregularVoFs(const IntVectSet& a_vofsToChange,
                       const EBGraph&    a_newGraph,
                       const BaseIVFAB<VolData>& a_grownData,
                       const EBGraph&    a_oldGraph)
    {
      m_implem->addFullIrregularVoFs(a_vofsToChange, a_newGraph, a_grownData, a_oldGraph);
    }


  ///multifluid angels dancing on the heads of pins.
  void
  addEmptyIrregularVoFs(const IntVectSet& a_vofsToChange,
                        const EBGraph&    a_newGraph)
    {
      m_implem->addEmptyIrregularVoFs(a_vofsToChange, a_newGraph);
    }


  ///
  void
  computeNormalsAndBoundaryAreas(const EBGraph& a_graph,
                                 const Box&     a_validRegion);

  BaseIVFAB<VolData>& getVolData()
    {
      return m_implem->getVolData();
    }

  const BaseIVFAB<VolData>& getVolData() const
    {
      return m_implem->getVolData();
    }

  ///
  static int preAllocatable()
    {
      return 2; // dyanmic allocatable.
    }


private:

  ///
  RefCountedPtr<EBDataImplem> m_implem;

  friend class EBIndexSpace;
  friend class EBISLevel;
};

/*******************************/
inline const Real& EBData::volFrac(const VolIndex& a_vof) const
{
  return m_implem->volFrac(a_vof);
}
/*******************************/
inline const Real& EBData::bndryArea(const VolIndex& a_vof, int face) const
{
  return m_implem->bndryArea(a_vof, face);
}
inline const Real& EBData::bndryArea(const VolIndex& a_vof) const
{
  return m_implem->bndryArea(a_vof);
}
/*******************************/
inline const RealVect& EBData::normal(const VolIndex& a_vof, int face) const
{
  return m_implem->normal(a_vof, face);
}
inline const RealVect& EBData::normal(const VolIndex& a_vof) const
{
  return m_implem->normal(a_vof);
}
/*******************************/
inline const RealVect& EBData::centroid(const VolIndex& a_vof) const
{
  return m_implem->centroid(a_vof);
}
/*******************************/
inline const RealVect& EBData::bndryCentroid(const VolIndex& a_vof, int face) const
{
  return m_implem->bndryCentroid(a_vof, face);
}
inline const RealVect& EBData::bndryCentroid(const VolIndex& a_vof) const
{
  return m_implem->bndryCentroid(a_vof);
}
/*******************************/
inline const RealVect& EBData::centroid(const FaceIndex& a_face) const
{
  return m_implem->centroid(a_face);
}
/*******************************/
inline const Real& EBData::areaFrac(const FaceIndex& a_face) const
{
  return m_implem->areaFrac(a_face);
}
/*******************************/

#include "NamespaceFooter.H"
#endif
