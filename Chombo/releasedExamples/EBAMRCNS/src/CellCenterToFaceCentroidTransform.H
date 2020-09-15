#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CELLCENTERTOFACECENTROIDTRANSFORM_H_
#define _CELLCENTERTOFACECENTROIDTRANSFORM_H_

#include "Stencils.H"
#include "EBIndexSpace.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "EBArith.H"
#include "VolIndex.H"
#include "RefCountedPtr.H"
#include "EBISBox.H"
#include "AggStencil.H"
#include "AggStencil.H"


///
/**
   Average data from cell centers to face centroids
*/
class CellCenterToFaceCentroidTransform
{
public:
  ///
  /**
  */
  CellCenterToFaceCentroidTransform(EBLevelGrid                              &  a_eblg,
                                    const IntVect                            &  a_ghostCellsCell,
                                    const IntVect                            &  a_ghostCellsFace,
                                    bool a_forceNoEBCF)
  {
    if(!a_forceNoEBCF)
      {
        //specifically, the coarse-fine interpolation will break
        MayDay::Error("this class does not deal with coarse-fine crossing EB");
      }
    m_eblg            =  a_eblg         ;   
    m_ghostCellsCell  =  a_ghostCellsCell;   
    m_ghostCellsFace  =  a_ghostCellsFace;   
    defineInternals();
  }

  ///
  virtual ~CellCenterToFaceCentroidTransform()
  {
  }

  ///
  /**
     get face centroid phi
     from cell center phi
     This version should go pretty fast as it uses aggstencil
  */
  void
  transform(LevelData<EBFluxFAB>        & a_phiFace,
            const LevelData<EBCellFAB>  & a_phiCell) const;
    
  ///
  /**
     get face centroid phi
     from cell center phi
     This is the debug version and should creep along
     but at least you can see things as they are going by.
  */
  void
  slowTransform(LevelData<EBFluxFAB>        & a_phiFace,
                const LevelData<EBCellFAB>  & a_phiCell) const;
    
    
protected:
  void defineInternals();
  EBLevelGrid   m_eblg;
  IntVect       m_ghostCellsCell;
  IntVect       m_ghostCellsFace;


  LayoutData<RefCountedPtr<AggStencil<EBCellFAB, EBFaceFAB> > > m_stencil[SpaceDim];


  void
  getStencils(Vector<RefCountedPtr< BaseStencil> >  & a_stencils, 
              Vector<RefCountedPtr< BaseIndex  > >  & a_baseDstFaces, 
              const Box                             & a_grid,
              const EBISBox                         & a_ebis,
              const DataIndex                       & a_dit,
              const int                             & a_idir) const;

  void
  getStencil(VoFStencil           & a_stencil,
             const FaceIndex      & a_face,
             const Box            & a_grid,
             const EBISBox        & a_ebis,
             const DataIndex      & a_dit) const;


private:
  
  //Weak construction is bad.   Even morally.
  CellCenterToFaceCentroidTransform()
  {
  }

  //lots of pointered data here.
  void operator=(const CellCenterToFaceCentroidTransform& a_opin)
  {
  }
  
  //lots of pointered data here.
  CellCenterToFaceCentroidTransform(const CellCenterToFaceCentroidTransform& a_opin)
  {
  }

};

#endif
