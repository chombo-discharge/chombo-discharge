#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Stencils.H"
#include "EBIndexSpace.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "EBFluxFactory.H"
#include "EBArith.H"
#include "VolIndex.H"
#include "RefCountedPtr.H"
#include "EBISBox.H"
#include "NeumannConductivityDomainBC.H"
#include "RegTransformF_F.H"
#include "CellCenterToFaceCentroidTransform.H"
/////
void
CellCenterToFaceCentroidTransform::
transform(LevelData<EBFluxFAB>        & a_phiFace,
          const LevelData<EBCellFAB>  & a_phiCell) const
{
  CH_assert(a_phiFace.ghostVect() == m_ghostCellsFace);
  CH_assert(a_phiCell.ghostVect() == m_ghostCellsCell);

  LevelData<EBCellFAB>& phicast = (LevelData<EBCellFAB>&)a_phiCell;
  phicast.exchange();
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB &  cellEBCF = phicast[dit()];
      for(int idir = 0; idir < SpaceDim; idir++)
        {

          int isrc = 0; int idst = 0; int nco = a_phiCell.nComp();
          EBFaceFAB &  faceEBFF = a_phiFace[dit()][idir];
          bool incrOnly = false;
          m_stencil[idir][dit()]->apply(faceEBFF, cellEBCF, isrc, idst, nco, incrOnly);
        }
    }
}
/////
void
CellCenterToFaceCentroidTransform::
slowTransform(LevelData<EBFluxFAB>        & a_phiFace,
              const LevelData<EBCellFAB>  & a_phiCell) const
{
  CH_assert(a_phiFace.ghostVect() == m_ghostCellsFace);
  CH_assert(a_phiCell.ghostVect() == m_ghostCellsCell);

  LevelData<EBCellFAB>& phicast = (LevelData<EBCellFAB>&)a_phiCell;
  phicast.exchange();
  IntVect ivdebug(D_DECL(16,128,0));
  int istop = 0;
  for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebis = m_eblg.getEBISL()[dit()];
      const Box & valid = m_eblg.getDBL()[dit()];
      if(valid.contains(ivdebug))
        {
          istop = 1;
        }
      else
        {
          istop = 0;
        }
      EBCellFAB &  cellEBCF = phicast[dit()];
      IntVectSet  ivsGrid(valid);
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          FaceIterator faceit(ivsGrid, ebis.getEBGraph(), idir, FaceStop::SurroundingWithBoundary);
          const Vector<FaceIndex>& facevec = faceit.getVector();
          for(int ivec = 0; ivec < facevec.size(); ivec++)
            {
              const FaceIndex& face = facevec[ivec];
              const IntVect ivlo = face.getVoF(Side::Lo).gridIndex();
              const IntVect ivhi = face.getVoF(Side::Hi).gridIndex();
              if((idir== 1)&&((ivlo == ivdebug) || (ivhi == ivdebug)))
                {
                  istop = 1;
                }
              VoFStencil sten;
              getStencil(sten,  facevec[ivec], valid, ebis, dit());
              for(int icomp = 0; icomp < cellEBCF.nComp(); icomp++)
                {
                  sten.setAllVariables(icomp);
                  Real faceVal = applyVoFStencil(sten, cellEBCF, icomp);
                  a_phiFace[dit()][idir](facevec[ivec], icomp) = faceVal;
                }
            }
        }
    }
}
/////
void
CellCenterToFaceCentroidTransform::
defineInternals()
{
  EBFluxFactory fluxFact(m_eblg.getEBISL());
  EBCellFactory cellFact(m_eblg.getEBISL());
  LevelData<EBCellFAB> cellPhiProxy(m_eblg.getDBL(), 1, m_ghostCellsCell, cellFact);
  LevelData<EBFluxFAB>    fluxProxy(m_eblg.getDBL(), 1, m_ghostCellsFace, fluxFact);

  for(int idir = 0; idir < SpaceDim; idir++)
    {
      m_stencil[idir].define(m_eblg.getDBL());
      for(DataIterator dit = m_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
        {
          const     Box& grid = m_eblg.getDBL()  [dit()];
          const EBISBox& ebis = m_eblg.getEBISL()[dit()];
          Vector< RefCountedPtr<BaseIndex  > > baseDstFaces;
          Vector< RefCountedPtr<BaseStencil> > baseSten;
          getStencils(baseSten, baseDstFaces, grid, ebis,  dit(), idir);

          EBCellFAB& cellProxy = cellPhiProxy[dit()];
          EBFaceFAB& faceProxy = fluxProxy[dit()][idir];
          m_stencil[idir][dit()] = RefCountedPtr<AggStencil <EBCellFAB, EBFaceFAB>  >
            (new AggStencil<EBCellFAB, EBFaceFAB >(baseDstFaces, baseSten, cellProxy, faceProxy));
        }
    }
}
/////
void
CellCenterToFaceCentroidTransform::
getStencils(Vector<RefCountedPtr< BaseStencil> >  & a_stencils, 
            Vector<RefCountedPtr< BaseIndex  > >  & a_baseDstFaces, 
            const Box                             & a_grid,
            const EBISBox                         & a_ebis,
            const DataIndex                       & a_dit,
            const int                             & a_idir) const 
{
  // to keep ghost cell values sane
  IntVectSet  ivsGrid(a_grid);
  FaceIterator faceit(ivsGrid, a_ebis.getEBGraph(), a_idir, FaceStop::SurroundingWithBoundary);
  const Vector<FaceIndex>& facevec = faceit.getVector();
  a_baseDstFaces.resize(facevec.size());
  a_stencils.resize(    facevec.size());
  for(int ivec = 0; ivec < facevec.size(); ivec++)
    {
      VoFStencil sten;
      getStencil(sten,  facevec[ivec], a_grid, a_ebis, a_dit);
      a_stencils    [ivec]  = RefCountedPtr<BaseStencil>(new VoFStencil(sten));
      a_baseDstFaces[ivec]  = RefCountedPtr<BaseIndex  >(new FaceIndex(facevec[ivec]));
    }
}
/////
void
CellCenterToFaceCentroidTransform::
getStencil(VoFStencil           & a_stencil,
           const FaceIndex      & a_face,
           const Box            & a_grid,
           const EBISBox        & a_ebis,
           const DataIndex      & a_dit) const
{
  //other things will break here if CF/EB crosses, notably the coarse-fine interpolation.
  //Still, since it costs nothing, I might as well include the correct cfivs
  const IntVectSet& cfivs = (*m_eblg.getCFIVS())[a_dit];
  //these are the faces and weights to get from the face center to the face centroid
  FaceStencil centroidSten = EBArith::getInterpStencil(a_face, cfivs, a_ebis, m_eblg.getDomain());
  Box domBox = m_eblg.getDomain().domainBox();
  a_stencil.clear();
  for(int isten = 0; isten < centroidSten.size(); isten++)
    {
      const FaceIndex& face   = centroidSten.face(isten);
      const Real     & facewt = centroidSten.weight(isten);
      //sum up the face-centered stencils*weights
      IntVect ivlo = face.getVoF(Side::Lo).gridIndex();
      IntVect ivhi = face.getVoF(Side::Hi).gridIndex();
      if(domBox.contains(ivlo) && a_ebis.isCovered(ivlo))
        {
          MayDay::Error("face input inconsistent with geometry i was orginally given");
        }
      if(domBox.contains(ivhi) && a_ebis.isCovered(ivhi))
        {
          MayDay::Error("face input inconsistent with geometry i was orginally given");
        }
      bool uselo = domBox.contains(ivlo) && (!cfivs.contains(ivlo));
      bool usehi = domBox.contains(ivhi) && (!cfivs.contains(ivhi));
      if((!uselo) && (!usehi))
        {
          MayDay::Error("construction error in transform");
        }
      VoFStencil faceCentStencil;
      if(uselo && usehi)
        {
          faceCentStencil.add(face.getVoF(Side::Lo), 0.5);
          faceCentStencil.add(face.getVoF(Side::Hi), 0.5);
        }
      else if(uselo)
        {
          faceCentStencil.add(face.getVoF(Side::Lo), 1.0);
        }
      else if(usehi)
        {
          faceCentStencil.add(face.getVoF(Side::Hi), 1.0); 
        }
      else
        {
          MayDay::Error("logic error");
        }
      faceCentStencil *= facewt;
      a_stencil += faceCentStencil;
    }
                                   
}
