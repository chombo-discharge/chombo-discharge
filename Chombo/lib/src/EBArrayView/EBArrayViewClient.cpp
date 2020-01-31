#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef CH_SPACEDIM

// DTGraves, March, 2000

// Provisionally hacked upon by BJKeen, July 2001 to add some simple
// face centered data viewing.

#include <cmath>
#include <cstdio>
#include "EBArrayViewClient.H"
#include "LayoutIterator.H"
#include "VoFIterator.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "BoxIterator.H"
#include "DatasetClient.H"
#include "FaceIterator.H"
#include "EBIndexSpace.H"
#include "NamespaceHeader.H"

const int MAXBUFSIZEEB  = 1024;
const char *defaultFormatEB = "%7.5e";
const char *defaultLabelEB1 = "EBChombo MultiArrayView";
const char *defaultLabelEB2 = "EBChombo ArrayView";
// -------------------------------------------------------------------
// -------------------------------------------------------------------

bool
LDEBCellFABView(const LevelData<EBCellFAB >* const debugLevelData,
                const Box* const a_domainPtr)
{
  DisjointBoxLayout inputDBL = debugLevelData->getBoxes();
  const Box& a_domain = *a_domainPtr;
  bool noeekflag = true;
  const LevelData<EBCellFAB >& datain = *debugLevelData;
  DataIterator dit = datain.dataIterator();

  Vector<Box> boxes;
  LayoutIterator lit = inputDBL.layoutIterator();
  for (lit.reset(); lit.ok(); ++lit)
    boxes.push_back(inputDBL.get(lit()));
  //assign all boxes to proc 0
  Vector<int> assign(boxes.size(), 0);
  DisjointBoxLayout tempDBL(boxes, assign);
  tempDBL.close();
  const EBIndexSpace* const ebisPtr = Chombo_EBIS::instance();
  int nghost = 0;
  IntVect ghostvect = datain.ghostVect();
  for (int idir = 0; idir < SpaceDim; idir++)
    nghost = Max(nghost, ghostvect[idir]);
  EBISLayout ebisl;
  ebisPtr->fillEBISLayout(ebisl, tempDBL, a_domain, nghost);
  EBCellFactory ebfact(ebisl);
  LevelData<EBCellFAB > tempLDF(tempDBL,
                                datain.nComp(),
                                datain.ghostVect(),
                                ebfact);

  DataIterator dit2 = tempLDF.dataIterator();
  for (dit2.reset(); dit2.ok(); ++dit2)
    tempLDF[dit2()].setVal(0.);

  datain.copyTo(debugLevelData->interval(),
                tempLDF,
                tempLDF.interval());

  LevelData<EBCellFAB >* leveldata = &tempLDF;
  LayoutData<EBCellFAB >* layoutdata =
    static_cast<LayoutData<EBCellFAB >*>(leveldata);

  if (procID() == 0)
    noeekflag = MultiEBCellFABView(layoutdata);
  return noeekflag;
}
/************/
/************/
bool
IVFABView(const BaseIVFAB<Real>* IVFAB_ptr)
{
  const BaseIVFAB<Real>& ivfab = *IVFAB_ptr;
  const int ncomp = ivfab.nComp();

  // find maximum number of vofs per cell
  const IntVectSet& ivs = ivfab.getIVS();
  const EBGraph& ebgraph = ivfab.getEBGraph();

  int max_vofs = 0;
  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit)
    {
      IntVect iv = vofit().gridIndex();
      Vector<VolIndex> vofs = ebgraph.getVoFs(iv);
      int vofsize= vofs.size();
      max_vofs = Max(max_vofs, vofsize);
    }

  BaseFab<Real> fab(ivs.minBox(), max_vofs*ncomp+1);
  fab.setVal(0.0);
  for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      const IntVect& iv = vof.gridIndex();
      int min_offset = int(fab(iv,max_vofs*ncomp));
      fab(iv,max_vofs*ncomp) += 1.0;
      for (int comp = 0; comp < ncomp; ++comp)
        {
          for (int offset = min_offset; offset < max_vofs; ++offset)
            {
              fab(iv,max_vofs*comp+offset) = ivfab(vof,comp);
            }
        }
    }
  Real** data_ptr = new Real* [fab.nComp()];
  for (int comp = 0; comp < fab.nComp(); ++comp)
    {
      data_ptr[comp] = fab.dataPtr(comp);
    }
  return(IrregViewRealPtrArrayNVarDims(data_ptr, fab.nComp(),
                                       fab.box().loVect(),
                                       fab.box().hiVect(),
                                       defaultFormatEB, "BaseIVFAB"));
}
/************/
/************/
bool
EBCellFABView(const EBCellFAB* a_fab)
{
  const int nvar = a_fab->nComp();
  BaseFab<Real> fab(a_fab->getRegion(), nvar+1);
  makeFabFromEBFAB(fab, *a_fab);

  Real** data_ptr = new Real* [nvar+1];
  for (int var = 0; var <= nvar; ++var)
    {
      data_ptr[var] = fab.dataPtr(var);
    }
  //  return(IrregViewRealPtrArrayNVarDims(data_ptr, fab.nComp(),
  return(IrregViewRealPtrArrayNVarDims(data_ptr, fab.nComp(),
                                       fab.box().loVect(),
                                       fab.box().hiVect(),
                                       defaultFormatEB, defaultLabelEB2));
}
/************/
//this just puts the number of vofs in an extra variable
//at the end of the list.  This is probably not what we want
//to do but it's simple and I think the spreadsheet keys
//on that to leave blank stuff with zeros in the last var
// (of course you could never tell that from the source...)
/************/
void  makeFabFromEBFAB(BaseFab<Real>& fab,const EBCellFAB& ebfab)
{
  const BaseIVFAB<Real>* IVFAB_ptr = &(ebfab.getMultiValuedFAB());
  const BaseFab<Real>* fab_ptr = &(ebfab.getSingleValuedFAB());
  const BaseIVFAB<Real>& ivfab = *IVFAB_ptr;
  const int nvar = ivfab.nComp();
  CH_assert(fab.nComp() == nvar+1);
  CH_assert (fab_ptr->nComp() == nvar);

  const EBISBox& ebgraph = ebfab.getEBISBox();
  fab.setVal(1.0);
  fab.copy(*fab_ptr, 0, 0, nvar);
  BoxIterator bit(fab.box());
  for (bit.begin(); bit.ok(); bit.next())
    {
      IntVect iv = bit();
      if (ebgraph.isCovered(iv))
        {
          fab(iv,nvar) = 0.0;
        }
      else
        {
          Vector<VolIndex> vecvof = ebgraph.getVoFs(iv);
          Real val = Real(vecvof.size());
          fab(iv,nvar) = val;
        }
    }
  //this will just use the last vof's value
  //in a multiply-valued cell.
  for (int var = 0; var < nvar; ++var)
    {
      const IntVectSet ivs = ivfab.getIVS();
      for (VoFIterator vofit(ivs, ebgraph.getEBGraph()); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          const IntVect& iv = vof.gridIndex();
          fab(iv,var) = ivfab(vof,var);
        }
    }
}
/************/
/************/
bool
IFFABView(const BaseIFFAB<Real>* IFFAB_ptr)
{
  //the old one was buggy for cell faces with
  //multiple faces and we no longer have faceiterator
  //so i am just going to trash
  //it until i need it.

  // I need it, so I'm hacking something up based on ivfabview.
  // It's kludgiferous, since it's done w/o a FaceIterator; once
  // this is ready we shall need a rewrite; currently it's done in the
  // 2x work mode -BJK.

  //dodge on,  ben
  //cout << "iffabview dodgily implemented by student!" << endl;

  const BaseIFFAB<Real>& iffab = *IFFAB_ptr;
  const int ncomp = iffab.nComp();
  const int dir = iffab.direction();
  const EBGraph& ebbox = iffab.getEBGraph();

  // Find the max number of faces to be associated to any gridindexwise direction.
  int max_faces = 0;
  IVSIterator ivit(iffab.getIVS());
  for ( ivit.reset(); ivit.ok(); ++ivit )
    {
      SideIterator sit;
      for ( sit.reset(); sit.ok(); ++sit )
        {
          Vector<FaceIndex> myfaces = ebbox.getAllFaces(ivit(), dir, sit() );
          int cursize = myfaces.size();
          max_faces = Max( max_faces, cursize );
        }
    }

  // Grow out the box to consistent-w-chomboly
  Box mybox = iffab.getIVS().minBox();
  mybox.growHi( dir, 1 );

  BaseFab<Real> fab(mybox, max_faces*ncomp+1);
  fab.setVal(0.0);

  // Assuming that the format for this fab is that the last component
  // of each space tells how many faces/whatnot live in that cell,
  // this seems to be the right thing to do.

  // The IVFABView code I'm copying from for this is very undercommented.

  for ( ivit.reset(); ivit.ok(); ++ivit )
    {
      SideIterator sit;
      for ( sit.reset(); sit.ok(); ++sit )
        {
          Vector<FaceIndex> faces = ebbox.getAllFaces( ivit(), dir, sit() );

          IntVect targiv;
          if ( Side::Hi == sit() )
            targiv = ivit() + BASISV(dir);
          else
            targiv = ivit();

          fab(targiv, max_faces*ncomp) = (int)(faces.size());

          // This part is admittedly quite weird. I'm writing it to
          // coincide with what IVFABView seems to do. What it does,
          // too, it seems to do in a weird way.

          for ( int k=0; k<faces.size(); k++ )
            for ( int comp=0; comp < ncomp; comp++ )
              for ( int j=k; j < max_faces; j++ )
                fab(targiv, max_faces*comp + j) = iffab(faces[k], comp);
        }
    }

  Real **data_ptr = new Real *[fab.nComp()];
  for (int comp = 0; comp < fab.nComp(); ++comp)
    {
      data_ptr[comp] = fab.dataPtr(comp);
    }
  return(IrregViewRealPtrArrayNVarDims(data_ptr, fab.nComp(),
                                       fab.box().loVect(),
                                       fab.box().hiVect(),
                                       defaultFormatEB, "BaseIFFAB"));
}
/************/
/************/
bool EBFaceFABView(const EBFaceFAB* a_fab)
{
  //cheese on,  ben
  //cout << "ebfacefabview cheezily implemented" << endl;
  // basically just make the corresponding IFFAB

  EBISBox ebbox = a_fab->getEBISBox();
  IntVectSet myivs(a_fab->getCellRegion());
  BaseIFFAB<Real> chzfab( myivs, ebbox.getEBGraph(), a_fab->direction(),
                          a_fab->nComp() );

  FaceIterator fit( myivs, ebbox.getEBGraph(), a_fab->direction(),
                    FaceStop::SurroundingWithBoundary );

  // copy over.
  for ( fit.reset(); fit.ok(); ++fit )
    {
      for ( int k=0; k<a_fab->nComp(); k++ )
        chzfab(fit(),k) = (*a_fab)(fit(), k);
    }

  return IFFABView( &chzfab );
}

// -------------------------------------------------------------------
// This complete voodoo (everything below this comment)
//  brought to you by Vince Beckner and Dave Modiano.
// -------------------------------------------------------------------
bool
IrregViewRealPtrArrayNVarDims(Real *data[], int nvar,    // size nvar
                              const int *lodim, const int *hidim,  // size BL_SPACEDIM
                              const char *format, const char *label)
{
  int sockfd;

  if ( ! CreateSocket(sockfd))
    {
      return false;
    }

  // --------------------------------------------------- send data label
  if ( ! SendString(sockfd, label))
    {
      return false;
    }

  // --------------------------------------------------- send format
  if ( ! SendString(sockfd, format))
    {
      return false;
    }

  // --------------------------------------------------- send isIrregular
  if (!SendString(sockfd, "true")) // is irregular
    {
      return false;
    }

  // --------------------------------------------------- send isMultiFab
  if ( ! SendString(sockfd, "false"))
  {  // not a MultiFab
    return false;
  }

  // --------------------------------------------------- send nElements
  // dont send nElements

  // --------------------------------------------------- send the data
  return (SendRealArray(sockfd, data, nvar, lodim, hidim));

}  // end of function

// -------------------------------------------------------------------
bool
MultiEBCellFABView(const LayoutData<EBCellFAB >* layoutdata)
{
  const char *format= defaultFormatEB;
  const char *label= defaultLabelEB1;
  int  sockfd;
  char buffer[MAXBUFSIZEEB];

  if ( ! CreateSocket(sockfd))
    {
      return false;
    }

  // --------------------------------------------------- send data label
  if ( ! SendString(sockfd, label))
    {
      return false;
    }

  // --------------------------------------------------- send format
  if (format == NULL)
    {
      if ( ! SendString(sockfd, defaultFormatEB))
        {
          return false;
        }
    }
  else
    {
    if ( ! SendString(sockfd, format))
      {
        return false;
      }
    }

  // --------------------------------------------------- send isIrregular
  if (!SendString(sockfd, "true")) // is  irregular
    {
      return false;
    }

  // --------------------------------------------------- send isMulti
  if ( ! SendString(sockfd, "true"))
    {  // this has multiple grids
      return false;
    }

  // --------------------------------------------------- send nElements
  int num_elems = layoutdata->boxLayout().numBoxes(procID());
  //  cout << ">>> sending nElements = " << num_elems << endl;
  sprintf(buffer, "%d", num_elems);
  if ( ! SendString(sockfd, buffer))
    {
      return false;
    }

  //ArrayViewData data(layoutdata);

  // --------------------------------------------------- send the data
  for (DataIterator it = layoutdata->dataIterator(); it.ok(); ++it)
    {
      // construct dataArray for this element
      const EBCellFAB& ebfab = layoutdata->operator[](it());
      int nvar = ebfab.nComp();
      int fabvar = nvar+1;
      BaseFab<Real>  fab(ebfab.getRegion(), fabvar);
      makeFabFromEBFAB(fab, ebfab);

      Real** dataArray = new Real * [fabvar];
      for (int d = 0; d < fabvar; d++)
        {  // build the array of Real *
          dataArray[d] = fab.dataPtr(d);  // dont assume contiguous
        }
      int lo_vect[SpaceDim];
      int hi_vect[SpaceDim];
      for (int d = 0; d < SpaceDim; ++d)
        {
          lo_vect[d] = fab.box().smallEnd(d);
          hi_vect[d] = fab.box().bigEnd(d);
        }
      /*
        if ( ! SendRealArray(sockfd, dataArray, nvar,
        (fab.box()).loVect(), (fab.box()).hiVect()))
      */
      if ( ! SendRealArray(sockfd, dataArray, fabvar, lo_vect, hi_vect) )
        {
          return false;
        }
      delete [] dataArray;
    }

  return true;

}  // end of function
#include "NamespaceFooter.H"

#endif // CH_SPACEDIM
