#ifndef __SetupFunctionssss__
#define __SetupFunctionssss__
#define PI 3.141592653589793
#include "NamespaceHeader.H"

#if CH_SPACEDIM==2
IntVect lodeb(15, 15);
IntVect hideb(21, 21);
void 
dumpBlob(BaseFab<Real>* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(6)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    Box area(lodeb, hideb);

    BaseFab<Real> & data = *dataPtr;
    Box databox = dataPtr->box();
    cout << "data region contains:" << endl;
    for(int j = hideb[1]; j >= lodeb[1]; j--)
    {
      for(int i = lodeb[0]; i <= hideb[0]; i++)
      {
        IntVect pt(i,j);
        if(databox.contains(pt))
        {
          cout << pt << ":" << data(pt, 0) << "  ";
        }
      }
      cout << endl;
    }
  }
}
void 
dumpEB1(EBCellFAB* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(6)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    Box area(lodeb, hideb);

    EBGraph graph = dataPtr->getEBISBox().getEBGraph();
    EBCellFAB & data = *dataPtr;
    Box databox = dataPtr->box();
    cout << "data region contains:" << endl;
    for(int j = hideb[1]; j >= lodeb[1]; j--)
    {
      for(int i = lodeb[0]; i <= hideb[0]; i++)
      {
        IntVect pt(i,j);
        if(databox.contains(pt))
        {
          Vector<VolIndex > vofs = graph.getVoFs(pt);
          cout << pt << ":";
          if(vofs.size() >  0)
          {
            cout  << data(vofs[0], 0) << "  ";
          }
          else
          {
            cout  << "------------" << "  ";
          }
        }
      }
      cout << endl;
    }
  }
}




void 
dumpFace1(EBFaceFAB* dataPtr)
{
  if(dataPtr != NULL)
  {
    cout    << setprecision(6)
            << setiosflags(ios::showpoint)
            << setiosflags(ios::scientific);

    Box area(lodeb, hideb);

    EBFaceFAB & data = *dataPtr;
    EBGraph graph = data.getEBISBox().getEBGraph();
    Box databox = data.getCellRegion();
    cout << "data region contains:" << endl;
    for(int j = hideb[1]; j >= lodeb[1]; j--)
    {
      for(int i = lodeb[0]; i <= hideb[0]; i++)
      {
        IntVect pt(i,j);
        if(databox.contains(pt))
        {
          Vector<VolIndex> vofs = graph.getVoFs(pt);
          cout << pt << ":";
          bool printLines = false;
          if(vofs.size() >  0)
          {
            Vector< FaceIndex>faces = graph.getFaces(vofs[0], data.direction(), Side::Lo);
            if(faces.size() == 0)
            {
              printLines = true;
            }
            else
            {
              cout  << data(faces[0], 0) << "  ";
            }
          }
          else
          {
            printLines = true;
          }
          if(printLines)
          {
            cout  << "------------" << "  ";
          }
        }
      }
      cout << endl;
    }
  }
}
#endif
//=================================================

void makeGrids(DisjointBoxLayout& a_grids,
               Real             & a_dx,
               const int        & a_nx,
               const int        & a_maxGrid)
{
  pout() << "making grids" << endl;

  IntVect domLo = IntVect::Zero;
  IntVect domHi  = (a_nx - 1)*IntVect::Unit;

  // EB and periodic do not mix
  ProblemDomain domain(domLo, domHi);

  Vector<Box> boxes;
  unsigned int blockfactor = 8;
  domainSplit(domain, boxes, a_maxGrid, blockfactor);
  
  Vector<int> procs;

  a_dx = 1.0/a_nx;
  LoadBalance(procs, boxes);
  a_grids = DisjointBoxLayout(boxes, procs, domain);

}

///
shared_ptr<BaseIF>  getImplicitFunction(Real  & a_geomCen,
                                        Real  & a_geomRad,
                                        int   & a_whichGeom)

{
  shared_ptr<BaseIF>  retval;
  ParmParse pp;
  
  a_geomCen = 0;
  a_geomRad = 1;
  pp.get("which_geom", a_whichGeom);
  if(a_whichGeom == 0)
  {
    pout() << "sphere" << endl;

    pp.get("geom_cen", a_geomCen);
    pp.get("geom_rad", a_geomRad);
    pout() << "geom_cen = " << a_geomCen       << endl;
    pout() << "geom_rad = " << a_geomRad       << endl;

    retval = shared_ptr<BaseIF>(new SphereIF(a_geomRad, a_geomCen*RealVect::Unit, false));
  }
  else if(a_whichGeom ==  1)
  {
    pout() << "plane" << endl;
    RealVect normal, startPt;
    vector<double> v_norm, v_start;
    pp.getarr("geom_normal", v_norm, 0, SpaceDim);
    pp.getarr("geom_start_pt", v_start, 0, SpaceDim);
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      normal[ idir] = v_norm[ idir];
      startPt[idir] = v_start[idir];
      pout() << "normal ["<< idir << "] = " << normal [idir]  << endl;
      pout() << "startPt["<< idir << "] = " << startPt[idir]  << endl;
    }
    retval =  shared_ptr<BaseIF>(new PlaneIF(normal,startPt, true));
  }
  else
  {
    MayDay::Error("bogus geometry");
  }
  return retval;
}
//=================================================
void defineGeometry(DisjointBoxLayout& a_grids,
                    Real             & a_dx,
                    Real             & a_geomCen,
                    Real             & a_geomRad,
                    Real             & a_blobCen,
                    Real             & a_blobRad,
                    int              & a_whichGeom,
                    int              & a_nx)
{
  pout() << "defining geometry" << endl;

  ParmParse pp;
  int maxGrid = 32;
    
  pp.get("nx"        , a_nx);
  pp.get("maxGrid", maxGrid);
  pp.get("blob_cen", a_blobCen);
  pp.get("blob_rad", a_blobRad);

  pout() << "nx       = " << a_nx     << endl;
  pout() << "maxGrid  = " << maxGrid  << endl;
  pout() << "blob_cen = " << a_blobCen       << endl;
  pout() << "blob_rad = " << a_blobRad       << endl;

  makeGrids(a_grids, a_dx, a_nx, maxGrid);
  Box domain = a_grids.physDomain().domainBox();

  RealVect origin = RealVect::Zero;

  pout() << "creating implicit function" << endl;
  shared_ptr<BaseIF>  impfunc = getImplicitFunction(a_geomCen, a_geomRad, a_whichGeom);

  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  GeometryShop workshop(*impfunc, 0, a_dx*RealVect::Unit);
  ebisPtr->define(domain, RealVect::Zero, a_dx, workshop, 1024, 0);
}
//=================================================

Real InitializeLinePt(IntVect   a_p,
                      Real      a_px,
                      Real      a_rad,
                      Real      a_dx)
{
  Real val = 0;
  Real xloc[SpaceDim];
  for(int idir = 0; idir < SpaceDim; idir++)
  {
    xloc[idir] = (a_p[idir] + 0.5)*a_dx;
    
  }
  
  Real funcval = std::abs(a_px - xloc[0]);
  
  Real rloc = funcval;
  
  if(rloc < a_rad && (xloc[1] < 0.75) )
  {
    Real carg = 0.5*PI*(rloc)/a_rad;
    Real cosval = cos(carg);
    val = std::abs(cosval);
  }
  return val;                                          // 
}

//=================================================

void InitializeStuff(LevelData<EBCellFAB>        & a_scala,
                     LevelData<EBCellFAB>        & a_ccVel,
                     LevelData<EBFluxFAB>        & a_fcVel,
                     const DisjointBoxLayout     &  a_grids,
                     const Real                  &  a_dx,
                     const Real                  &  a_geomCen,
                     const Real                  &  a_geomRad,
                     const Real                  &  a_blobCen,
                     const Real                  &  a_blobRad,
                     const Real                  &  a_maxVelMag,
                     const Real                  &  a_maxVelRad)
{

  ParmParse pp;

  pout() << "setting constant velocity " << endl;
  pout() << "calling initializeline for scalar" << endl;
  vector<Real> geom_normal;
  pp.getarr("geom_normal", geom_normal, 0, SpaceDim);
  RealVect constVel = RealVect::Zero;
  constVel[0] =  geom_normal[1];
  constVel[1] = -geom_normal[0];
  DataIterator dit = a_grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    EBCellFAB& scala = a_scala[dit[ibox]];
    EBCellFAB& ccVel = a_ccVel[dit[ibox]];
    EBFluxFAB& fcVel = a_fcVel[dit[ibox]];
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      ccVel.setVal(idir, constVel[idir]);
      fcVel[idir].setVal(constVel[idir]);
    }
    scala.setVal(0.);
    Box grid = a_grids[dit[ibox]];
    EBGraph ebgraph = scala.getEBISBox().getEBGraph();
    IntVectSet ivs(grid);
    for(VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit)
    {
      IntVect iv = vofit().gridIndex();
      scala(vofit(), 0) = InitializeLinePt(iv, a_blobCen, a_blobRad, a_dx);
    }
  }
}


#include "NamespaceFooter.H"

#endif
