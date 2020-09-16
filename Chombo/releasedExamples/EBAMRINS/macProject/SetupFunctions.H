#ifndef __SetupFunctionssss__
#define __SetupFunctionssss__
#define PI 3.141592653589793
#include "NamespaceHeader.H"

#if CH_SPACEDIM==2
IntVect lodeb(54, 14);
IntVect hideb(56, 17);
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

    retval = shared_ptr<BaseIF>(new SphereIF(a_geomRad, a_geomCen*RealVect::Unit, true));
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
                    int              & a_whichGeom,
                    int              & a_nx)
{
  pout() << "defining geometry" << endl;

  ParmParse pp;
  int maxGrid = 32;
    
  pp.get("nx"        , a_nx);
  pp.get("maxGrid", maxGrid);

  pout() << "nx       = " << a_nx     << endl;
  pout() << "maxGrid  = " << maxGrid  << endl;

  makeGrids(a_grids, a_dx, a_nx, maxGrid);
  Box domain = a_grids.physDomain().domainBox();

  RealVect origin = RealVect::Zero;

  pout() << "creating implicit function" << endl;
  shared_ptr<BaseIF>  impfunc = getImplicitFunction(a_geomCen, a_geomRad, a_whichGeom);

  RefCountedPtr<BaseIF> refcp(&(*impfunc));
  refcp.neverDelete();
  EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
  WrappedGShop workshop(refcp, RealVect::Zero, a_dx, a_grids.physDomain(), 0, 0);
  ebisPtr->define(domain, RealVect::Zero, a_dx, workshop, 1024, -1);
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
#define DIM CH_SPACEDIM
//=================================================

void InitializeVFace(Real&        a_vel,
                     IntVect      a_p,
                     Real         a_geomCen,
                     Real         a_geomRad,
                     Real         a_maxVelMag,
                     Real         a_maxVelRad,
                     Real         a_dx,
                     unsigned int a_dir)
{
  Real rlocsq = 0;
  Real xrel[DIM];
  for(int idir = 0; idir < DIM; idir++)
  {
    xrel[idir] = (a_p[idir] + 0.5)*a_dx - a_geomCen;
    rlocsq += xrel[idir]*xrel[idir];
  }

  Real velmag = sin(PI*rlocsq);
  Real vel[DIM];
#if DIM==2
  vel[0] =  velmag*xrel[1];
  vel[1] = -velmag*xrel[0];
#else
  vel[0] =  velmag*( xrel[1] + xrel[2]);
  vel[1] =  velmag*(-xrel[0] - xrel[2]);
  vel[2] = -velmag*( xrel[0] - xrel[1]);
#endif
  a_vel = vel[a_dir];
}
//=================================================

void InitializeVel(LevelData<EBFluxFAB>        &  a_fcVel,
                   const EBLevelGrid           &  a_eblg,
                   const Real                  &  a_dx,
                   const Real                  &  a_geomCen,
                   const Real                  &  a_geomRad,
                   const Real                  &  a_maxVelMag,
                   const Real                  &  a_maxVelRad)
{

  pout() << "setting intitial velocity " << endl;
  DisjointBoxLayout grids = a_eblg.getDBL();
  EBISLayout        ebisl = a_eblg.getEBISL();
  DataIterator dit = grids.dataIterator();
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box grid = grids[dit[ibox]];
    IntVectSet ivs(grid);
    EBFluxFAB& fcVel = a_fcVel[dit[ibox]];

    for(int idir = 0; idir < SpaceDim; idir++)
    {
      fcVel[idir].setVal(0);
      EBGraph ebgraph = ebisl[dit[ibox]].getEBGraph();

      for(FaceIterator faceit(ivs, ebgraph, idir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
      {
        IntVect iv = faceit().gridIndex(Side::Hi);
        Real velval;
        InitializeVFace(velval,
                        iv,
                        a_geomCen,
                        a_geomRad,
                        a_maxVelMag,
                        a_maxVelRad,
                        a_dx,
                        idir);
        fcVel[idir](faceit(), 0) = velval;
        
      }
    }
  }
}


#include "NamespaceFooter.H"

#endif
