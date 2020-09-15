#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BCG_Integrator.H"
#include "SGMultigrid.H"
#include "GodunovAdvectionOp.H"
#include "Proto.H"
#include "FluxBox.H"
#include "ProtoInterface.H"

using ProtoCh::aliasBoxData;
using ProtoCh::getProtoBox;
using Proto::GodunovAdvectionOp;
using Proto::Shift;
////
void
BCG_Integrator::
averageVelocityToFaces(LevelData<FluxBox>   & a_velface,
                       LevelData<FArrayBox> & a_velcell)
{
  CH_TIME("bcg:aveveltofaces");
  enforceBoundaryConditions(a_velcell);

  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box patch  = m_grids[dit[ibox]];
    Bx  cellbx = getProtoBox(patch);
    BoxData<double, DIM> velcellbd;
    aliasBoxData<double, DIM>(velcellbd, a_velcell[dit[ibox]]);

    for(int idir = 0; idir < DIM; idir++)
    {
      const Box& facebx =  cellbx.extrude(idir);
      BoxData<double, 1> velcomp  = slice(velcellbd, idir);
      BoxData<double, 1> velfacebd;

      aliasBoxData<double, 1>(velfacebd, a_velface[dit[ibox]][idir]);
      velfacebd |= m_velToFaceSten[idir](velcomp, facebx);
    }
  }
}
////
void
BCG_Integrator::
averageGradientToCell(LevelData<FArrayBox>   & a_gradcell,
                      LevelData<FluxBox>     & a_gradface)
{
  CH_TIME("bcg:avegradtocells");
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box patch  = m_grids[dit[ibox]];
    Bx  cellbx = getProtoBox(patch);
    BoxData<double, DIM> cellbd;
    aliasBoxData<double, DIM>(cellbd, a_gradcell[dit[ibox]]);
    for(int idir = 0; idir < DIM; idir++)
    {
      BoxData<double, DIM> facebd;
      aliasBoxData<double, DIM>(facebd, a_gradface[dit[ibox]][idir]);
      BoxData<double, 1> facecomp  = slice(facebd, idir);
      BoxData<double, 1> cellcomp  = slice(cellbd, idir);
      cellcomp  |= m_gradToCellSten[idir](facecomp, cellbx);
    }
  }
}
////
void
BCG_Integrator::
MACDivergence(LevelData<FArrayBox> & a_divergence,
              LevelData<FluxBox>   & a_velface)
{
  CH_TIME("bcg:macdiv");
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box patch  = m_grids[dit[ibox]];
    Bx  cellbx = getProtoBox(patch);
    BoxData<double, 1> divbd;
    aliasBoxData<double, 1>(divbd, a_divergence[dit[ibox]]);
    divbd.setVal(0.);
    for(int idir = 0; idir < DIM; idir++)
    {
      BoxData<double, 1> velbd;
      aliasBoxData<double, 1>(velbd, a_velface[dit[ibox]][idir]);
      divbd  += m_macDivergeSten[idir](velbd, cellbx);
    }
  }
}
////
void
BCG_Integrator::
MACGradient(LevelData<FluxBox>   &  a_macGrad,
            LevelData<FArrayBox> &  a_phicc)
{
  CH_TIME("bcg:macgrad");
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box patch  = m_grids[dit[ibox]];
    Bx  cellbx = getProtoBox(patch);

    BoxData<double, 1> phibd;
    aliasBoxData<double, 1>(phibd, a_phicc[dit[ibox]]);

    for(int faceDir = 0; faceDir < DIM; faceDir++)
    {
      BoxData<double, DIM> gradbd;
      aliasBoxData<double, DIM>(gradbd, a_macGrad[dit[ibox]][faceDir]);
      const Box& facebx =  cellbx.extrude(faceDir);
      
      for(int gradDir = 0; gradDir < DIM; gradDir++)
      {
        BoxData<double, 1> gradcomp  = slice(gradbd, gradDir);
        gradcomp |= m_macGradientSten[faceDir][gradDir](phibd, facebx);
      }
    }
  }
}

//divides face centered velocity into vortical compoenent and a gradient of a scalar
void
BCG_Integrator::
MACProject(LevelData<FluxBox> & a_velocity,
           LevelData<FluxBox> & a_gradpres)
{
  CH_TIME("bcg:macproject");
  double alpha = 0; double beta = 1; //solving Poisson
  LevelData<FArrayBox> divergence(m_grids, 1, IntVect::Zero);
  LevelData<FArrayBox>     scalar(m_grids, 1, m_nghost*IntVect::Unit);

  MACDivergence(divergence, a_velocity);
  double sumdiv = norm(divergence, Interval(0, 0), -3);
  std::cout << "sum of divergence = " << sumdiv << std::endl;

  solveElliptic(scalar, divergence, alpha, beta, string("projection:"));

  MACGradient(a_gradpres, scalar);
  
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    for(int idir = 0; idir < DIM; idir++)
    {
      a_velocity[dit[ibox]][idir] -= a_gradpres[dit[ibox]][idir];
    }  
  }  
}
//divides cell centered velocity into vortical compoenent and a gradient of a scalar
void
BCG_Integrator::
ccProject(LevelData<FArrayBox>& a_velocity,
          LevelData<FArrayBox>& a_gradpres)
{
  CH_TIME("bcg:ccproject");
  LevelData<FluxBox>  faceVel(m_grids,   1, IntVect::Zero);
  LevelData<FluxBox> faceGrad(m_grids, DIM, IntVect::Zero);
  
  averageVelocityToFaces(faceVel, a_velocity);
  MACProject(faceVel, faceGrad);
  averageGradientToCell(a_gradpres, faceGrad);
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    a_velocity[dit[ibox]] -= a_gradpres[dit[ibox]];
  }  
}
////
void
BCG_Integrator::
getNuLaplU(BoxData<double, 1> & a_source, 
           const BoxData<double, 1> & a_scalarCell,
           const Bx           & a_grid,
           double a_coeff)
{
  CH_TIME("bcg::nulapphi");
  //compue nu*lapl u
  Stencil<double> lap_operator = a_coeff*Stencil<double>::Laplacian();
  a_source |= lap_operator(a_scalarCell, a_grid);
}
////
void
BCG_Integrator::
getUDotDelU(LevelData<FArrayBox> & a_udelu,
            LevelData<FArrayBox> & a_velocity,
            const double         & a_dt)
{
  CH_TIME("bcg::getudelu");
  //using velocity averaged to faces as advection velocity
  //just used for upwinding
  LevelData<FluxBox>   advectVel(m_grids,   1, IntVect::Zero);
  LevelData<FluxBox>    faceGrad(m_grids, DIM, IntVect::Zero);
  LevelData<FluxBox>    faceVelo(m_grids, DIM, IntVect::Zero);

  averageVelocityToFaces(advectVel, a_velocity);

  GodunovAdvectionOp advOp(m_dx);
  a_velocity.exchange(Interval(0, DIM-1), m_exchangeCopier);
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box patch  = m_grids[dit[ibox]];
    Bx  cellbx = getProtoBox(patch);
    int doingvel = 1;
    BoxData<double, DIM> velobd, udelubd;
    aliasBoxData<double, DIM>(velobd , a_velocity[dit[ibox]]);
    aliasBoxData<double, DIM>(udelubd,    a_udelu[dit[ibox]]);
    for(int velComp = 0; velComp < DIM; velComp++)
    {
      BoxData<double,1>  scalarCell = slice(velobd, velComp);
      BoxData<double,1>  scalarFace[DIM];
      BoxData<double,1>  advectFace[DIM];
      for(int faceDir = 0; faceDir < DIM; faceDir++)
      {
        BoxData<double, DIM> facevelobd, facegradbd;
        aliasBoxData<double, DIM>(facevelobd , faceVelo[dit[ibox]][faceDir]);
        aliasBoxData<double, DIM>(facegradbd , faceGrad[dit[ibox]][faceDir]);
        aliasBoxData<double,   1>(advectFace[faceDir], advectVel[dit[ibox]][faceDir]);
        scalarFace[faceDir] = slice(facevelobd, velComp);
      }

      //needs to be grown by one because it has to exist on the other side
      //of the box face for Godunov algorithm to work
      Bx grownBx = cellbx.grow(1);
      BoxData<double, 1> source(cellbx);
      getNuLaplU(source, scalarCell, grownBx, m_viscosity);
      advOp.advectToFaces(scalarFace, advectFace, scalarCell, 
                          source, velobd, cellbx, doingvel, a_dt);
    }
  }
  
  //copy normal compoenents to advection velocity so we can mac project
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box patch  = m_grids[dit[ibox]];
    for(int faceDir = 0; faceDir < DIM; faceDir++)
    {
      BoxData<double, DIM> facevelobd;
      BoxData<double, 1  > advectvelbd;
      aliasBoxData<double, DIM>(facevelobd  ,  faceVelo[dit[ibox]][faceDir]);
      aliasBoxData<double,   1>(advectvelbd , advectVel[dit[ibox]][faceDir]);
      BoxData<double, 1> velcompbd = slice(facevelobd, faceDir);
      velcompbd.copyTo(advectvelbd);
    }
  }



  MACProject(advectVel, faceGrad);
  //now subtract the gradient of the scalar off all components 
  //(advectvel was corrected in mac project and just contains normal compoenents)
  //we use advectvel again below
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    for(int faceDir = 0; faceDir < DIM; faceDir++)
    {
      faceVelo[dit[ibox]][faceDir] -= faceGrad[dit[ibox]][faceDir];
    }
  }

//compute the divergence of the flux to get udelu
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box patch  = m_grids[dit[ibox]];
    Bx  cellbx = getProtoBox(patch);

    BoxData<double, DIM> udelubd;
    aliasBoxData<double, DIM>(udelubd , a_udelu[dit[ibox]]);
    for(int velDir = 0; velDir < DIM; velDir++)
    {
      BoxData<double,DIM> facevelobd[DIM];
      BoxData<double,1>      velcomp[DIM];
      BoxData<double,1>  advectvelbd[DIM];
      for(int faceDir = 0; faceDir < DIM; faceDir++)
      {
        aliasBoxData<double, DIM>( facevelobd[faceDir] ,  faceVelo[dit[ibox]][faceDir]);
        aliasBoxData<double,   1>(advectvelbd[faceDir] , advectVel[dit[ibox]][faceDir]);
        velcomp[faceDir] =   slice(facevelobd[faceDir], velDir);
      }

      BoxData<double, 1> divFcomp = slice(udelubd, velDir);
      advOp.divUPhi(divFcomp, advectvelbd, velcomp, cellbx);
    }
  }
}

////
void
BCG_Integrator::
advanceSolution(LevelData<FArrayBox>& a_velocity, 
                LevelData<FArrayBox>& a_gradpres, 
                const double        & a_dt)
{
  CH_TIME("bcg::advancesolution");
  //get advective derivative
  LevelData<FArrayBox> udelu(m_grids, DIM, IntVect::Zero);
  getUDotDelU(udelu, a_velocity, a_dt);

  LevelData<FArrayBox> ustar(m_grids, DIM, m_nghost*IntVect::Unit);
  DataIterator dit = m_grids.dataIterator();
  if(m_viscosity > 1.0e-16)
  {

    LevelData<FArrayBox> halfnulaplu(m_grids, DIM, IntVect::Zero);
    LevelData<FArrayBox>         rhs(m_grids, DIM, IntVect::Zero);

    
#pragma omp parallel
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      Box patch  = m_grids[dit[ibox]];
      Bx  cellbx = getProtoBox(patch);

      BoxData<double, DIM> velobd, laplbd, rhsbd, udelubd, gradpbd;
      aliasBoxData<double, DIM>(velobd ,  a_velocity[dit[ibox]]);
      aliasBoxData<double, DIM>(laplbd , halfnulaplu[dit[ibox]]);
      aliasBoxData<double, DIM>( rhsbd ,         rhs[dit[ibox]]);
      aliasBoxData<double, DIM>(udelubd,       udelu[dit[ibox]]);
      aliasBoxData<double, DIM>(gradpbd,  a_gradpres[dit[ibox]]);

      for(int velDir = 0; velDir < DIM; velDir++)
      {
        BoxData<double,1>  velcomp = slice(velobd, velDir);
        BoxData<double,1>  lapcomp = slice(laplbd, velDir);
        getNuLaplU(lapcomp, velcomp, cellbx, 0.5*m_viscosity);
      }

      //form c-n rhs
      rhsbd.setVal(0.);
      rhsbd -= udelubd;
      rhsbd -= gradpbd;
      rhsbd += laplbd;
      rhsbd *= a_dt;
      rhsbd += velobd;
    }

    //solve for u*
    double alpha = 1;
    double beta =  -0.5*m_viscosity*a_dt;
    for(int velDir = 0; velDir < DIM; velDir++)
    {
      AliasDataFactory<FArrayBox> velfact(&ustar,  Interval(velDir, velDir));
      AliasDataFactory<FArrayBox> rhsfact(&  rhs,  Interval(velDir, velDir));
      LevelData<FArrayBox> velscal(m_grids, 1, m_nghost*IntVect::Unit, velfact);
      LevelData<FArrayBox> rhsscal(m_grids, 1, m_nghost*IntVect::Unit, rhsfact);

      solveElliptic(velscal, rhsscal, alpha, beta, string("viscous_solve:"));
    }

  }
  else
  {
#pragma omp parallel
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      ustar[dit[ibox]].setVal(0.);
      ustar[dit[ibox]] -=      udelu[dit[ibox]];
      ustar[dit[ibox]] += a_gradpres[dit[ibox]];

      ustar[dit[ibox]] *= a_dt;
      ustar[dit[ibox]] += a_velocity[dit[ibox]];
    }
  }
  //now we need to project the solution onto its divergence-free subspace.
  //we include the usual tricks to improve stability
  //w = vel + dt*gph^n-1/2
  //u = P(w)
  //grad p^n+1/2 = (1/dt)(I-P)w;
  ccProject(ustar, a_gradpres);
  ustar.copyTo(a_velocity);
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    a_gradpres[dit[ibox]] /= a_dt;
  }

}
///
void
BCG_Integrator::
solveElliptic(LevelData<FArrayBox> & a_phi,
              LevelData<FArrayBox> & a_rhs,
              const   double       & a_alpha, 
              const   double       & a_beta,
              const string         & a_solvename)
{
  CH_TIME("bcg::solveelliptic");
  using std::cout;
  using std::endl;
  int numsmooth = 4;
  int maxiter = 27;
  double tol = 1.0e-9;
  SGMultigrid solver(a_alpha, a_beta, m_dx, m_grids);
  SGMultigrid::s_numSmoothUp   = numsmooth;
  SGMultigrid::s_numSmoothDown = numsmooth;

  LevelData<FArrayBox> res(m_grids, 1, IntVect::Zero);

  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    a_phi[dit[ibox]].setVal(0.);
  }

  int iter = 0;
  double rhsmax = norm(a_rhs, Interval(0,0), 0);

  double resStart = std::max(rhsmax, tol);
  double resIter  = rhsmax;
  double resIterOld = resIter;
  cout << a_solvename << "iter = " << iter << ", ||resid|| = " << resIter << endl;
  double mintol = 1.0e-10;
  double ratio = 2;
  double minrat = 1.05;
  while((resIter > mintol) && (resIter > tol*resStart) && (iter <  maxiter) && (ratio >minrat))
  {
    solver.vCycle(a_phi, a_rhs);
    solver.residual(res, a_phi, a_rhs);
  
    iter++;
    resIter = norm(res, Interval(0,0), 0);
    ratio = resIterOld/resIter;
    resIterOld = resIter;
    cout << "iter = " << iter << ", ||resid|| = " << resIter << endl;
  }
  a_phi.exchange(m_exchangeCopier);
}
////
void
BCG_Integrator::
defineStencils()
{
  CH_TIME("bcg::definestencils");
  double dx = m_dx;
  //stencils to average velocities from cells to faces and to increment the divergence 
  // of a mac velocity
  for(int faceDir = 0; faceDir < DIM; faceDir++)
  {
    m_velToFaceSten[faceDir]  = (0.5)*Shift::Zeros() +  (0.5)*Shift::Basis(faceDir,-1.0);    
    m_gradToCellSten[faceDir] = (0.5)*Shift::Zeros() +  (0.5)*Shift::Basis(faceDir, 1.0);    
    m_macDivergeSten[faceDir] = (-1.0/dx)*Shift::Zeros() + (1.0/dx)*Shift::Basis(faceDir, 1.0);    
    m_macGradientSten[faceDir][faceDir] = (1.0/dx)*Shift::Zeros() + (-1.0/dx)*Shift::Basis(faceDir, -1.0);    
    for(int gradDir = 0; gradDir < DIM; gradDir++)
    {
      if(faceDir != gradDir)
      {
        Point hihi = Point::Basis(gradDir, 1.0);
        Point hilo = Point::Basis(gradDir,-1.0);
        Point lolo = Point::Basis(gradDir,-1.0) + Point::Basis(faceDir, -1);
        Point lohi = Point::Basis(gradDir, 1.0) + Point::Basis(faceDir, -1);
        m_macGradientSten[faceDir][gradDir] = 
          ( 0.25/dx)*Shift(hihi) + ( 0.25/dx)*Shift(hilo)  +
          (-0.25/dx)*Shift(lolo) + (-0.25/dx)*Shift(lohi);
      }
    }
  }

///vorticity stencils
#if DIM==2
  m_vortsten[0] = (-0.5/dx)*Shift::Basis(1, 1.0) +  (0.5/dx)*Shift::Basis(1, -1.0);
  m_vortsten[1] =  (0.5/dx)*Shift::Basis(0, 1.0) + (-0.5/dx)*Shift::Basis(0, -1.0);
#else
  m_vortsten[0][2] = (-0.5/dx)*Shift::Basis(1, 1.0) +  (0.5/dx)*Shift::Basis(1, -1.0);
  m_vortsten[1][2] =  (0.5/dx)*Shift::Basis(0, 1.0) + (-0.5/dx)*Shift::Basis(0, -1.0);
  m_vortsten[2][2] = (0.0)*Shift::Zeros();

  m_vortsten[0][1] = (-0.5/dx)*Shift::Basis(2, 1.0) +  (0.5/dx)*Shift::Basis(2, -1.0);
  m_vortsten[2][1] =  (0.5/dx)*Shift::Basis(0, 1.0) + (-0.5/dx)*Shift::Basis(0, -1.0);
  m_vortsten[1][1] = (0.0)*Shift::Zeros();

  m_vortsten[2][0] =  (0.5/dx)*Shift::Basis(1, 1.0) + (-0.5/dx)*Shift::Basis(1, -1.0);
  m_vortsten[1][0] =  (0.5/dx)*Shift::Basis(2, 1.0) +  (0.5/dx)*Shift::Basis(2, -1.0);
  m_vortsten[0][0] = (0.0)*Shift::Zeros();
#endif
}

//cheerfully stolen from the euler example
void
BCG_Integrator::
enforceBoundaryConditions(LevelData<FArrayBox>& a_vel)
{

  CH_TIME("bcg::enforcebcs");
  a_vel.exchange(m_exchangeCopier);  
}
///
void 
BCG_Integrator::
computeVorticity(LevelData<FArrayBox>        & a_vorticity,
                 LevelData<FArrayBox>        & a_velocity)
{
  CH_TIME("bcg::vorticity");
  enforceBoundaryConditions(a_velocity);
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Box patch  = m_grids[dit[ibox]];
    Bx  cellbx = getProtoBox(patch);
    BoxData<double, DIM> velobd;
    aliasBoxData<double, DIM>(velobd ,  a_velocity[dit[ibox]]);
    a_vorticity[dit[ibox]].setVal(0.);
#if DIM==2
    BoxData<double, 1> vortbd;
    aliasBoxData<double, 1>(vortbd ,  a_vorticity[dit[ibox]]);
    for(int idir = 0; idir < DIM; idir++)
    {
      BoxData<double, 1> velcomp  = slice(velobd, idir);
      vortbd += m_vortsten[idir](velcomp, cellbx);
    }
#else
    BoxData<double, DIM> vortbd;
    aliasBoxData<double, DIM>(vortbd ,  a_vorticity[dit[ibox]]);
    for(int veldir = 0; veldir < DIM; veldir++)
    {
      for(int vortdir = 0; vortdir < DIM; vortdir++)
      {
        if(veldir != vortdir)
        {
          BoxData<double, 1> velcomp   = slice(velobd , veldir);
          BoxData<double, 1> vortcomp  = slice(vortbd , vortdir);
          vortcomp += m_vortsten[veldir][vortdir](velcomp, cellbx);
        }
      }
    }
#endif
  }//end of for loop
}
