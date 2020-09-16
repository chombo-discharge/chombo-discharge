#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BCG_INTEGRATOR_H_
#define _BCG_INTEGRATOR_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"

#include "Proto.H"


using std::cout;
typedef Proto::Box Bx;
using   Proto::Point;
using   Proto::BoxData;
using   Proto::slice;
using   Proto::Stencil;
using     std::cout;
using     std::endl;
using std::string;
///
/**
 */
class BCG_Integrator
{
public:

  /// Default constructor.

  BCG_Integrator(const DisjointBoxLayout & a_grids,
                 const double            & a_dx,
                 const double            & a_viscosity,
                 int a_nghost):m_grids(a_grids), m_dx(a_dx), m_viscosity(a_viscosity), m_nghost(a_nghost)
  {
    m_exchangeCopier.exchangeDefine(a_grids, a_nghost*IntVect::Unit);
    defineStencils();
  }

  //divides cell centered velocity into vortical compoenent and a gradient of a scalar
  void
  ccProject(LevelData<FArrayBox>& a_velocity,
            LevelData<FArrayBox>& a_gradpres);

  ///
  /**
     Given velocity at n and gradient of pressure at n-1/2, advance both using BCG
   */
  void advanceSolution(LevelData<FArrayBox>& a_velocity, 
                       LevelData<FArrayBox>& a_gradpres, 
                       const double        & a_dt);

  ///
  void 
  computeVorticity(LevelData<FArrayBox>& vorticity,
                   LevelData<FArrayBox>& velocity);


  ///public for testing
  void
  getUDotDelU(LevelData<FArrayBox> & a_udelu,
              LevelData<FArrayBox> & a_velocity,
              const double & a_dt);

  ///  Destructor.
  ~BCG_Integrator(){};


private:
  //weak construction forbidden
  BCG_Integrator();  


  ///fills ghost cells of velocity with periodic image
  void
  enforceBoundaryConditions(LevelData<FArrayBox>& a_velocity);


  //divides face centered velocity into vortical compoenent and a gradient of a scalar
  void
  MACProject(LevelData<FluxBox>& a_velocity,
             LevelData<FluxBox>& a_gradpres);


  //computes nu*lapl(velcomp)
  void
  getSourceTerm(LevelData<FArrayBox> & a_source,
                LevelData<FArrayBox> & a_vel, 
                unsigned int a_comp);
  


  void
  averageVelocityToFaces(LevelData<FluxBox>   & a_velface,
                         LevelData<FArrayBox> & a_velcell);


  void
  MACDivergence(LevelData<FArrayBox> & a_divergence,
                LevelData<FluxBox>   & a_velface);
                 
  void
  MACGradient(LevelData<FluxBox>   &  a_macGrad,
              LevelData<FArrayBox> &  a_phicc);
                 
  void
  averageGradientToCell(LevelData<FArrayBox> & a_gradcell,
                        LevelData<FluxBox>   & a_gradface);
  void
  solveElliptic(LevelData<FArrayBox> & a_phi,
                LevelData<FArrayBox> & a_rhs,
                const double& a_alpha, const double& a_beta, 
                const string& a_basename);

  void
  getNuLaplU(BoxData<double, 1> & a_source, 
             const BoxData<double, 1> & a_scalarCell,
             const Bx           & a_grid,
             double a_coeff);

  DisjointBoxLayout  m_grids;
  double             m_dx;
  double             m_viscosity;
  int m_nghost;

  void  defineStencils();

  Stencil<double> m_velToFaceSten[DIM];
  Stencil<double> m_gradToCellSten[DIM];
  Stencil<double> m_macDivergeSten[DIM];
  Stencil<double> m_macGradientSten[DIM][DIM];
  Copier m_exchangeCopier;
#if DIM==2
  Stencil<double> m_vortsten[2];
#else
  Stencil<double> m_vortsten[3][3];
#endif

};
#endif //end include guard
