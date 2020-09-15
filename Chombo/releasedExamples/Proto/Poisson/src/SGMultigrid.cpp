#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SGMultigrid.H"
#include "Proto.H"
#include "ProtoInterface.H"
#include "BoxLayout.H"

int  SGMultigrid::s_numSmoothDown   = 2;
int  SGMultigrid::s_numSmoothUp     = 2;
int  SGMultigrid::s_numSmoothBottom = 16;
typedef Proto::Box Bx;
using Proto::BoxData;
using Proto::Stencil;
using Proto::Shift;
/****/
void
SGMultigridLevel::
applyOp(LevelData<FArrayBox >       & a_lph,
        const LevelData<FArrayBox > & a_phi)
                    
{
  enforceBoundaryConditions(a_phi);
  {
  CH_TIME("sgmg::applyop_no_comm");
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel for
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Bx grid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    BoxData<double,1> lphbd, phibd;
    ProtoCh::aliasBoxData<double, 1>(lphbd, a_lph[dit[ibox]]);
    ProtoCh::aliasBoxData<double, 1>(phibd, a_phi[dit[ibox]]);

    lphbd |= m_negoperator(phibd,   grid);
    lphbd *= -1.0;
  }
  ch_flops() += 3*m_numPointsThisProc*m_negoperator.size();
  }
}
/****/
void
SGMultigrid::
applyOp(LevelData<FArrayBox >       & a_lph,
        const LevelData<FArrayBox > & a_phi)
{
  return m_finest->applyOp(a_lph, a_phi);
}
/***/
SGMultigrid::
SGMultigrid(const double & a_alpha,
            const double & a_beta,
            const double & a_dx,
            const DisjointBoxLayout   & a_grids,
            bool a_useDenseStencil)
{
  m_finest = std::shared_ptr<SGMultigridLevel>(new SGMultigridLevel(a_alpha, a_beta, a_dx, a_grids, a_useDenseStencil));
}
/***/
void
SGMultigrid::
residual(LevelData<FArrayBox >       & a_res,
         const LevelData<FArrayBox > & a_phi,
         const LevelData<FArrayBox > & a_rhs)
{
  CH_TIME("sgmg::resid");
  return m_finest->residual(a_res, a_phi, a_rhs);
}
/***/
void
SGMultigrid::
vCycle(LevelData<FArrayBox >       & a_phi,
       const LevelData<FArrayBox > & a_rhs)
{
  CH_TIME("sgmg::vcycle");
  return m_finest->vCycle(a_phi, a_rhs);
}
/***/
void
SGMultigridLevel::
getMultiColors()
{
// Color offsets are grouped into "red"=even number of nonzeros (first 2^(DIM-1)) 
// and "black= odd number of nonzeros (the rest).
#if DIM==2
  s_colors[0] = Point::Zeros();//(0,0)
  s_colors[1] = Point::Ones();//(1,1)
  s_colors[2] = Point::Zeros() + Point::Basis(1);//(0,1)
  s_colors[3] = Point::Zeros() + Point::Basis(0);//(1,0)
#elif DIM==3
  s_colors[0] = Point::Zeros();//(0,0,0)
  s_colors[1] = Point::Zeros() + Point::Basis(0) + Point::Basis(1);//(1,1,0)
  s_colors[2] = Point::Zeros() + Point::Basis(1) + Point::Basis(2);//(0,1,1)
  s_colors[3] = Point::Zeros() + Point::Basis(0) + Point::Basis(2);//(1,0,1)
  s_colors[4] = Point::Zeros() + Point::Basis(1);//(0,1,0)
  s_colors[5] = Point::Zeros() + Point::Basis(0);//(1,0,0)
  s_colors[6] = Point::Zeros() + Point::Basis(2);//(0,0,1)
  s_colors[7] = Point::Ones();//(1,1,1)
#else
  compiler_error_this_code_is_only_written_for_dim_2_or_3();
#endif
}
/***/
SGMultigridLevel::
SGMultigridLevel(const double & a_alpha,
                 const double & a_beta,
                 const double & a_dx,
                 const DisjointBoxLayout & a_grids,
                 bool a_useDenseStencil)
{
  m_alpha     = a_alpha;
  m_beta      = a_beta;
  m_dx        = a_dx;
  m_grids     = a_grids;
  m_useDenseStencil = a_useDenseStencil;
  m_exchangeCopier.exchangeDefine(m_grids, IntVect::Unit);
  m_resid.define(m_grids,1, IntVect::Zero);
  m_numPointsThisProc = numPointsThisProc(a_grids);

  defineStencils();

  defineCoarserObjects();
}
/***/
void
SGMultigridLevel::
defineCoarserObjects()
{
  CH_TIME("sgmglevel::defineCoarser");
  if(m_grids.coarsenable(4))
  {
    DisjointBoxLayout gridsCoar;
    coarsen(gridsCoar,  m_grids, 2);
    m_residC.define(gridsCoar, 1, IntVect::Zero);
    m_deltaC.define(gridsCoar, 1, IntVect::Unit);
    m_coarser = std::shared_ptr<SGMultigridLevel>(new SGMultigridLevel(m_alpha, m_beta, 2*m_dx, gridsCoar, m_useDenseStencil));
    m_hasCoarser = true;
  }
  else
  {
    m_hasCoarser= false;
  }
}
/***/
void
SGMultigridLevel::
defineStencils()
{
  CH_TIME("sgmglevel::definestencils");
  getMultiColors();

  //always need -lapl(phi) so store that.
  if((m_useDenseStencil) && (DIM > 1) && (DIM < 4))
  {
#if DIM==2
    m_negoperator = (-m_alpha)*Shift(Point::Zeros()) + (-m_beta/(m_dx*m_dx))*(Stencil<double>::Laplacian_9());
#else
    m_negoperator = (-m_alpha)*Shift(Point::Zeros()) + (-m_beta/(m_dx*m_dx))*(Stencil<double>::Laplacian_27());
#endif
  }
  else
  {
    m_negoperator = (-m_alpha)*Shift(Point::Zeros()) + (-m_beta/(m_dx*m_dx))*(Stencil<double>::Laplacian());
  }
  
//  cout << "2/dx/dx = " << 2./m_dx/m_dx << ", 1/dx/dx=" << 1./m_dx/m_dx << ", neg operator = ";
//  m_negoperator.print();
  double safety = 1.0;
  double diag = m_alpha + (m_beta*(-2.*DIM)/(m_dx*m_dx));
  m_lambda = safety/diag;

      
  m_restrict = Stencil<double>();
  double numpts = double(MG_NUM_COLORS);
  for(int icolor = 0; icolor < MG_NUM_COLORS; icolor++)
  {
    m_relaxOpPhi[icolor] = (m_lambda)*m_negoperator;
    m_relaxOpPhi[icolor] *= (1.0)*Shift(s_colors[icolor]);

    m_relaxOpRhs[icolor] = (m_lambda)*Shift(s_colors[icolor]);
    m_updateOpPhi[icolor] = (1.0)*Shift(Point::Zeros());

    m_relaxOpPhi[icolor].destRatio() = Point::Ones(1);
    m_relaxOpPhi[icolor].srcRatio()  = Point::Ones(2);

    m_updateOpPhi[icolor].destRatio() = Point::Ones(2);
    m_updateOpPhi[icolor].srcRatio()  = Point::Ones(1);
    m_updateOpPhi[icolor].destShift() = s_colors[icolor];

    m_relaxOpRhs[icolor].destRatio() = Point::Ones(1);
    m_relaxOpRhs[icolor].srcRatio()  = Point::Ones(2);

    m_prolong[icolor]  =  (1.0)*Shift(Point::Zeros());
    m_prolong[icolor].destRatio() = Point::Ones(2);
    m_prolong[icolor].destShift() = s_colors[icolor];

    m_restrict += (1.0/numpts)*Shift(s_colors[icolor]);
  }    
  m_restrict.srcRatio() = Point::Ones(2);
}
/****/
void
SGMultigridLevel::
residual(LevelData<FArrayBox >       & a_res,
         const LevelData<FArrayBox > & a_phi,
         const LevelData<FArrayBox > & a_rhs)
                    
{
  enforceBoundaryConditions(a_phi);
  {
  CH_TIME("sgmglevel::residual");
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel for
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Bx grid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    BoxData<double, 1> resbd, phibd, rhsbd;
    ProtoCh::aliasBoxData<double, 1>(resbd, a_res[dit[ibox]]);
    ProtoCh::aliasBoxData<double, 1>(rhsbd, a_rhs[dit[ibox]]);
    ProtoCh::aliasBoxData<double, 1>(phibd, a_phi[dit[ibox]]);
    resbd |= m_negoperator(phibd, grid);
    resbd += rhsbd;
  }
  }
}
/****/
void
SGMultigridLevel::
relax(LevelData<FArrayBox >       & a_phi,
      const LevelData<FArrayBox > & a_rhs)
{
  CH_TIME("sgmglevel::relax");
  // GSRB. As implemented here, only correct for second-order 5 / 7 point operators. 
  // To go to higher order, need to use full multicolor algorithm. 
  DataIterator dit = m_grids.dataIterator();
  if(m_useDenseStencil)
  {
    
    for(int icolor =0; icolor < MG_NUM_COLORS; icolor++)
    {
      enforceBoundaryConditions(a_phi);
      {
        CH_TIME("smglevel::relax_one_multi_no_comm");
#pragma omp parallel for
        for(int ibox = 0; ibox < dit.size(); ibox++)
        {
          Bx grid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
          BoxData<double,1> rhsbd, phibd;
          ProtoCh::aliasBoxData<double, 1>(rhsbd, a_rhs[dit[ibox]]);
          ProtoCh::aliasBoxData<double, 1>(phibd, a_phi[dit[ibox]]);
          Box coarDom = grid.coarsen(2);
          BoxData<double,1> phisrc(coarDom);
          //needs to be coarse domain because of the whole gsrb thing            
          phisrc |= m_relaxOpPhi[icolor](phibd,   coarDom);
          phisrc += m_relaxOpRhs[icolor](rhsbd,   coarDom);
          phibd +=  m_updateOpPhi[icolor](phisrc,  coarDom);
        }
      }
    }
  }
  else
  {
    // Loop over red, black colors.
    for (int evenOdd=0;evenOdd < 2 ; evenOdd++)
    {
      enforceBoundaryConditions(a_phi);
      {
        CH_TIME("smglevel::relax_one_redblack_no_comm");
#pragma omp parallel for
        for(int ibox = 0; ibox < dit.size(); ibox++)
        {
          Bx grid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
          BoxData<double,1> rhsbd, phibd;
          ProtoCh::aliasBoxData<double, 1>(rhsbd, a_rhs[dit[ibox]]);
          ProtoCh::aliasBoxData<double, 1>(phibd, a_phi[dit[ibox]]);
          Box coarDom = grid.coarsen(2);
          BoxData<double,1> phisrc(coarDom);
          // loop over 2^(DIM-1) coarsened domains in each color. 
          for(int icolor = evenOdd*MG_NUM_COLORS/2; icolor < evenOdd*MG_NUM_COLORS/2 + MG_NUM_COLORS/2; icolor++)
          {
            //needs to be coarse domain because of the whole gsrb thing            
            phisrc |= m_relaxOpPhi[icolor](phibd,   coarDom);
            phisrc += m_relaxOpRhs[icolor](rhsbd,   coarDom);
            phibd +=  m_updateOpPhi[icolor](phisrc,  coarDom);
          }
        }
      }
    }
  }
  unsigned long long int stensize = m_relaxOpPhi.size() + m_relaxOpRhs.size() + m_updateOpPhi.size();
  ch_flops() += 3*m_numPointsThisProc*stensize;
}
/****/
void
SGMultigridLevel::
restrictResidual(LevelData<FArrayBox >       & a_resc,
                 const LevelData<FArrayBox > & a_res)
{
  CH_TIME("sgmglevel::restrict");
  DataIterator dit = m_grids.dataIterator();
#pragma omp parallel for
  for(int ibox = 0; ibox < dit.size(); ibox++)
  {
    Bx grid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
    BoxData<double,1> rescbd, resfbd;
    ProtoCh::aliasBoxData<double, 1>(rescbd, a_resc[dit[ibox]]);
    ProtoCh::aliasBoxData<double, 1>(resfbd,  a_res[dit[ibox]]);
    //called by the coarser mg level
    rescbd |= m_restrict(resfbd,  grid);
  }
  ch_flops() += 3*m_numPointsThisProc*(m_restrict.size());
}
/****/
void
SGMultigridLevel::
prolongIncrement(LevelData<FArrayBox >      & a_phi,
                 const LevelData<FArrayBox >& a_delta)
{
  CH_TIME("sgmglevel::prolong");
  //called by the coarser mg level
  DataIterator dit = m_grids.dataIterator();
  int nbox = dit.size();
  for(int icolor = 0; icolor < MG_NUM_COLORS; icolor++)
  {
#pragma omp parallel for
    for(int ibox = 0; ibox < nbox; ibox++)
    {
      Bx grid = ProtoCh::getProtoBox(m_grids[dit[ibox]]);
      BoxData<double,1> phibd, deltabd;
      ProtoCh::aliasBoxData<double, 1>(phibd  ,   a_phi[dit[ibox]]);
      ProtoCh::aliasBoxData<double, 1>(deltabd, a_delta[dit[ibox]]);
      phibd +=  m_prolong[icolor](deltabd, grid);
    }
  }
  //prolong size = 1 per fine grid point
  ch_flops() += 3*m_numPointsThisProc*MG_NUM_COLORS;

}
/****/
void 
SGMultigridLevel::
vCycle(LevelData<FArrayBox >         & a_phi,
       const LevelData<FArrayBox >   & a_rhs)
{
  CH_TIME("sgmglevel::vcycle");
  for(int irelax = 0; irelax < SGMultigrid::s_numSmoothDown; irelax++)
  {
    relax(a_phi,a_rhs); 
  }

  if (m_hasCoarser)
  {
    residual(m_resid,a_phi,a_rhs);                      
    m_coarser->restrictResidual(m_residC,m_resid);

    DataIterator dit = m_grids.dataIterator();
#pragma omp parallel for
    for(int ibox = 0; ibox < dit.size(); ibox++)
    {
      m_deltaC[dit[ibox]].setVal(0.);
    }

    m_coarser->vCycle(m_deltaC,m_residC);
    m_coarser->prolongIncrement(a_phi,m_deltaC);
  }
  else
  {
    for(int irelax = 0; irelax < SGMultigrid::s_numSmoothBottom; irelax++)
    {
      relax(a_phi,a_rhs);
    }
  }

  for(int irelax = 0; irelax < SGMultigrid::s_numSmoothUp; irelax++)
  {
    relax(a_phi,a_rhs);
  }

}
/****/
