#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "BoxIterator.H"
#include "VoFIterator.H"

#include "EBArith.H"
#include "PolyGeom.H"

#include "InflowOutflowViscousTensorEBBC.H"
#include "CH_Timer.H"
#include "EBViscousTensorOp.H"
#include "EBAMRCNSParams.H"
#include "NamespaceHeader.H"

InflowOutflowViscousTensorEBBC::
InflowOutflowViscousTensorEBBC()
{
}


InflowOutflowViscousTensorEBBC::
InflowOutflowViscousTensorEBBC(const Real& a_dirichletStart,
                               const int & a_flowDir,
                               const ProblemDomain& a_domain,
                               const EBISLayout&    a_layout,
                               const RealVect&      a_dx):
  DirichletViscousTensorEBBC(a_domain, a_layout, a_dx, &IntVect::Zero, &IntVect::Zero)
{
  m_dirichletStart = a_dirichletStart;
  m_flowDir = a_flowDir;

  //this is to set the dirichlet bits correctly for when I use them
  setValue(0.);
}

InflowOutflowViscousTensorEBBC::
~InflowOutflowViscousTensorEBBC()
{
}
/*****/
void
InflowOutflowViscousTensorEBBC::
define(const LayoutData<IntVectSet>& a_cfivs,
       const Real&                   a_factor)
{
  DirichletViscousTensorEBBC::define(a_cfivs, a_factor);
  const DisjointBoxLayout& dbl = m_ebisl.getDisjointLayout();

  LayoutData<IntVectSet>   irregIVS;
  LayoutData<VoFIterator > vofItIrreg;
  vofItIrreg.define(dbl); // vofiterator cache
  irregIVS.define(dbl);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      irregIVS[dit()] = m_ebisl[dit()].getIrregIVS(dbl.get(dit()));
      vofItIrreg[dit()].define(irregIVS[dit()],m_ebisl[dit()].getEBGraph());
    }

  //the Dirichlet stencils have been made.  set to null if not sticky yet
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {

      for (vofItIrreg[dit()].reset(); vofItIrreg[dit()].ok(); ++vofItIrreg[dit()])
        {
          const VolIndex& vof = vofItIrreg[dit()]();
          RealVect normal = m_ebisl[dit()].normal(vof);
          RealVect bdCent = m_ebisl[dit()].bndryCentroid(vof);
          bdCent *= m_dx;
          RealVect point  = EBArith::getVofLocation(vof, m_dx, RealVect::Zero);
          point += bdCent;

          bool sticky = (point[m_flowDir] > m_dirichletStart);
          if(!sticky)
            {
              for (int ivar = 0; ivar < SpaceDim; ivar++)

                {
                  m_fluxWeight [ivar][dit()](vof,0) = 0.;
                  m_fluxStencil[ivar][dit()](vof,0) *= 0.0;
                }
            }
        }
    }
  m_isDefined = true;
}

#include "NamespaceFooter.H"
