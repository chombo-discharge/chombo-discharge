/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzElectrostaticEBBC.cpp
  @brief  Implementation of CD_MFHelmholtzElectrostaticEBBC.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzElectrostaticEBBC.H>
#include <CD_LeastSquares.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzElectrostaticEBBC::MFHelmholtzElectrostaticEBBC(const int                               a_phase,
                                                           const ElectrostaticEbBc&                a_electrostaticBCs,
                                                           const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC)
  : MFHelmholtzEBBC(a_phase, a_jumpBC)
{
  CH_TIME("MFHelmholtzElectrostaticEBBC::MFHelmholtzElectrostaticEBBC");

  m_order  = -1;
  m_weight = -1;

  this->setDomainDropOrder(-1);

  m_electrostaticBCs = a_electrostaticBCs;
}

MFHelmholtzElectrostaticEBBC::~MFHelmholtzElectrostaticEBBC()
{
  CH_TIME("MFHelmholtzElectrostaticEBBC::~MFHelmholtzElectrostaticEBBC");
}

void
MFHelmholtzElectrostaticEBBC::setOrder(const int a_order)
{
  CH_TIME("MFHelmholtzElectrostaticEBBC::setOrder(int)");

  CH_assert(a_order > 0);

  m_order = a_order;
}

void
MFHelmholtzElectrostaticEBBC::setWeight(const int a_weight)
{
  CH_TIME("MFHelmholtzElectrostaticEBBC::setWeight(int)");

  CH_assert(a_weight >= 0);

  m_weight = a_weight;
}

void
MFHelmholtzElectrostaticEBBC::setDomainDropOrder(const int a_domainSize)
{
  CH_TIME("MFHelmholtzElectrostaticEBBC::setDomainDropOrder()");

  m_domainDropOrder = a_domainSize;
}

void
MFHelmholtzElectrostaticEBBC::defineSinglePhase()
{
  CH_TIME("MFHelmholtzElectrostaticEBBC::defineSinglePhase()");

  CH_assert(m_order > 0);
  CH_assert(m_weight >= 0);

  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const ProblemDomain&     domain = m_eblg.getDomain();

  // Drop order if we must
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (domain.size()[dir] <= m_domainDropOrder) {
      m_order = 1;
    }
  }

  const DataIterator& dit = dbl.dataIterator();

  const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box         box     = dbl[din];
    const EBISBox&    ebisbox = m_eblg.getEBISL()[din];
    const IntVectSet& ivs     = ebisbox.getIrregIVS(box);

    BaseIVFAB<Real>&       weights  = m_boundaryWeights[din];
    BaseIVFAB<VoFStencil>& stencils = m_gradPhiStencils[din];

    VoFIterator& singlePhaseVofs = m_jumpBC->getSinglePhaseVofs(m_phase, din);

    auto kernel = [&](const VolIndex& vof) -> void {
      const Real areaFrac = ebisbox.bndryArea(vof);

      int order = -1;

      bool foundStencil = false;
      bool dropOrder    = false;

      std::pair<Real, VoFStencil> pairSten;

      // Drop stencil order if this cell is not a valid grid cell (i.e., one that lies on the AMR grids and is not covered by a finer grid)
      if (!(m_validCells.isNull())) {
        if ((*m_validCells)[din](vof.gridIndex(), 0) == false) {
          dropOrder = false;
        }
      }
      else {
        dropOrder = false;
      }

#warning "Dev code in MFHelmholtzElectrostaticBC in stencil drop order"
#warning "Must remove dev2d.inputs and dev3d.inputs"
      // Try semi-circle first.
      order = dropOrder ? 1 : m_order;
      while (!foundStencil && order > 0) {
        foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten,
                                                                vof,
                                                                VofUtils::Neighborhood::SemiCircle,
                                                                din,
                                                                order,
                                                                m_weight);
        order--;

        // Check if stencil reaches too far across CF
        if (foundStencil) {
          foundStencil = this->isStencilValidCF(pairSten.second, din);
        }

        // If the stencil consists of only irregular cells, drop order.
            // If the stencil consists of only irregular cells, drop order.
	    bool allIrregular = true;

            const VoFStencil& sten = pairSten.second;
            for (int i = 0; i < sten.size(); i++) {
              const VolIndex& vof = sten.vof(i);

              if (ebisbox.isRegular(vof.gridIndex())) {
		allIrregular = false;
              }
            }

	    //            const bool tooIrregular = numIrregular > 0;

            if (foundStencil && allIrregular && order > 0) {
              //	  std::cout << "dropping that motherfucking order" << std::endl;
              foundStencil = false;
            }	
        // int               numIrregular = 0;
	
        // const VoFStencil& sten         = pairSten.second;
        // for (int i = 0; i < sten.size(); i++) {
        //   const VolIndex& vof = sten.vof(i);

        //   if (ebisbox.isIrregular(vof.gridIndex())) {
	//     numIrregular++;
        //   }
        // }

        // const bool tooIrregular = numIrregular > 0;

        // if (foundStencil && tooIrregular && order> 0) {
	//   //	  std::cout << "dropping that motherfucking order" << std::endl;
        //   foundStencil = false;
        // }
      }

      // Try quadrant if that didn't work.
      order = dropOrder ? 1 : m_order;
      while (!foundStencil && order > 0) {
        foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten,
                                                                vof,
                                                                VofUtils::Neighborhood::Quadrant,
                                                                din,
                                                                order,
                                                                m_weight);
        order--;

        // Check if stencil reaches too far across CF
        if (foundStencil) {
          foundStencil = this->isStencilValidCF(pairSten.second, din);
        }
      }

      // Last ditch effort: Try a full radius
      order = dropOrder ? 1 : m_order;
      while (!foundStencil && order > 0) {
        foundStencil = this->getLeastSquaresBoundaryGradStencil(pairSten,
                                                                vof,
                                                                VofUtils::Neighborhood::Radius,
                                                                din,
                                                                order,
                                                                m_weight);
        order--;

        // Check if stencil reaches too far across CF
        if (foundStencil) {
          foundStencil = this->isStencilValidCF(pairSten.second, din);
        }
      }

      if (foundStencil) {
        weights(vof, m_comp)  = pairSten.first;
        stencils(vof, m_comp) = pairSten.second;

        // Stencil and weight must also be scaled by the B-coefficient, dx (because it's used in kappa*Div(F)) and the area fraction.
        weights(vof, m_comp) *= areaFrac / m_dx;
        stencils(vof, m_comp) *= areaFrac / m_dx;
      }
      else {
        // Dead cell. No flux
        const std::string baseErr = "MFHelmholtzElectrostaticEBBC::defineSinglePhase - dead cell on domain = ";
        const std::string vofErr  = " on vof = ";
        const std::string impErr  = " (this may cause multigrid divergence)";

        //        std::cout << baseErr << m_eblg.getDomain() << vofErr << vof << impErr << std::endl;

        weights(vof, m_comp) = 0.0;
        stencils(vof, m_comp).clear();
      }
    };

    BoxLoops::loop(singlePhaseVofs, kernel);
  }
}

void
MFHelmholtzElectrostaticEBBC::applyEBFluxSinglePhase(VoFIterator&           a_singlePhaseVofs,
                                                     EBCellFAB&             a_Lphi,
                                                     const EBCellFAB&       a_phi,
                                                     const BaseIVFAB<Real>& a_Bcoef,
                                                     const DataIndex&       a_dit,
                                                     const Real&            a_beta,
                                                     const bool&            a_homogeneousPhysBC) const
{
  CH_TIME(
    "MFHelmholtzElectrostaticEBBC::applyEBFluxSinglePhase(VoFIterator, EBCellFAB, EBCellFAB, DataIndex, Real, bool)");

  // Apply the stencil for computing the contribution to kappaDivF. Note divF is sum(faces) B*grad(Phi)/dx and that this
  // is the contribution from the EB face. B/dx is already included in the stencils and boundary weights, but beta is not.

  // This is a safeguard against a corner case where we fail to correctly represent the interface region and also have
  // no electrodes. In this case we simply bypass the flux calculation. Note that this normally happens when the user
  // only has dielectrics AND some of the dielectric cells become incorrectly represented at the coarser multigrid levels.
  const bool hasElectrode = m_electrostaticBCs.getBcs().size() > 0;

  // Do single phase cells
  if (!a_homogeneousPhysBC && hasElectrode) {
    auto kernel = [&](const VolIndex& vof) -> void {
      const RealVect pos   = this->getBoundaryPosition(vof, a_dit);
      const Real     value = this->getElectrodePotential(pos);
      const Real     Bcoef = a_Bcoef(vof, m_comp);

      a_Lphi(vof, m_comp) += a_beta * Bcoef * value * m_boundaryWeights[a_dit](vof, m_comp);
    };

    BoxLoops::loop(a_singlePhaseVofs, kernel);
  }

  return;
}

#include <CD_NamespaceFooter.H>
