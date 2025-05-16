/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBParticleMeshImplem.H
  @brief  Implementatoin of CD_EBParticleMesh.H
  @author Robert Marskar
*/

#ifndef CD_EBParticleMeshImplem_H
#define CD_EBParticleMeshImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBParticleMesh.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

template <class P, const Real& (P::*particleScalarField)() const>
void
EBParticleMesh::deposit(const List<P>&       a_particleList,
                        EBCellFAB&           a_rho,
                        const DepositionType a_depositionType,
                        const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit");

  const Interval variables(0, 0);

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();

    this
      ->depositParticle(a_rho, m_probLo, m_dx, curPosition, &curStrength, variables, a_depositionType, a_forceIrregNGP);
  }
}

template <class P, Real (P::*particleScalarField)() const>
void
EBParticleMesh::deposit(const List<P>&       a_particleList,
                        EBCellFAB&           a_rho,
                        const DepositionType a_depositionType,
                        const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit");

  const Interval variables(0, 0);

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();

    this
      ->depositParticle(a_rho, m_probLo, m_dx, curPosition, &curStrength, variables, a_depositionType, a_forceIrregNGP);
  }
}

template <class P, const Real& (P::*particleScalarField)() const>
void
EBParticleMesh::deposit2(const List<P>&       a_particleList,
                         EBCellFAB&           a_rho,
                         const DepositionType a_depositionType,
                         const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit2");

  const Interval variables(0, 0);

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();

    this->depositParticle2(a_rho,
                           m_probLo,
                           m_dx,
                           curPosition,
                           &curStrength,
                           variables,
                           a_depositionType,
                           a_forceIrregNGP);
  }
}

template <class P, Real (P::*particleScalarField)() const>
void
EBParticleMesh::deposit2(const List<P>&       a_particleList,
                         EBCellFAB&           a_rho,
                         const DepositionType a_depositionType,
                         const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit2");

  const Interval variables(0, 0);

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();

    this->depositParticle2(a_rho,
                           m_probLo,
                           m_dx,
                           curPosition,
                           &curStrength,
                           variables,
                           a_depositionType,
                           a_forceIrregNGP);
  }
}

template <class P, const Real& (P::*particleScalarField)() const>
void
EBParticleMesh::deposit4(const List<P>&       a_particleList,
                         EBCellFAB&           a_rho,
                         const DepositionType a_depositionType,
                         const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit4");

  const Interval variables(0, 0);

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();

    this->depositParticle4(a_rho,
                           m_probLo,
                           m_dx,
                           curPosition,
                           &curStrength,
                           variables,
                           a_depositionType,
                           a_forceIrregNGP);
  }
}

template <class P, Real (P::*particleScalarField)() const>
void
EBParticleMesh::deposit4(const List<P>&       a_particleList,
                         EBCellFAB&           a_rho,
                         const DepositionType a_depositionType,
                         const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit4");

  const Interval variables(0, 0);

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();

    this->depositParticle4(a_rho,
                           m_probLo,
                           m_dx,
                           curPosition,
                           &curStrength,
                           variables,
                           a_depositionType,
                           a_forceIrregNGP);
  }
}

template <class P, const RealVect& (P::*particleVectorField)() const>
void
EBParticleMesh::deposit(const List<P>&       a_particleList,
                        EBCellFAB&           a_rho,
                        const DepositionType a_depositionType,
                        const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit");

  // TLDR: This is a jack-of-all-trades deposition function. The user will use this function to supply a pointer to the field that will be
  //       deposited. As per API, this function must be of the type 'const Real& myParticleClass::myDepositionField() const'

  const Interval variables(0, SpaceDim - 1);

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const RealVect& curStrength = (curParticle.*particleVectorField)();

    this->depositParticle(a_rho,
                          m_probLo,
                          m_dx,
                          curPosition,
                          curStrength.dataPtr(),
                          variables,
                          a_depositionType,
                          a_forceIrregNGP);
  }
}

template <class P, RealVect (P::*particleVectorField)() const>
void
EBParticleMesh::deposit(const List<P>&       a_particleList,
                        EBCellFAB&           a_rho,
                        const DepositionType a_depositionType,
                        const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit");

  const Interval variables(0, SpaceDim - 1);

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const RealVect& curStrength = (curParticle.*particleVectorField)();

    this->depositParticle(a_rho,
                          m_probLo,
                          m_dx,
                          curPosition,
                          curStrength.dataPtr(),
                          variables,
                          a_depositionType,
                          a_forceIrregNGP);
  }
}

template <class P, const RealVect& (P::*particleVectorField)() const>
void
EBParticleMesh::deposit2(const List<P>&       a_particleList,
                         EBCellFAB&           a_rho,
                         const DepositionType a_depositionType,
                         const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit2");

  const Interval variables(0, SpaceDim - 1);

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const RealVect& curStrength = (curParticle.*particleVectorField)();

    this->depositParticle2(a_rho,
                           m_probLo,
                           m_dx,
                           curPosition,
                           curStrength.dataPtr(),
                           variables,
                           a_depositionType,
                           a_forceIrregNGP);
  }
}

template <class P, RealVect (P::*particleVectorField)() const>
void
EBParticleMesh::deposit2(const List<P>&       a_particleList,
                         EBCellFAB&           a_rho,
                         const DepositionType a_depositionType,
                         const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit2");

  const Interval variables(0, SpaceDim - 1);

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const RealVect& curStrength = (curParticle.*particleVectorField)();

    this->depositParticle2(a_rho,
                           m_probLo,
                           m_dx,
                           curPosition,
                           curStrength.dataPtr(),
                           variables,
                           a_depositionType,
                           a_forceIrregNGP);
  }
}

template <class P, const RealVect& (P::*particleVectorField)() const>
void
EBParticleMesh::deposit4(const List<P>&       a_particleList,
                         EBCellFAB&           a_rho,
                         const DepositionType a_depositionType,
                         const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit4");

  const Interval variables(0, SpaceDim - 1);

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const RealVect& curStrength = (curParticle.*particleVectorField)();

    this->depositParticle4(a_rho,
                           m_probLo,
                           m_dx,
                           curPosition,
                           curStrength.dataPtr(),
                           variables,
                           a_depositionType,
                           a_forceIrregNGP);
  }
}

template <class P, RealVect (P::*particleVectorField)() const>
void
EBParticleMesh::deposit4(const List<P>&       a_particleList,
                         EBCellFAB&           a_rho,
                         const DepositionType a_depositionType,
                         const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::deposit4");

  const Interval variables(0, SpaceDim - 1);

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const RealVect& curStrength = (curParticle.*particleVectorField)();

    this->depositParticle4(a_rho,
                           m_probLo,
                           m_dx,
                           curPosition,
                           curStrength.dataPtr(),
                           variables,
                           a_depositionType,
                           a_forceIrregNGP);
  }
}

template <class P, Real& (P::*particleScalarField)()>
void
EBParticleMesh::interpolate(List<P>&             a_particleList,
                            const EBCellFAB&     a_meshScalarField,
                            const DepositionType a_interpType,
                            const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::interpolate(Real)");

  CH_assert(a_meshScalarField.nComp() == 1);

  const Interval variables(0, 0);

  Box validBox = m_domain.domainBox();

  switch (a_interpType) {
  case DepositionType::NGP: {
    validBox = m_domain.domainBox();

    break;
  }
  case DepositionType::CIC: {
    validBox = grow(validBox, -1);

    break;
  }
  case DepositionType::TSC: {
    validBox = grow(validBox, -2);

    break;
  }
  case DepositionType::W4: {
    validBox = grow(validBox, -3);

    break;
  }
  default: {
    MayDay::Error("EBParticleMesh::interpolate - logic bust");
  }
  }

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    P&              curParticle = lit();
    const RealVect& curPosition = curParticle.position();

    Real& curParticleField = (curParticle.*particleScalarField)();

    this->interpolateParticle(&curParticleField,
                              a_meshScalarField,
                              validBox,
                              m_probLo,
                              m_dx,
                              curPosition,
                              variables,
                              a_interpType,
                              a_forceIrregNGP);
  }
}

template <class P, RealVect& (P::*particleVectorField)()>
void
EBParticleMesh::interpolate(List<P>&             a_particleList,
                            const EBCellFAB&     a_meshVectorField,
                            const DepositionType a_interpType,
                            const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::interpolate(RealVect)");

  CH_assert(a_meshVectorField.nComp() == SpaceDim);

  const Interval variables(0, SpaceDim - 1);

  // TLDR: This is a jack-of-all-trades interpolation function. The user will use this function to supply a pointer to the field that will be
  //       interpolated to. As per API, this function must be of the type 'RealVect& myParticleClass::myVectorVariable()'

  Box validBox = m_domain.domainBox();

  switch (a_interpType) {
  case DepositionType::NGP: {
    validBox = m_domain.domainBox();

    break;
  }
  case DepositionType::CIC: {
    validBox = grow(validBox, -1);

    break;
  }
  case DepositionType::TSC: {
    validBox = grow(validBox, -2);

    break;
  }
  case DepositionType::W4: {
    validBox = grow(validBox, -3);

    break;
  }
  default: {
    MayDay::Error("EBParticleMesh::interpolate - logic bust");
  }
  }

  for (ListIterator<P> lit(a_particleList); lit; ++lit) {
    P&              curParticle = lit();
    const RealVect& curPosition = curParticle.position();

    RealVect& curParticleField = (curParticle.*particleVectorField)();

    this->interpolateParticle(curParticleField.dataPtr(),
                              a_meshVectorField,
                              validBox,
                              m_probLo,
                              m_dx,
                              curPosition,
                              variables,
                              a_interpType,
                              a_forceIrregNGP);
  }
}

inline void
EBParticleMesh::depositParticle(EBCellFAB&           a_rho,
                                const RealVect&      a_probLo,
                                const RealVect&      a_dx,
                                const RealVect&      a_position,
                                const Real*          a_strength,
                                const Interval       a_components,
                                const DepositionType a_depositionType,
                                const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::depositParticle");

  const int startComp = a_components.begin();
  const int endComp   = a_components.end();

  CH_assert(a_rho.nComp() >= endComp - 1);

  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.

  // Nifty lambda for converting RealVect position to IntVect (lower-left corner). Note that the input vector
  // is the displacement from a_probLo and not the physical position.
  const auto& dx           = a_dx;
  auto        getCellIndex = [&dx](const RealVect a_rv) -> IntVect {
    return IntVect(D_DECL(std::floor(a_rv[0] / dx[0]), std::floor(a_rv[1] / dx[1]), std::floor(a_rv[2] / dx[2])));
  };

  // Get grid cell corresponding to this particle.
  const IntVect particleIndex = getCellIndex(a_position - a_probLo);

  // Assertion -- particle must live on this patch.
  CH_assert(m_region.contains(particleIndex));

  // Factors needed for the kernels. The little lambda will give us the distance
  // between the cell we're looking at and the particle position (in units of the grid resolution)
  const Real invVol = 1. / (std::pow(a_dx[0], SpaceDim) * 1E-4);

  auto particleDisplacement = [&](const IntVect& iv) -> RealVect {
    return (a_probLo + (RealVect(iv) + 0.5 * RealVect::Unit) * a_dx - a_position) / a_dx;
  };

  // Get regular data.
  FArrayBox& rho = a_rho.getFArrayBox();

  // We can force NGP deposition in cut-cells if we want.
  if (m_ebisbox.isIrregular(particleIndex) && a_forceIrregNGP) {
    for (int comp = startComp; comp <= endComp; comp++) {
      rho(particleIndex, comp) += a_strength[comp] * invVol;
    }
  }
  else {
    switch (a_depositionType) {
    case DepositionType::NGP: {
      for (int comp = startComp; comp <= endComp; comp++) {
        rho(particleIndex, comp) += a_strength[comp] * invVol;
      }
      break;
    }
    case DepositionType::CIC: {
      // Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
      // the box corresponding to the cloud.
      const IntVect loIndex = getCellIndex(a_position - a_probLo - 0.5 * a_dx);
      const Box     cicBox  = Box(loIndex, loIndex + IntVect::Unit);

      // This is the cloud-in-cell deposition kernel.
      auto cicKernel = [&](const IntVect& iv) -> void {
        Real weight = invVol;

        // L is the distance between the grid cell and the particle.
        const RealVect L = particleDisplacement(iv);
        for (int dir = 0; dir < SpaceDim; dir++) {
          weight *= (1. - std::abs(L[dir]));
        }

        for (int comp = startComp; comp <= endComp; comp++) {
          rho(iv, comp) += weight * a_strength[comp];
        }
      };

      // Add mass to cells.
      BoxLoops::loop(cicBox, cicKernel);

      break;
    }
    case DepositionType::TSC: {

      // Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
      // the box corresponding to the cloud.
      const IntVect loIndex = getCellIndex(a_position - a_probLo - a_dx);
      const Box     tscBox  = Box(loIndex, loIndex + 2 * IntVect::Unit);

      // This is the triangle-shaped cloud kernel.
      auto tscKernel = [&](const IntVect& iv) -> void {
        Real weight = invVol;

        const RealVect L = particleDisplacement(iv);
        for (int dir = 0; dir < SpaceDim; dir++) {
          const Real l = std::abs(L[dir]);

          if (l < 0.5) {
            weight *= 0.75 - l * l;
          }
          else {
            weight *= 0.5 * (1.5 - l) * (1.5 - l);
          }
        }

        for (int comp = startComp; comp <= endComp; comp++) {
          rho(iv, comp) += weight * a_strength[comp];
        }
      };

      // Run the kernel.
      BoxLoops::loop(tscBox, tscKernel);

      break;
    }
    case DepositionType::W4: {

      // Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
      // the box corresponding to the cloud.
      const IntVect loIndex = getCellIndex(a_position - a_probLo - 1.5 * a_dx);
      const Box     w4Box   = Box(loIndex, loIndex + 3 * IntVect::Unit);

      // This is the fourth order kernel.
      auto w4Kernel = [&](const IntVect& iv) -> void {
        Real weight = invVol;

        // L is the distance between the grid cell and the particle.
        const RealVect L = particleDisplacement(iv);
        for (int dir = 0; dir < SpaceDim; dir++) {
          const Real l = std::abs(L[dir]);

          if (l < 1.0) {
            weight *= 1.0 - 2.5 * l * l + 1.5 * l * l * l;
          }
          else {
            weight *= 0.5 * (2. - l) * (2. - l) * (1. - l);
          }
        }

        for (int comp = startComp; comp <= endComp; comp++) {
          rho(iv, comp) += weight * a_strength[comp];
        }
      };

      // Run the kernel.
      BoxLoops::loop(w4Box, w4Kernel);

      break;
    }
    default: {
      MayDay::Error("EBParticleMesh::depositParticle - logic bust, unknown particle deposition.");
      break;
    }
    }
  }
}

inline void
EBParticleMesh::depositParticle2(EBCellFAB&           a_rho,
                                 const RealVect&      a_probLo,
                                 const RealVect&      a_dx,
                                 const RealVect&      a_position,
                                 const Real*          a_strength,
                                 const Interval       a_components,
                                 const DepositionType a_depositionType,
                                 const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::depositParticle2");

  const int startComp = a_components.begin();
  const int endComp   = a_components.end();

  CH_assert(a_rho.nComp() >= endComp - 1);

  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.
  //
  //       Note that this is the version which deposits with particle widths that are 4 times the "usual" width. Currently, only NGP
  //       and CIC is supported.

  // Nifty lambda for converting RealVect position to IntVect (lower-left corner). Note that the input vector
  // is the displacement from a_probLo and not the physical position.
  const auto& dx           = a_dx;
  auto        getCellIndex = [&dx](const RealVect a_rv) -> IntVect {
    return IntVect(D_DECL(std::floor(a_rv[0] / dx[0]), std::floor(a_rv[1] / dx[1]), std::floor(a_rv[2] / dx[2])));
  };

  // Get grid cell corresponding to this particle.
  const IntVect particleIndex = getCellIndex(a_position - a_probLo);

  // Assertion -- particle must live on this patch.
  CH_assert(m_region.contains(particleIndex));

  // Factors needed for the kernels. The little lambda will give us the distance
  // between the cell we're looking at and the particle position (in units of the grid resolution)
  const Real invVol = 1. / (std::pow(a_dx[0], SpaceDim) * 1.E-4);

  auto particleDisplacement = [&](const IntVect& iv) -> RealVect {
    return (a_probLo + (RealVect(iv) + 0.5 * RealVect::Unit) * a_dx - a_position) / a_dx;
  };

  // Get regular data.
  FArrayBox& rho = a_rho.getFArrayBox();

  // Force NGP in cut-cells if we want.
  if (m_ebisbox.isIrregular(particleIndex) && a_forceIrregNGP) {
    for (int comp = startComp; comp <= endComp; comp++) {
      rho(particleIndex, comp) += a_strength[comp] * invVol;
    }
  }
  else {
    switch (a_depositionType) {
    case DepositionType::NGP: {
      for (int comp = startComp; comp <= endComp; comp++) {
        rho(particleIndex, comp) += a_strength[comp] * invVol;
      }

      break;
    }
    case DepositionType::CIC: {
      // Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
      // the box corresponding to the cloud.
      const IntVect loIndex = getCellIndex(a_position - a_probLo - 1.0 * a_dx);
      const Box     cic2Box = Box(loIndex, loIndex + 2 * IntVect::Unit);

      // This is the kernel for a CIC particle with with 2*dx
      auto cic2Kernel = [&](const IntVect& iv) -> void {
        Real weight = invVol;

        // L is the distance between the grid cell and the particle.
        const RealVect L = particleDisplacement(iv);
        for (int dir = 0; dir < SpaceDim; dir++) {
          const Real l = std::abs(L[dir]);

          if (l > 0.5) {
            weight *= 0.5 * (1.5 - l);
          }
          else {
            weight *= 0.5;
          }
        }

        for (int comp = startComp; comp <= endComp; comp++) {
          rho(iv, comp) += a_strength[comp] * weight;
        }
      };

      BoxLoops::loop(cic2Box, cic2Kernel);

      break;
    }
    default: {
      MayDay::Error(
        "EBParticleMesh::depositParticle2 - Invalid deposition type - only NGP and CIC supported for this deposition method. TSC/W4 have not been worked out.");
      break;
    }
    }
  }
}

inline void
EBParticleMesh::depositParticle4(EBCellFAB&           a_rho,
                                 const RealVect&      a_probLo,
                                 const RealVect&      a_dx,
                                 const RealVect&      a_position,
                                 const Real*          a_strength,
                                 const Interval       a_components,
                                 const DepositionType a_depositionType,
                                 const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::depositParticle2");

  const int startComp = a_components.begin();
  const int endComp   = a_components.end();

  CH_assert(a_rho.nComp() >= endComp - 1);

  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.
  //
  //       Note that this is the version which deposits with particle widths that are 4 times the "usual" width. Currently, only NGP
  //       and CIC is supported.

  // Nifty lambda for converting RealVect position to IntVect (lower-left corner). Note that the input vector
  // is the displacement from a_probLo and not the physical position.
  const auto& dx           = a_dx;
  auto        getCellIndex = [&dx](const RealVect a_rv) -> IntVect {
    return IntVect(D_DECL(std::floor(a_rv[0] / dx[0]), std::floor(a_rv[1] / dx[1]), std::floor(a_rv[2] / dx[2])));
  };

  // Get grid cell corresponding to this particle.
  const IntVect particleIndex = getCellIndex(a_position - a_probLo);

  // Assertion -- particle must live on this patch.
  CH_assert(m_region.contains(particleIndex));

  // Factors needed for the kernels. The little lambda will give us the distance
  // between the cell we're looking at and the particle position (in units of the grid resolution)
  const Real invVol = 1.0 / (std::pow(a_dx[0], SpaceDim)*1.E-4);

  auto particleDisplacement = [&](const IntVect& iv) -> RealVect {
    return (a_probLo + (RealVect(iv) + 0.5 * RealVect::Unit) * a_dx - a_position) / a_dx;
  };

  // Get regular data.
  FArrayBox& rho = a_rho.getFArrayBox();

  // Force NGP in cut-cells if we want.
  if (m_ebisbox.isIrregular(particleIndex) && a_forceIrregNGP) {
    for (int comp = startComp; comp <= endComp; comp++) {
      rho(particleIndex, comp) += a_strength[comp] * invVol;
    }
  }
  else {
    switch (a_depositionType) {
    case DepositionType::NGP: {
      for (int comp = startComp; comp <= endComp; comp++) {
        rho(particleIndex, comp) += a_strength[comp] * invVol;
      }
      break;
    }
    case DepositionType::CIC: {
      // Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
      // the box corresponding to the cloud.
      const IntVect loIndex = getCellIndex(a_position - a_probLo - 2.0 * a_dx);
      const Box     cic4Box = Box(loIndex, loIndex + 4 * IntVect::Unit);

      // This is the kernel for a CIC particle with with 4*dx
      auto cic4Kernel = [&](const IntVect& iv) -> void {
        Real weight = invVol;

        // L is the distance between the grid cell and the particle -- in units of dx.
        const RealVect L = particleDisplacement(iv);
        for (int dir = 0; dir < SpaceDim; dir++) {
          const Real l = std::abs(L[dir]);

          if (l >= 1.5) {
            weight *= 0.25 * (2.5 - l);
          }
          else {
            weight *= 0.25;
          }
        }

        for (int comp = startComp; comp <= endComp; comp++) {
          rho(iv, comp) += a_strength[comp] * weight;
        }
      };

      BoxLoops::loop(cic4Box, cic4Kernel);

      break;
    }
    default: {
      MayDay::Error(
        "EBParticleMesh::depositParticle4 - Invalid deposition type - only NGP and CIC supported for this deposition method. TSC/W4 have not been worked out.");
      break;
    }
    }
  }
}

inline void
EBParticleMesh::interpolateParticle(Real*                a_particleField,
                                    const EBCellFAB&     a_meshField,
                                    const Box&           a_validBox,
                                    const RealVect&      a_probLo,
                                    const RealVect&      a_dx,
                                    const RealVect&      a_position,
                                    const Interval&      a_interval,
                                    const DepositionType a_interpType,
                                    const bool           a_forceIrregNGP) const
{
  CH_TIME("EBParticleMesh::interpolateParticle");

  const int startComp = a_interval.begin();
  const int endComp   = a_interval.end();

  CH_assert(a_meshField.nComp() >= endComp - 1);

  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.

  // Nifty lambda for converting RealVect position to IntVect (lower-left corner). Note that the input vector
  // is the displacement from a_probLo and not the physical position.
  const auto& dx = a_dx;

  auto getCellIndex = [&dx](const RealVect a_rv) -> IntVect {
    return IntVect(D_DECL(std::floor(a_rv[0] / dx[0]), std::floor(a_rv[1] / dx[1]), std::floor(a_rv[2] / dx[2])));
  };

  // Get grid cell corresponding to this particle.
  const IntVect particleIndex = getCellIndex(a_position - a_probLo);

  // Assertion -- particle must live on this patch.
  CH_assert(m_region.contains(particleIndex));

  // Nifty little lambda for computing the displacement between the particle and a grid cell.
  auto particleDisplacement = [&](const IntVect& iv) -> RealVect {
    return (a_probLo + (RealVect(iv) + 0.5 * RealVect::Unit) * a_dx - a_position) / a_dx;
  };

  // Get regular data.
  const FArrayBox& meshField = a_meshField.getFArrayBox();

  // Irregular cells always do an NGP deposit to prevent clouds leaking into the other side.
  if ((m_ebisbox.isIrregular(particleIndex) && a_forceIrregNGP) || !(a_validBox.contains(particleIndex))) {
    for (int comp = startComp; comp <= endComp; comp++) {
      a_particleField[comp] = meshField(particleIndex, comp);
    }
  }
  else if (m_ebisbox.isCovered(particleIndex)) { // Need to set to something.
    for (int comp = startComp; comp <= endComp; comp++) {
      a_particleField[comp] = 0.0;
    }
  }
  else {
    switch (a_interpType) {
    case DepositionType::NGP: // Hook for nearest-grid-point interpolation.
    {
      for (int comp = startComp; comp <= endComp; comp++) {
        a_particleField[comp] = meshField(particleIndex, comp);
      }

      break;
    }
    case DepositionType::CIC: // Hook for cloud-in-cell interpolation.
    {
      // Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
      // the box corresponding to the cloud. This region contains the cells where we gather the field.
      const IntVect loIndex = getCellIndex(a_position - a_probLo - 0.5 * a_dx);
      const Box     cicBox  = Box(loIndex, loIndex + IntVect::Unit);

      // This is the cloud-in-cell interpolation kernel. It computes the weight from the current cell and
      // gathers the force.
      auto cicKernel = [&](const IntVect& iv) -> void {
        Real weight = 1.0;

        // L is the distance between the grid cell and the particle.
        const RealVect L = particleDisplacement(iv);
        for (int dir = 0; dir < SpaceDim; dir++) {
          weight *= (1. - std::abs(L[dir]));
        }

        for (int comp = startComp; comp <= endComp; comp++) {
          a_particleField[comp] += weight * meshField(iv, comp);
        }
      };

      // Run kernel and gather field from mesh.
      for (int comp = startComp; comp <= endComp; comp++) {
        a_particleField[comp] = 0.0;
      }

      BoxLoops::loop(cicBox, cicKernel);

      break;
    }
    case DepositionType::TSC: // Hook for triangle-shaped cloud interpolation.
    {
      // Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
      // the box corresponding to the cloud. This region contains the cells where we gather the field.
      const IntVect loIndex = getCellIndex(a_position - a_probLo - 1.0 * a_dx);
      const Box     tscBox  = Box(loIndex, loIndex + 2 * IntVect::Unit);

      // This is the triangle-shaped-cloud interpolation kernel. It computes the weight from the current cell and
      // gather the force.
      auto tscKernel = [&](const IntVect& iv) -> void {
        Real weight = 1.0;

        // L is the distance between the grid cell and the particle.
        const RealVect L = particleDisplacement(iv);
        for (int dir = 0; dir < SpaceDim; dir++) {
          const Real& l = std::abs(L[dir]);

          if (l < 0.5) {
            weight *= 0.75 - l * l;
          }
          else {
            weight *= 0.5 * (1.5 - l) * (1.5 - l);
          }
        }

        for (int comp = startComp; comp <= endComp; comp++) {
          a_particleField[comp] += weight * meshField(iv, comp);
        }
      };

      // Run kernel and gather field from mesh.
      for (int comp = startComp; comp <= endComp; comp++) {
        a_particleField[comp] = 0.0;
      }

      BoxLoops::loop(tscBox, tscKernel);

      break;
    }
    case DepositionType::W4: // Hook for fourth order interpolation.
    {
      // Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
      // the box corresponding to the cloud. This region contains the cells where we gather the field.
      const IntVect loIndex = getCellIndex(a_position - a_probLo - 1.5 * a_dx);
      const Box     w4Box   = Box(loIndex, loIndex + 3 * IntVect::Unit);

      // This is the fourth order kernel.
      auto w4Kernel = [&](const IntVect& iv) -> void {
        Real weight = 1.0;

        // L is the distance between the grid cell and the particle.
        const RealVect L = particleDisplacement(iv);
        for (int dir = 0; dir < SpaceDim; dir++) {
          const Real l = std::abs(L[dir]);

          if (l < 1.0) {
            weight *= 1.0 - 2.5 * l * l + 1.5 * l * l * l;
          }
          else {
            weight *= 0.5 * (2. - l) * (2. - l) * (1. - l);
          }
        }

        for (int comp = startComp; comp <= endComp; comp++) {
          a_particleField[comp] += weight * meshField(iv, comp);
        }
      };

      // Run kernel and gather field from mesh.
      for (int comp = startComp; comp <= endComp; comp++) {
        a_particleField[comp] = 0.0;
      }

      BoxLoops::loop(w4Box, w4Kernel);

      break;
    }
    default:
      MayDay::Error("EBParticleMesh::interpolateParticle(RealVect) - Invalid interpolation type requested.");
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
