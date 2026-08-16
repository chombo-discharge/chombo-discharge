/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ItoMergeParticle.H
 * @brief  Declaration of the reduced Ito particle used for merging and for regrid transport
 * @author Robert Marskar
 */

#ifndef CD_ITOMERGEPARTICLE_H
#define CD_ITOMERGEPARTICLE_H

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Reduced Ito particle: everything an ItoParticle carries that is not recomputed from the mesh.
 * @details This is *the* Ito particle representation for merging, and the holder the bulk particles ride
 * across a regrid in. Position, weight and the particleID/rankID metadata are container-owned columns, so
 * the payload is a single scalar -- the average energy -- and the whole particle is 53 B in 3-D against
 * ItoParticle's 149 B.
 *
 * Nothing else in ItoParticle needs to survive either operation. Its mobility, diffusion and velocity
 * columns are interpolated from the mesh (interpolateMobilities()/interpolateDiffusion()/
 * interpolateVelocities()); its scratch columns are scratch; and old_x/y/z is rewritten wholesale at the
 * top of every time step (ItoKMCGodunovStepper::setOldPositions()) and read back inside the same advance,
 * so no regrid can land between the write and the read. The energy is the one genuine state candidate --
 * it is read on a live path, ItoSolver::interpolateMobilities() feeding ItoSpecies::mobility(Real) -- so
 * it is carried, and dropping the rest is lossless.
 *
 * @warning The old_x/y/z argument above expires the day anything writes those columns outside
 * ItoKMCGodunovStepper::advance(); ItoKMCStepper::computeEdotJSource() reads them and would be such a
 * caller if it were ever wired up (it has no callers anywhere in the tree today).
 */
struct ItoMergeParticle
{
  ParticleReal energy = 0.0; ///< Average particle energy -- the only payload this particle carries.
};

/**
 * @brief ParticleTraits for ItoMergeParticle: a single energy payload column.
 */
template <>
struct ParticleTraits<ItoMergeParticle>
{
  /**
   * @brief The one payload column.
   */
  static constexpr auto columns = std::make_tuple(&ItoMergeParticle::energy);
};

#include <CD_NamespaceFooter.H>

#endif
