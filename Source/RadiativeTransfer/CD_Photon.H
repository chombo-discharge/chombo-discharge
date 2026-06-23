/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Photon.H
  @brief  SoA payload for Monte Carlo radiative-transfer photons
  @author Robert Marskar
*/

#ifndef CD_PHOTON_H
#define CD_PHOTON_H

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_GenericParticle.H>
#include <CD_NamespaceHeader.H>

/**
  @brief SoA payload for Monte Carlo radiative-transfer photons.
  @details kappa is the mean absorption coefficient (inverse mean free path); velocity is the photon
  propagation velocity (per component, D_DECL). Weight and position are container-owned columns, NOT
  payload members. All payload columns (kappa and velocity) are checkpointed; position and weight are
  always checkpointed by the container, so no h5Columns subset is declared.
*/
struct Photon
{
  ParticleReal kappa = 0.0; ///< Mean absorption coefficient (inverse mean free path).

  ParticleReal D_DECL(velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0); ///< Velocity (per component).
};

/**
  @brief ParticleTraits specialization listing the payload columns of Photon.
*/
template <>
struct ParticleTraits<Photon>
{
  static constexpr auto columns = std::make_tuple(&Photon::kappa,
                                                  D_DECL(&Photon::velocityX, &Photon::velocityY, &Photon::velocityZ));
};

/**
  @brief Lightweight AoS photon carrying position, weight, velocity, and absorption coefficient.
  @details Bridge type used where photons must own their own position+weight in AoS form, which the SoA
  Photon payload (position/weight are container-owned) cannot: the per-cell ItoKMC reaction/reconcile
  scratch lists. real<0> is the weight, real<1> is the absorption coefficient (kappa), vect<0> is the
  velocity, and the position is the GenericParticle base position. Mirrors the pre-SoA Photon API so the
  reaction code maps over with minimal churn; results are appended into the SoA McPhoto container.
*/
class PhotonSeed : public GenericParticle<2, 1>
{
public:
  /**
    @brief Default constructor. Zero weight/kappa/velocity/position.
  */
  inline PhotonSeed();

  /**
    @brief Full constructor.
    @param[in] a_position Photon position.
    @param[in] a_velocity Photon velocity.
    @param[in] a_kappa    Mean absorption coefficient.
    @param[in] a_weight   Particle weight (a_weight=1 implies a physical photon).
  */
  inline PhotonSeed(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_weight);

  /**
    @brief Destructor.
  */
  inline virtual ~PhotonSeed();

  /**
    @brief Define function. Sets all parameters.
    @param[in] a_position Photon position.
    @param[in] a_velocity Photon velocity.
    @param[in] a_kappa    Mean absorption coefficient.
    @param[in] a_weight   Particle weight (a_weight=1 implies a physical photon).
  */
  inline void
  define(const RealVect& a_position, const RealVect& a_velocity, const Real& a_kappa, const Real a_weight);

  /**
    @brief Get the photon weight.
    @return Reference to the weight.
  */
  inline Real&
  weight();

  /**
    @brief Get the photon weight.
    @return Const reference to the weight.
  */
  inline const Real&
  weight() const;

  /**
    @brief Get the photon absorption coefficient.
    @return Reference to kappa.
  */
  inline Real&
  kappa();

  /**
    @brief Get the photon absorption coefficient.
    @return Const reference to kappa.
  */
  inline const Real&
  kappa() const;

  /**
    @brief Get the photon velocity.
    @return Reference to the velocity.
  */
  inline RealVect&
  velocity();

  /**
    @brief Get the photon velocity.
    @return Const reference to the velocity.
  */
  inline const RealVect&
  velocity() const;
};

#include <CD_NamespaceFooter.H>

#include <CD_PhotonImplem.H>

#endif
