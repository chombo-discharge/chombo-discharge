/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_CdrPlasmaPhysics.H
   @brief  Declaration of the Physics::CdrPlasma::CdrPlasmaPhysics interface class
   @author Robert Marskar
*/

#ifndef CD_CDRPLASMAPHYSICS_H
#define CD_CDRPLASMAPHYSICS_H

// Chombo includes
#include <RealVect.H>
#include <RefCountedPtr.H>
#include <LoHiSide.H>

// Our includes
#include <CD_CdrSpecies.H>
#include <CD_RtSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /**
       @brief Abstract interface for specifying plasma kinetics in the CdrPlasma physics module.
       @details This class defines the coupling functions required by CdrPlasmaStepper to advance
       a system of convection-diffusion-reaction (CDR) equations for charged and neutral species
       together with radiative transfer equations (RTE) for photon densities, all in a
       self-consistent electric field obtained from Poisson's equation.

       Concrete subclasses must implement:
       - advanceReactionNetwork: compute CDR and RTE source terms from local state.
       - computeCdrDriftVelocities: drift velocities for each CDR species.
       - computeCdrDiffusionCoefficients: diffusion coefficients for each CDR species.
       - computeCdrElectrodeFluxes: EB boundary fluxes on electrode surfaces.
       - computeCdrDielectricFluxes: EB boundary fluxes on dielectric surfaces.
       - computeCdrDomainFluxes: fluxes through the computational domain boundary.
       - initialSigma: initial surface charge density on dielectrics.
       - computeAlpha: Townsend ionization coefficient (used by cell taggers).
       - computeEta: Townsend attachment coefficient (used by cell taggers).

       Species are registered by filling m_cdrSpecies and m_rtSpecies in the subclass constructor.
       The ordering of species in all input/output vectors throughout this interface matches the
       ordering in those two containers.
    */
    class CdrPlasmaPhysics
    {
    public:
      /**
         @brief Default constructor. Does nothing.
      */
      CdrPlasmaPhysics()
      {}

      /**
         @brief Destructor. Does nothing.
      */
      virtual ~CdrPlasmaPhysics()
      {}

      /**
         @brief Parse runtime options. Default implementation does nothing.
      */
      virtual void
      parseRuntimeOptions() {};

      /**
         @brief Return the number of plot variables provided by this physics class.
         @details CdrPlasmaStepper calls this to pre-allocate output storage. The returned value must
         match the length of the vectors returned by getPlotVariableNames() and getPlotVariables().
         @return Number of extra plot variables. Default is 0.
      */
      virtual int
      getNumberOfPlotVariables() const
      {
        return 0;
      }

      /**
         @brief Return the names of the extra plot variables.
         @details The length of the returned vector must match getNumberOfPlotVariables(). The
         ordering must be consistent with getPlotVariables().
         @return Vector of variable names.
      */
      virtual Vector<std::string>
      getPlotVariableNames() const
      {
        return Vector<std::string>(this->getNumberOfPlotVariables(), std::string("empty data"));
      }

      /**
         @brief Evaluate the extra plot variables at a single grid point.
         @details Called by CdrPlasmaStepper when writing plot files. The length of the returned
         vector must match getNumberOfPlotVariables().
         @param[in] a_cdrDensities  CDR species densities at this point.
         @param[in] a_cdrGradients  CDR species density gradients at this point.
         @param[in] a_rteDensities  RTE photon densities at this point.
         @param[in] a_E             Electric field vector at this point.
         @param[in] a_position      Physical position.
         @param[in] a_dx            Grid cell spacing.
         @param[in] a_dt            Time step used in the last advance.
         @param[in] a_time          Current simulation time.
         @param[in] a_kappa         Volume fraction of the grid cell.
         @return Vector of plot variable values.
      */
      virtual Vector<Real>
      getPlotVariables(const Vector<Real>&     a_cdrDensities,
                       const Vector<RealVect>& a_cdrGradients,
                       const Vector<Real>&     a_rteDensities,
                       const RealVect&         a_E,
                       const RealVect&         a_position,
                       const Real              a_dx,
                       const Real              a_dt,
                       const Real              a_time,
                       const Real              a_kappa) const
      {
        return Vector<Real>(this->getNumberOfPlotVariables(), 1.0);
      }

      /**
         @brief Compute the Townsend ionization coefficient alpha.
         @details Used primarily by cell tagging classes to identify high-ionization regions.
         @param[in] a_E        Electric field magnitude (V/m).
         @param[in] a_position Physical position.
         @return Townsend ionization coefficient (1/m).
      */
      virtual Real
      computeAlpha(const Real a_E, const RealVect& a_position) const = 0;

      /**
         @brief Compute the Townsend attachment coefficient eta.
         @details Used primarily by cell tagging classes to identify attachment-dominated regions.
         @param[in] a_E        Electric field magnitude (V/m).
         @param[in] a_position Physical position.
         @return Townsend attachment coefficient (1/m).
      */
      virtual Real
      computeEta(const Real a_E, const RealVect& a_position) const = 0;

      /**
         @brief Advance the reaction network over a time interval and return source terms.
         @details The advance is assumed to take the form phi^{k+1} = phi^k + S * dt. Implementations
         may either fill a_cdrSources and a_rteSources directly (explicit rule), or perform a
         fully implicit sub-step and back out the equivalent source S from the result.
         @param[out] a_cdrSources   CDR source terms (same ordering as m_cdrSpecies).
         @param[out] a_rteSources   RTE source terms (same ordering as m_rtSpecies).
         @param[in]  a_cdrDensities CDR species densities at this point.
         @param[in]  a_cdrGradients CDR species density gradients at this point.
         @param[in]  a_rteDensities RTE photon densities at this point.
         @param[in]  a_E            Electric field vector at this point.
         @param[in]  a_pos          Physical position.
         @param[in]  a_dx           Grid cell spacing.
         @param[in]  a_dt           Time step size.
         @param[in]  a_time         Current simulation time.
         @param[in]  a_kappa        Volume fraction of the grid cell.
      */
      virtual void
      advanceReactionNetwork(Vector<Real>&           a_cdrSources,
                             Vector<Real>&           a_rteSources,
                             const Vector<Real>&     a_cdrDensities,
                             const Vector<RealVect>& a_cdrGradients,
                             const Vector<Real>&     a_rteDensities,
                             const RealVect&         a_E,
                             const RealVect&         a_pos,
                             const Real              a_dx,
                             const Real              a_dt,
                             const Real              a_time,
                             const Real              a_kappa) const = 0;

      /**
         @brief Compute drift velocities for all CDR species.
         @param[in] a_time         Current simulation time.
         @param[in] a_pos          Physical position.
         @param[in] a_E            Electric field vector.
         @param[in] a_cdrDensities CDR species densities.
         @return Drift velocity vectors in the same order as m_cdrSpecies.
      */
      virtual Vector<RealVect>
      computeCdrDriftVelocities(const Real          a_time,
                                const RealVect&     a_pos,
                                const RealVect&     a_E,
                                const Vector<Real>& a_cdrDensities) const = 0;

      /**
         @brief Compute diffusion coefficients for all CDR species.
         @param[in] a_time         Current simulation time.
         @param[in] a_pos          Physical position.
         @param[in] a_E            Electric field vector.
         @param[in] a_cdrDensities CDR species densities.
         @return Scalar diffusion coefficients in the same order as m_cdrSpecies.
      */
      virtual Vector<Real>
      computeCdrDiffusionCoefficients(const Real          a_time,
                                      const RealVect&     a_pos,
                                      const RealVect&     a_E,
                                      const Vector<Real>& a_cdrDensities) const = 0;

      /**
         @brief Compute CDR species fluxes on electrode-gas EB interfaces.
         @details Used as a boundary condition on electrode EB cells. A positive flux is directed
         into the gas phase (outward from the electrode).
         @param[in] a_time            Current simulation time.
         @param[in] a_pos             EB centroid position.
         @param[in] a_normal          Outward unit normal pointing into the gas phase.
         @param[in] a_E               Electric field at the EB centroid.
         @param[in] a_cdrDensities    CDR densities at the EB centroid.
         @param[in] a_cdrVelocities   Normal component of CDR drift velocities at the EB.
         @param[in] a_cdrGradients    Normal gradient of CDR densities at the EB.
         @param[in] a_rteFluxes       Normal component of RTE fluxes at the EB.
         @param[in] a_extrapCdrFluxes CDR fluxes extrapolated from the gas-phase interior.
         @return Boundary fluxes in the same order as m_cdrSpecies.
      */
      virtual Vector<Real>
      computeCdrElectrodeFluxes(const Real          a_time,
                                const RealVect&     a_pos,
                                const RealVect&     a_normal,
                                const RealVect&     a_E,
                                const Vector<Real>& a_cdrDensities,
                                const Vector<Real>& a_cdrVelocities,
                                const Vector<Real>& a_cdrGradients,
                                const Vector<Real>& a_rteFluxes,
                                const Vector<Real>& a_extrapCdrFluxes) const = 0;

      /**
         @brief Compute CDR species fluxes on dielectric-gas EB interfaces.
         @details Used as a boundary condition on dielectric EB cells. A positive flux is directed
         into the gas phase (outward from the dielectric).
         @param[in] a_time            Current simulation time.
         @param[in] a_pos             EB centroid position.
         @param[in] a_normal          Outward unit normal pointing into the gas phase.
         @param[in] a_E               Electric field at the EB centroid.
         @param[in] a_cdrDensities    CDR densities at the EB centroid.
         @param[in] a_cdrVelocities   Normal component of CDR drift velocities at the EB.
         @param[in] a_cdrGradients    Normal gradient of CDR densities at the EB.
         @param[in] a_rteFluxes       Normal component of RTE fluxes at the EB.
         @param[in] a_extrapCdrFluxes CDR fluxes extrapolated from the gas-phase interior.
         @return Boundary fluxes in the same order as m_cdrSpecies.
      */
      virtual Vector<Real>
      computeCdrDielectricFluxes(const Real          a_time,
                                 const RealVect&     a_pos,
                                 const RealVect&     a_normal,
                                 const RealVect&     a_E,
                                 const Vector<Real>& a_cdrDensities,
                                 const Vector<Real>& a_cdrVelocities,
                                 const Vector<Real>& a_cdrGradients,
                                 const Vector<Real>& a_rteFluxes,
                                 const Vector<Real>& a_extrapCdrFluxes) const = 0;

      /**
         @brief Compute CDR species fluxes through computational domain boundaries.
         @details A positive flux is directed into the domain.
         @param[in] a_time            Current simulation time.
         @param[in] a_pos             Position on the domain face.
         @param[in] a_dir             Coordinate direction (0=x, 1=y, 2=z).
         @param[in] a_side            Low or high side of the domain.
         @param[in] a_E               Electric field at the domain face.
         @param[in] a_cdrDensities    CDR densities at the domain face.
         @param[in] a_cdrVelocities   Normal component of CDR drift velocities at the domain face.
         @param[in] a_cdrGradients    Normal component of CDR density gradients at the domain face.
         @param[in] a_rteFluxes       Normal component of RTE fluxes at the domain face.
         @param[in] a_extrapCdrFluxes CDR fluxes extrapolated from the interior.
         @return Boundary fluxes in the same order as m_cdrSpecies.
      */
      virtual Vector<Real>
      computeCdrDomainFluxes(const Real           a_time,
                             const RealVect&      a_pos,
                             const int            a_dir,
                             const Side::LoHiSide a_side,
                             const RealVect&      a_E,
                             const Vector<Real>&  a_cdrDensities,
                             const Vector<Real>&  a_cdrVelocities,
                             const Vector<Real>&  a_cdrGradients,
                             const Vector<Real>&  a_rteFluxes,
                             const Vector<Real>&  a_extrapCdrFluxes) const = 0;

      /**
         @brief Return the initial surface charge density on dielectric surfaces.
         @param[in] a_time Current simulation time (typically 0).
         @param[in] a_pos  Physical position on the dielectric surface.
         @return Initial surface charge density (C/m^2 in 3D, C/m in 2D).
      */
      virtual Real
      initialSigma(const Real a_time, const RealVect& a_pos) const = 0;

      /**
         @brief Return the list of CDR species.
         @return Const reference to the CDR species vector.
      */
      const Vector<RefCountedPtr<CdrSpecies>>&
      getCdrSpecies() const
      {
        return m_cdrSpecies;
      }

      /**
         @brief Return the list of RTE species.
         @return Const reference to the RTE species vector.
      */
      const Vector<RefCountedPtr<RtSpecies>>&
      getRtSpecies() const
      {
        return m_rtSpecies;
      }

      /**
         @brief Return the number of CDR species.
         @return Number of CDR species.
      */
      int
      getNumCdrSpecies() const
      {
        return m_cdrSpecies.size();
      }

      /**
         @brief Return the number of RTE equations.
         @return Number of RTE species.
      */
      int
      getNumRtSpecies() const
      {
        return m_rtSpecies.size();
      }

    protected:
      /**
         @brief CDR species definitions. Subclasses populate this in their constructor.
      */
      Vector<RefCountedPtr<CdrSpecies>> m_cdrSpecies;

      /**
         @brief RTE species definitions. Subclasses populate this in their constructor.
      */
      Vector<RefCountedPtr<RtSpecies>> m_rtSpecies;

      /**
         @brief Number of CDR species.
      */
      int m_numCdrSpecies;

      /**
         @brief Number of RTE species.
      */
      int m_numRtSpecies;
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
