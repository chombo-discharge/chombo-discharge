/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaMorrowBourdon.H
  @brief  Declaration of the Morrow-Lowke plasma model for air
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaMorrowBourdon_H
#define CD_CdrPlasmaMorrowBourdon_H

// Our includes
#include <CD_CdrPlasmaPhysics.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {
  
    /*!
      @brief Class that implements the Morrow-Lowke discharge model for air.
      @details This class uses the Morrow-Lowke kinetic scheme for air and the three-group Bourdon model for the RTE equations. 
    */
    class CdrPlasmaMorrowBourdon : public CdrPlasmaPhysics {
    public:

      CdrPlasmaMorrowBourdon();
      ~CdrPlasmaMorrowBourdon();

      void advanceReactionNetwork(Vector<Real>&          a_particle_sources,
				  Vector<Real>&          a_Photon_sources,
				  const Vector<Real>     a_particle_densities,
				  const Vector<RealVect> a_particle_gradients,
				  const Vector<Real>     a_Photon_densities,
				  const RealVect         a_E,
				  const RealVect         a_pos,
				  const Real             a_dx,
				  const Real             a_dt,
				  const Real             a_time,
				  const Real             a_kappa) const;


      Vector<Real> computeCdrDiffusionCoefficients(const Real         a_time,
						   const RealVect     a_pos,
						   const RealVect     a_E,
						   const Vector<Real> a_cdr_densities) const;


      Vector<RealVect> computeCdrDriftVelocities(const Real         a_time,
						 const RealVect     a_pos,
						 const RealVect     a_E,
						 const Vector<Real> a_cdr_densities) const;


      Vector<Real> computeCdrElectrodeFluxes(const Real         a_time,
					     const RealVect     a_pos,
					     const RealVect     a_normal,
					     const RealVect     a_E,
					     const Vector<Real> a_cdr_densities,
					     const Vector<Real> a_cdr_velocities,
					     const Vector<Real> a_cdr_gradients,
					     const Vector<Real> a_rte_fluxes,
					     const Vector<Real> a_extrap_cdr_fluxes) const;

      Vector<Real> computeCdrDielectricFluxes(const Real         a_time,
					      const RealVect     a_pos,
					      const RealVect     a_normal,
					      const RealVect     a_E,
					      const Vector<Real> a_cdr_densities,
					      const Vector<Real> a_cdr_velocities,
					      const Vector<Real> a_cdr_gradients,
					      const Vector<Real> a_rte_fluxes,
					      const Vector<Real> a_extrap_cdr_fluxes) const;

      Vector<Real> computeCdrDomainFluxes(const Real           a_time,
					  const RealVect       a_pos,
					  const int            a_dir,
					  const Side::LoHiSide a_side,
					  const RealVect       a_E,
					  const Vector<Real>   a_cdr_densities,
					  const Vector<Real>   a_cdr_velocities,
					  const Vector<Real>   a_cdr_gradients,
					  const Vector<Real>   a_rte_fluxes,
					  const Vector<Real>   a_extrap_cdr_fluxes) const;

      Real initialSigma(const Real a_time, const RealVect a_pos) const;



      class Electron : public CdrSpecies {
      public:

	Electron();
	~Electron();

	Real initialData(const RealVect a_pos, const Real a_time) const;
    
	Real m_uniform_density;
	Real m_seed_density;
	Real m_seed_radius;
	RealVect m_seed_pos;
      };

      class PositiveSpecies : public CdrSpecies {
      public:
	PositiveSpecies();
	~PositiveSpecies();

	Real initialData(const RealVect a_pos, const Real a_time) const;

	Real m_uniform_density;
	Real m_seed_density;
	Real m_seed_radius;
	RealVect m_seed_pos;
      };

      class negative_species : public CdrSpecies {
      public:
	negative_species();
	~negative_species();

	Real initialData(const RealVect a_pos, const Real a_time) const;

	Real m_uniform_density;
	Real m_seed_density;
	Real m_seed_radius;
	RealVect m_seed_pos;
      };

      class PhotonOne : public RtSpecies {
      public:
	PhotonOne();
	~PhotonOne();
	Real getAbsorptionCoefficient(const RealVect a_pos) const;
	Real getLambda() const { return m_lambda;}
	Real getA()      const { return m_A;}
	Real getPO2()    const { return m_pO2; }

      protected:
	Real m_lambda;
	Real m_A;
	Real m_pO2;
      };

      class PhotonTwo : public RtSpecies {
      public:
	PhotonTwo();
	~PhotonTwo();
	Real getAbsorptionCoefficient(const RealVect a_pos) const;
	Real getLambda() const { return m_lambda;}
	Real getA()      const { return m_A;}
	Real getPO2()    const { return m_pO2; }

      protected:
	Real m_lambda;
	Real m_A;
	Real m_pO2;
      };

      class PhotonThree : public RtSpecies {
      public:
	PhotonThree();
	~PhotonThree();
	Real getAbsorptionCoefficient(const RealVect a_pos) const;
	Real getLambda() const { return m_lambda;}
	Real getA()      const { return m_A;}
	Real getPO2()    const { return m_pO2; }

      protected:
	Real m_lambda;
	Real m_A;
	Real m_pO2;
      };


      Vector<int> m_wallBc;

      int m_nelec_idx; // Electron index
      int m_nplus_idx; // Positive ion index
      int m_nminu_idx; // Negative ion index


      int m_Photon1_idx;
      int m_Photon2_idx;
      int m_Photon3_idx;
      int m_noise_octaves;
      Real m_background_rate;
      Real m_townsend2_conductor;
      Real m_townsend2_dielectric;
      Real m_electrode_yield;
      Real m_dielectric_yield;
      Real m_dielectric_work;
      Real m_N;
      Real m_fracO2;
      Real m_fracN2;
      Real m_pq;
      Real m_p;
      Real m_exc_eff;
      Real m_photo_eff;
      Real m_temp;


      Vector<Real> computeCathodeFlux(const Vector<Real> a_extrapolated_fluxes,
				      const Vector<Real> a_ion_densities,
				      const Vector<Real> a_ion_velocities,
				      const Vector<Real> a_Photon_fluxes,
				      const RealVect     a_E,
				      const RealVect     a_pos,
				      const RealVect     a_normal,
				      const Real         a_time) const;

      Vector<Real> computeAnodeFlux(const Vector<Real> a_extrapolated_fluxes,
				    const Vector<Real> a_ion_densities,
				    const Vector<Real> a_ion_velocities,
				    const Vector<Real> a_Photon_fluxes,
				    const RealVect     a_E,
				    const RealVect     a_pos,
				    const RealVect     a_normal,
				    const Real         a_time) const;



      Real computeAlpha(const RealVect a_E) const;
      Real computeEta(const RealVect a_E) const;
      Real computeEta2(const RealVect a_E) const;
      Real computeEta3(const RealVect a_E) const;
      Real computeBeta(const RealVect a_E) const;
      Real computeDe(const RealVect a_E) const;
  
      RealVect computeVe(const RealVect a_E) const;
      RealVect computeVp(const RealVect a_E) const;
      RealVect computeVn(const RealVect a_E) const;

      void initSpecies();
  
      void parseGas();
      void parsePhotoi();
      void parseSEE();
      void parseBc();

    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
