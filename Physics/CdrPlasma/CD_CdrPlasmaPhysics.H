/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaPhysics.H
  @brief  Declaration of basic physics coupling functions for Physics/CdrPlasma
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaPhysics_H
#define CD_CdrPlasmaPhysics_H

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
  
    /*!
      @brief Abstract base class for specifying plasma kinetics
    */
    class CdrPlasmaPhysics {
    public:

      /*!
	@brief Constructor
      */
      CdrPlasmaPhysics(){
	m_eed_solve     = false;
	m_eed_index = 0;
      }

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaPhysics(){}

      /*!
	@brief Get number of plot variables for this physics class. 
	@details This is used by CdrPlasmaStepper for pre-allocating data that will be put in a plot file. 
      */
      virtual int getNumberOfPlotVariables() const {
	return 0;
      }

      /*!
	@brief Get plot variable names. The names and positions between this routine and getPlotVariables must be consistent!
      */
      virtual Vector<std::string> getPlotVariableNames() const {
	return Vector<std::string>(this->getNumberOfPlotVariables(), std::string("empty data"));
      }

      /*!
	@brief Provide plot variables. This is used by CdrPlasmaStepper when writing plot files. 
	@details The plot variables should have length this->getNumberOfPlotVariables (as should getPlotVariableNames)
	@param[in] a_position Physical coordinates
	@param[in] a_time     Time. 
      */
      virtual Vector<Real> getPlotVariables(const Vector<Real> a_cdrDensities,
					    const Vector<Real> a_rteDensities,
					    const RealVect     a_E,
					    const RealVect     a_position,
					    const Real         a_time) const {
	return Vector<Real> (this->getNumberOfPlotVariables(), 1.0);
      }

      /*!
	@brief Compute alpha. Should return Townsend ionization coefficient. 
	@details This function is mostly used for the cell tagging classes
	@param[in] a_E        Electric field.
	@param[in] a_position Physical coordinates
      */
      virtual Real computeAlpha(const Real a_E, const RealVect a_position) const  = 0;
      
      /*!
	@brief Routine intended for advancing a reaction network over a time a_dt. 
	@details This routine assumes that the subsequent advance is in the form phi^(k+1) = phi^k + S*a_dt. Thus, this routine exists such that users can EITHER
	fill a_cdrSources and a_rteSources directly with an explicit rule, OR they can perform a fully implicit advance within this routine and set S from that. 
	@param[out] a_cdrSources    Source terms for CDR equations.
	@param[out] a_rteSources    Source terms for RTE equations.
	@param[in]  a_cdrDensities  Grid-based density for particle species.
	@param[in]  a_cdrGradients  Grid-based gradients for particle species.
	@param[in]  a_rteDensities  Grid-based densities for photons.
	@param[in]  a_E             Electric field.
	@param[in]  a_pos           Position in space.
	@param[in]  a_dx            Grid resolution. 
	@param[in]  a_dt            Advanced time.
	@param[in]  a_time          Current time.
	@param[in]  a_kappa         Grid cell unit volume. 
      */
      virtual void advanceReactionNetwork(Vector<Real>&          a_cdrSources,
					  Vector<Real>&          a_rteSources,
					  const Vector<Real>     a_cdrDensities,
					  const Vector<RealVect> a_cdrGradients,
					  const Vector<Real>     a_rteDensities,
					  const RealVect         a_E,
					  const RealVect         a_pos,
					  const Real             a_dx,
					  const Real             a_dt,
					  const Real             a_time,
					  const Real             a_kappa) const  = 0;

      /*!
	@brief Compute velocities for the CDR equations
	@param[in] a_time         Time
	@param[in] a_pos          Position
	@param[in] a_E            Electric field
	@param[in] a_cdrDensities CDR densities
	@return Returns the drift velocities for each CDR species. The vector ordering is the same as m_CdrSpecies. 
      */
      virtual Vector<RealVect> computeCdrDriftVelocities(const Real         a_time,
							 const RealVect     a_pos,
							 const RealVect     a_E,
							 const Vector<Real> a_cdrDensities) const  = 0;

      /*!
	@brief Compute diffusion coefficients for the CDR equations. 
	@param[in] a_time         Time
	@param[in] a_pos          Position
	@param[in] a_E            Electric field
	@param[in] a_cdrDensities CDR densities
	@return Returns the diffusion coefficients for each CDR species. The vector ordering is the same as m_CdrSpecies. 
      */
      virtual Vector<Real> computeCdrDiffusionCoefficients(const Real         a_time,
							   const RealVect     a_pos,
							   const RealVect     a_E,
							   const Vector<Real> a_cdrDensities) const  = 0;


      /*!
	@brief Compute CDR fluxes on electrode-gas interfaces. This is used as a boundary condition in the CDR equations. 
	@param[in] a_time            Time
	@param[in] a_pos             Position
	@param[in] a_normal          Boundary normal vector. This points into the gas phase. 
	@param[in] a_E               Electric field
	@param[in] a_cdrVelocities   CDR velocities. Normal component only. 
	@param[in] a_cdrDensities    CDR densities. 
	@param[in] a_cdrGradients    Normal gradients of cdr densities
	@param[in] a_rteFluxes       RTE fluxes (normal component only)
	@param[in] a_extrapCdrFluxes Extrapolated fluxes from the gas side. 
	@return Returns the flux on an electrode interface cell. The vector ordering must be the same as m_CdrSpecies. 
      */
      virtual Vector<Real> computeCdrElectrodeFluxes(const Real         a_time,
						     const RealVect     a_pos,
						     const RealVect     a_normal,
						     const RealVect     a_E,
						     const Vector<Real> a_cdrDensities,
						     const Vector<Real> a_cdrVelocities,
						     const Vector<Real> a_cdrGradients,
						     const Vector<Real> a_rteFluxes,
						     const Vector<Real> a_extrapCdrFluxes) const  = 0;

      /*!
	@brief Compute CDR fluxes on dielectric-gas interfaces. This is used as a boundary condition in the CDR equations. 
	@param[in] a_time            Time
	@param[in] a_pos             Position
	@param[in] a_normal          Normal vector. This points into the gas phase. 
	@param[in] a_E Electric      field
	@param[in] a_cdrDensities    CDR densities (on the EB)
	@param[in] a_cdrVelocities   Normal component of CDR velocities (on the EB).
	@param[in] a_cdrGradients    Normal gradients of cdr densities
	@param[in] a_rteFluxes       RTE fluxes (normal component only)
	@param[in] a_extrapCdrFluxes Extrapolated fluxes from the gas side. 
	@return Returns the flux on a dielectric interface cell. The vector ordering must be the same as m_CdrSpecies. 
      */
      virtual Vector<Real> computeCdrDielectricFluxes(const Real         a_time,
						      const RealVect     a_pos,
						      const RealVect     a_normal,
						      const RealVect     a_E,
						      const Vector<Real> a_cdrDensities,
						      const Vector<Real> a_cdrVelocities,
						      const Vector<Real> a_cdrGradients,
						      const Vector<Real> a_rteFluxes,
						      const Vector<Real> a_extrapCdrFluxes) const  = 0;

      /*!
	@brief Compute CDR fluxes through domain sides. This is used as a boundary condition in the CDR equations.
	@param[in] a_time            Time
	@param[in] a_pos             Position
	@param[in] a_dir             Direction (0 = x, 1=y etc)
	@param[in] a_side            Side (low or high side)
	@param[in] a_E               Electric field
	@param[in] a_cdrDensities    CDR densities. 
	@param[in] a_cdrVelocities   CDR velocities (normal component only). 
	@param[in] a_cdrGradients    CDR gradients (normal component only)
	@param[in] a_rteFluxes       RTE fluxes (normal component only)
	@param[in] a_extrapCdrFluxes Extrapolated fluxes from the gas side. 
      */
      virtual Vector<Real> computeCdrDomainFluxes(const Real           a_time,
						  const RealVect       a_pos,
						  const int            a_dir,
						  const Side::LoHiSide a_side,
						  const RealVect       a_E,
						  const Vector<Real>   a_cdrDensities,
						  const Vector<Real>   a_cdrVelocities,
						  const Vector<Real>   a_cdrGradients,
						  const Vector<Real>   a_rteFluxes,
						  const Vector<Real>   a_extrapCdrFluxes) const  = 0;

      /*!
	@brief Set the initial surface charge
	@param[in] a_time Time
	@param[in] a_pos Position
      */
      virtual Real initialSigma(const Real a_time, const RealVect a_pos) const  = 0;      

      /*!
	@brief Initialize electron energy density
      */
      virtual Real initEed(const RealVect a_pos, const Real a_time, const RealVect a_E){
	return 0.0;
      }

      /*!
	@brief Get all ions
      */
      const Vector<RefCountedPtr<CdrSpecies> >& getCdrSpecies() const { 
	return m_CdrSpecies; 
      }

      /*!
	@brief Get all Photons
      */
      const Vector<RefCountedPtr<RtSpecies> >& getRtSpecies() const {
	return m_RtSpecies;
      }

      /*!
	@brief Return number of ion equations
      */
      int getNumCdrSpecies() const{
	return m_CdrSpecies.size();
      }

      /*!
	@brief Return number of RTE equations
      */
      int getNumRtSpecies() const {
	return m_RtSpecies.size();
      }

      /*!
	@brief Check if solver solves for eed
      */
      bool solveEed(){
	return m_eed_solve;
      }

      /*!
	@brief Get eed solver index
      */
      int getEedIndex(){
	return m_eed_index;
      }

      void setDt(const Real a_dt){
	m_dt = a_dt;
      }

      void setFinestDx(const Real a_finest_dx){
	m_finest_dx = a_finest_dx;
      }

  
    protected:
      bool m_eed_solve;     // Doing EED solve or not?
  
      int m_eed_index;      // EED solver index
  
      Real m_dt;            // FHD kinetics need access to time steps
      Real m_finest_dx;     // Finest dx

      /*!
	@brief List of species
      */
      Vector<RefCountedPtr<CdrSpecies> > m_CdrSpecies;

      /*!
	@brief List of optical transitions between species
      */
      Vector<RefCountedPtr<RtSpecies> > m_RtSpecies;

      /*!
	@brief Number of species
      */
      int m_numCdrSpecies;

      /*!
	@brief Number of RTE species
      */
      int m_numRtSpecies;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
