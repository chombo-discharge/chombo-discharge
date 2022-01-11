/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaStepper.H
  @brief  Declaration of a CdrPlasmaStepper, an abstract class for evolving cdr plasma equations in time across an amr mesh
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaStepper_H
#define CD_CdrPlasmaStepper_H

// Std includes
#include <functional>

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_CdrPlasmaPhysics.H>
#include <CD_AmrMesh.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_CdrLayout.H>
#include <CD_RtLayout.H>
#include <CD_FieldSolver.H>
#include <CD_SigmaSolver.H>
#include <CD_TimeStepper.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {
  
    /*!
      @brief Abstract class for evolving plasma equations in time
    */
    class CdrPlasmaStepper : public TimeStepper {
    protected:

      /*!
	@brief Just a silly namespace for centering stuff...
      */
      enum class Centering {
	CellCenter,
	CellCentroid,
	EbCentroid,
      };
      
    public:

      /*!
	@brief Enum for switching between source term computations. 'Interpolated' interpolates the cell-centered data
	to cell centroids. Kappa-weighted assumes centroid-centered data, and Upwind uses the upwind formulation in Villa et. al. 
      */
      enum class SourceTermComputation {
	Interpolated,
	InterpolatedStable,	
	CellAverage,
	Upwind
      };

      /*!
	@brief Default function that is passed to mapped boundary conditions
      */
      static Real s_constant_one(const RealVect a_pos);

      // New functions for new Driver interface
#ifdef CH_USE_HDF5
      virtual void writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const;
#endif
#ifdef CH_USE_HDF5      
      virtual void readCheckpointData(HDF5Handle& a_handle, const int a_lvl);
#endif
      virtual int getNumberOfPlotVariables() const;
      virtual void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const;
      virtual void writeJ(EBAMRCellData& a_output, int& a_comp) const;
      virtual void writePhysics(EBAMRCellData& a_output, int& a_comp) const;      
      virtual void postCheckpointSetup();
      virtual void printStepReport() override;
      virtual void registerOperators() override;
      virtual void allocate() override;
      virtual void postInitialize() override;
      virtual void registerRealms() override;
      virtual void postRegrid() override;

      /*!
	@brief Constructor
      */
      CdrPlasmaStepper();

      /*!
	@brief Constructor
      */
      CdrPlasmaStepper(RefCountedPtr<CdrPlasmaPhysics>& a_physics);

      /*!
	@brief Destructor
      */
      ~CdrPlasmaStepper();

      /*!
	@brief Advance method, advances equations.
	@param[in] a_dt Time step
	@returns Actual time step used (can be different from a_dt for adaptive methods)
      */
      virtual Real advance(const Real a_dt) = 0;

      /*!
	@brief Check if RTE solvers are stationary
      */
      virtual bool stationary_rte();

      /*!
	@brief Allocate internal storage
      */
      virtual void allocateInternals() = 0;

      /*!
	@brief Deallocate internals
      */
      virtual void deallocateInternals() = 0;

      /*!
	@brief Parse options.
      */
      virtual void parseOptions() = 0;

      /*!
	@brief Parse runtime options
      */
      virtual void parseRuntimeOptions() = 0;

      /*!
	@brief New caching function for regridding
      */
      virtual void preRegrid(const int a_lbase, const int a_finestLevel);

      /*!
	@brief preRegrid internal storage
      */
      virtual void preRegridInternals(const int a_lbase, const int a_finestLevel);

      /*!
	@brief Deallocation function
      */
      virtual void deallocate(){
	this->deallocateInternals();
	this->deallocateSolverInternals();
      }

      /*!
	@brief Compute the cell-centered conductivity
      */
      virtual void computeCellConductivity(EBAMRCellData& a_cellConductivity) const;

      /*!
	@brief Compute the face-centered conductivity from a cell-centered conductivity.
      */
      virtual void computeFaceConductivity(EBAMRFluxData&       a_conductivityFace,
					   EBAMRIVData&         a_conductivityEB,
					   const EBAMRCellData& a_conductivityCell) const;

      /*!
	@brief Set up a semi-implicit Poisson solver.
	@details This is the version that computes the conductivity first.
      */
      virtual void setupSemiImplicitPoisson(const Real a_dt);

      /*!
	@brief Set up a semi-implicit Poisson solver.
	@details This is the version that uses pre-existing conductivities.
	@param[in] a_factor Factor to multiply conductivities with (they could be pre-scaled by dt/eps0 in which case a_factor = 1).
      */
      virtual void setupSemiImplicitPoisson(const EBAMRFluxData& a_conductivityFace,
					    const EBAMRIVData&   a_conductivityEB,
					    const Real           a_factor);

      /*!
	@brief Solve the poisson equation. This calls computeSpaceChargeDensity() 
	@details This assumes that the data contain in CdrLayout is cell-centered
      */
      virtual bool solvePoisson();

      /*!
	@brief General poisson solver routine. 
	@details You have the option of supplying storage for a right-hand side if you want to. The densities can be either cell-centered or centroid-centered, If they are cell-centered, we will extrapolate the right-hand side to the centroids first. 
      */
      virtual bool solvePoisson(MFAMRCellData&                a_potential,
				MFAMRCellData&                a_rhs,
				const Vector<EBAMRCellData*>  a_densities,
				const EBAMRIVData&            a_sigma,
				const Centering               a_centering = Centering::CellCenter);

      /*!
	@brief Advance the reaction network. 
      */
      virtual void advanceReactionNetwork(const Real a_time, const Real a_dt);
  
      /*!
	@brief Compute reaction network sources. First computes the CDR gradients and then calls the other version
      */
      virtual void advanceReactionNetwork(Vector<EBAMRCellData*>&       a_particle_sources,
					  Vector<EBAMRCellData*>&       a_Photon_sources,
					  const Vector<EBAMRCellData*>& a_particle_densities,
					  const Vector<EBAMRCellData*>& a_Photon_densities,
					  const EBAMRCellData&          a_E,
					  const Real&                   a_time,
					  const Real&                   a_dt);

      /*!
	@brief Compute reaction network sources. Used for time-dependent FHD models
      */
      virtual void advanceReactionNetwork(Vector<EBAMRCellData*>&       a_particle_sources,
					  Vector<EBAMRCellData*>&       a_Photon_sources,
					  const Vector<EBAMRCellData*>& a_particle_densities,
					  const Vector<EBAMRCellData*>& a_particle_gradients,
					  const Vector<EBAMRCellData*>& a_Photon_densities,
					  const EBAMRCellData&          a_E,
					  const Real&                   a_time,
					  const Real&                   a_dt);


      virtual void advanceReactionNetwork(Vector<LevelData<EBCellFAB>* >&       a_particle_sources,
					  Vector<LevelData<EBCellFAB>* >&       a_Photon_sources,
					  const Vector<LevelData<EBCellFAB>* >& a_particle_densities,
					  const Vector<LevelData<EBCellFAB>* >& a_particle_gradients,
					  const Vector<LevelData<EBCellFAB>* >& a_Photon_densities,
					  const LevelData<EBCellFAB>&           a_E,
					  const Real&                           a_time,
					  const Real&                           a_dt,
					  const int                             a_level);
      /*!
	@brief Compute reaction network sources. Patch version.
      */
      virtual void advanceReactionNetworkRegularCells(Vector<EBCellFAB*>&       a_particle_sources,
						      Vector<EBCellFAB*>&       a_Photon_sources,
						      const Vector<EBCellFAB*>& a_particle_densities,
						      const Vector<EBCellFAB*>& a_particle_gradients,
						      const Vector<EBCellFAB*>& a_Photon_densities,
						      const EBCellFAB&          a_E,
						      const Real&               a_time,
						      const Real&               a_dt,
						      const Real&               a_dx,
						      const Box&                a_box);

      /*!
	@brief Compute reaction network sources. Fast patch version. 
      */
      virtual void advanceReactionNetworkRegularCellsFast(Vector<EBCellFAB*>&       a_particle_sources,
							  Vector<EBCellFAB*>&       a_Photon_sources,
							  const Vector<EBCellFAB*>& a_particle_densities,
							  const Vector<EBCellFAB*>& a_particle_gradients,
							  const Vector<EBCellFAB*>& a_Photon_densities,
							  const EBCellFAB&          a_E,
							  const Real&               a_time,
							  const Real&               a_dt,
							  const Real&               a_dx,
							  const Box&                a_box);

      virtual void advanceReactionNetworkRegularCellsFast2D(Vector<EBCellFAB*>&       a_particle_sources,
							    Vector<EBCellFAB*>&       a_Photon_sources,
							    const Vector<EBCellFAB*>& a_particle_densities,
							    const Vector<EBCellFAB*>& a_particle_gradients,
							    const Vector<EBCellFAB*>& a_Photon_densities,
							    const EBCellFAB&          a_E,
							    const Real&               a_time,
							    const Real&               a_dt,
							    const Real&               a_dx,
							    const Box&                a_box);

      virtual void advanceReactionNetworkRegularCellsFast3D(Vector<EBCellFAB*>&       a_particle_sources,
							    Vector<EBCellFAB*>&       a_Photon_sources,
							    const Vector<EBCellFAB*>& a_particle_densities,
							    const Vector<EBCellFAB*>& a_particle_gradients,
							    const Vector<EBCellFAB*>& a_Photon_densities,
							    const EBCellFAB&          a_E,
							    const Real&               a_time,
							    const Real&               a_dt,
							    const Real&               a_dx,
							    const Box&                a_box);

      /*!
	@brief Compute reaction network sources. Irregular patch version.
      */
      virtual void advanceReactionNetworkIrreg(Vector<EBCellFAB*>&          a_particle_sources,
					       Vector<EBCellFAB*>&          a_Photon_sources,
					       const Vector<EBCellFAB*>&    a_particle_densities,
					       const Vector<EBCellFAB*>&    a_particle_gradients,
					       const Vector<EBCellFAB*>&    a_particle_velocities,
					       const Vector<EBCellFAB*>&    a_Photon_densities,
					       const BaseIVFAB<VoFStencil>& a_interp_stencils,
					       const EBCellFAB&             a_E,
					       const Real&                  a_time,
					       const Real&                  a_dt,
					       const Real&                  a_dx,
					       const Box&                   a_box,
					       const int                    a_lvl,
					       const DataIndex&             a_dit);

      /*!
	@brief Compute reaction network sources. Irregular patch version.
      */
      virtual void advanceReactionNetworkIrregInterp(Vector<EBCellFAB*>&          a_particle_sources,
						     Vector<EBCellFAB*>&          a_Photon_sources,
						     const Vector<EBCellFAB*>&    a_particle_densities,
						     const Vector<EBCellFAB*>&    a_particle_gradients,
						     const Vector<EBCellFAB*>&    a_particle_velocities,
						     const Vector<EBCellFAB*>&    a_Photon_densities,
						     const BaseIVFAB<VoFStencil>& a_interp_stencils,
						     const EBCellFAB&             a_E,
						     const Real&                  a_time,
						     const Real&                  a_dt,
						     const Real&                  a_dx,
						     const Box&                   a_box,
						     const int                    a_lvl,
						     const DataIndex&             a_dit);

      virtual void advanceReactionNetworkIrregKappa(Vector<EBCellFAB*>&          a_particle_sources,
						    Vector<EBCellFAB*>&          a_Photon_sources,
						    const Vector<EBCellFAB*>&    a_particle_densities,
						    const Vector<EBCellFAB*>&    a_particle_gradients,
						    const Vector<EBCellFAB*>&    a_Photon_densities,
						    const BaseIVFAB<VoFStencil>& a_interp_stencils,
						    const EBCellFAB&             a_E,
						    const Real&                  a_time,
						    const Real&                  a_dt,
						    const Real&                  a_dx,
						    const Box&                   a_box,
						    const int                    a_lvl,
						    const DataIndex&             a_dit);
      

      /*!
	@brief Compute reaction network sources. Uses upwinding of reaction terms. 
      */
      virtual void advanceReactionNetworkIrregUpwind(Vector<EBCellFAB*>&          a_particle_sources,
						     Vector<EBCellFAB*>&          a_Photon_sources,
						     const Vector<EBCellFAB*>&    a_particle_densities,
						     const Vector<EBCellFAB*>&    a_particle_gradients,
						     const Vector<EBCellFAB*>&    a_particle_velocities,
						     const Vector<EBCellFAB*>&    a_Photon_densities,
						     const BaseIVFAB<VoFStencil>& a_interp_stencils,
						     const EBCellFAB&             a_E,
						     const Real&                  a_time,
						     const Real&                  a_dt,
						     const Real&                  a_dx,
						     const Box&                   a_box,
						     const int                    a_lvl,
						     const DataIndex&             a_dit);      


      /*!
	@brief Compute the time step and how it was restricted
      */
      virtual void computeDt(Real& a_dt, TimeCode& a_timeCode) = 0;

      /*!
	@brief Compute the cell-centered electric field on both phases
      */
      virtual void computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential);

      /*!
	@brief Compute the cell-centered electric field on a specific phase using whatever is available in m_fieldSolver
      */
      virtual void computeElectricField(EBAMRCellData& a_E, const phase::which_phase a_phase);

      /*!
	@brief Compute the cell-centered electric field on a specific phase
      */
      virtual void computeElectricField(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential);

      /*!
	@brief Compute the face-centered electric field (only normal components)
      */
      virtual void computeElectricField(EBAMRFluxData& a_E_face, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell);
  
      /*!
	@brief Compute the face-centered electric field (only normal components)
      */
      virtual void computeElectricField(EBAMRIVData& a_E_eb, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell);

      /*!
	@brief Compute the maximum of the electric field
      */
      virtual void computeMaxElectricField(Real& a_Emax, const phase::which_phase a_phase);

      /*!
	@brief Compute the current density
      */
      virtual void computeJ(EBAMRCellData& a_J) const;

      /*!
	@brief Compute the current density. Level version. 
      */
      virtual void computeJ(LevelData<EBCellFAB>& a_J, const int a_lvl) const;

  
      /*!
	@brief Compute diffusion things using whatever is available in the solvers
      */
      virtual void computeCdrDiffusion();

      /*!
	@brief Compute diffusion things using whatever is available in the solvers, except for the fields
      */
      virtual void computeCdrDiffusion(const EBAMRCellData& a_E_cell, const EBAMRIVData& a_E_eb);

      /*!
	@brief Compute diffusion coefficients
	@details The way we do this is that we compute the diffusion coefficients on cell centers, and then we average that
	onto cell faces. 
      */
      virtual void computeCdrDiffusionFace(Vector<EBAMRFluxData*>&       a_diffusionCoefficient_face,
					   const Vector<EBAMRCellData*>& a_cdr_densities,
					   const EBAMRCellData&          a_E,
					   const Real&                   a_time);

      /*!
	@brief Compute diffusion coefficients on cell centers. Calls the level version. 
      */
      virtual void computeCdrDiffusionCell(Vector<EBAMRCellData>&        a_diffusionCoefficient_cell,
					   const Vector<EBAMRCellData*>& a_cdr_densities,
					   const EBAMRCellData&          a_E,
					   const Real&                   a_time);

      /*!
	@brief Compute diffusion coefficients on cell centers. Regular cells. 
      */
      virtual void computeCdrDiffusionCellRegular(Vector<EBCellFAB*>&       a_diffusionCoefficient_cell,
						  const Vector<EBCellFAB*>& a_cdr_densities,
						  const EBCellFAB&          a_E,
						  const Box                 a_box,
						  const Real                a_dx,
						  const Real                a_time);

      /*!
	@brief Compute diffusion coefficients on cell centers. Fast version that calls dimension versions
      */
      virtual void computeCdrDiffusionCellRegularFast(Vector<EBCellFAB*>&       a_diffusionCoefficient_cell,
						      const Vector<EBCellFAB*>& a_cdr_densities,
						      const EBCellFAB&          a_E,
						      const Box                 a_box,
						      const Real                a_dx,
						      const Real                a_time);

      /*!
	@brief Compute diffusion coefficients on cell centers. Fast version that calls dimension versions
      */
      virtual void computeCdrDiffusionCellRegularFast2D(Vector<EBCellFAB*>&       a_diffusionCoefficient_cell,
							const Vector<EBCellFAB*>& a_cdr_densities,
							const EBCellFAB&          a_E,
							const Box                 a_box,
							const Real                a_dx,
							const Real                a_time);

      /*!
	@brief Compute diffusion coefficients on cell centers. Fast version that calls dimension versions
      */
      virtual void computeCdrDiffusionCellRegularFast3D(Vector<EBCellFAB*>&       a_diffusionCoefficient_cell,
							const Vector<EBCellFAB*>& a_cdr_densities,
							const EBCellFAB&          a_E,
							const Box                 a_box,
							const Real                a_dx,
							const Real                a_time);

      /*!
	@brief Compute diffusion coefficients on cell centers. Regular cells. 
      */
      virtual void computeCdrDiffusionCellIrregular(Vector<EBCellFAB*>&           a_diffusionCoefficient_cell,
						    const Vector<EBCellFAB*>&    a_cdr_densities,
						    const EBCellFAB&             a_E,
						    const Box                    a_box,
						    const Real                   a_dx,
						    const BaseIVFAB<VoFStencil>& a_interp_stencils,
						    const Real&                  a_time,
						    const int                    a_lvl,
						    const DataIndex&             a_dit);
  
      /*!
	@brief Compute diffusion coefficients. Level version. 
      */
      virtual void computeCdrDiffusionCell(Vector<LevelData<EBCellFAB>* >&       a_diffusionCoefficient_cell,
					   const Vector<LevelData<EBCellFAB>* >& a_cdr_densities,
					   const LevelData<EBCellFAB>&           a_E,
					   const int                             a_lvl,
					   const Real&                           a_time);

      /*!
	@brief Compute diffusion coefficients on the EB. 
      */
      virtual void computeCdrDiffusionEb(Vector<EBAMRIVData*>&       a_diffusionCoefficient_face,
					 const Vector<EBAMRIVData*>& a_cdr_densities,
					 const EBAMRIVData&          a_E,
					 const Real&                 a_time);

      /*!
	@brief Compute diffusion coefficients on the EB. Level version. 
      */
      virtual void computeCdrDiffusionEb(Vector<LevelData<BaseIVFAB<Real> >* >&       a_ebDiffusionCoefficient,
					 const Vector<LevelData<BaseIVFAB<Real> >* >& a_cdr_densities,
					 const LevelData<BaseIVFAB<Real> >&           a_E,
					 const Real&                                  a_time,
					 const int                                    a_lvl);

      /*!
	@brief Compute cdr_fluxes
      */
      virtual void computeCdrFluxes(Vector<LevelData<BaseIVFAB<Real> >*>&       a_fluxes,
				    const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_cdr_fluxes,
				    const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_cdr_densities,
				    const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_cdr_velocities,
				    const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_cdr_gradients,
				    const Vector<LevelData<BaseIVFAB<Real> >*>& a_extrap_rte_fluxes,
				    const LevelData<BaseIVFAB<Real> >&          a_field,
				    const Real&                                 a_time,
				    const int                                   a_lvl);

      virtual void computeCdrFluxes(Vector<EBAMRIVData*>&       a_fluxes,
				    const Vector<EBAMRIVData*>& a_extrap_cdr_fluxes,
				    const Vector<EBAMRIVData*>& a_extrap_cdr_densities,
				    const Vector<EBAMRIVData*>& a_extrap_cdr_velocities,
				    const Vector<EBAMRIVData*>& a_extrap_cdr_gradients,
				    const Vector<EBAMRIVData*>& a_extrap_rte_fluxes,
				    const EBAMRIVData&          a_field,
				    const Real&                 a_time);

      /*!
	@brief Compute cdr fluxes on domain faces
	@note Everything that comes in here should have been computed on the domain faces
      */
      virtual void computeCdrDomainFluxes(Vector<EBAMRIFData*>&       a_fluxes,
					  const Vector<EBAMRIFData*>& a_extrap_cdr_fluxes,
					  const Vector<EBAMRIFData*>& a_extrap_cdr_densities,
					  const Vector<EBAMRIFData*>& a_extrap_cdr_velocities,
					  const Vector<EBAMRIFData*>& a_extrap_cdr_gradients,
					  const Vector<EBAMRIFData*>& a_extrap_rte_fluxes,
					  const EBAMRIFData&          a_field,
					  const Real&                 a_time);

      /*!
	@brief Compute cdr fluxes on domain faces. Level version
      */
      virtual void computeCdrDomainFluxes(Vector<LevelData<DomainFluxIFFAB>*>        a_fluxes,
					  const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_cdr_fluxes,
					  const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_cdr_densities,
					  const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_cdr_velocities,
					  const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_cdr_gradients,
					  const Vector<LevelData<DomainFluxIFFAB>*>& a_extrap_rte_fluxes,
					  const LevelData<DomainFluxIFFAB>&          a_E,
					  const Real&                                a_time,
					  const int                                  a_lvl);

      /*!
	@brief Compute the gradient at the boundary (we also take the dot product w.r.t. to the normal vector)
	@note This version computes the gradients from phi
      */
      virtual void computeGradientsAtEb(Vector<EBAMRIVData*>&         a_grad,
					const phase::which_phase&     a_phase,
					const Vector<EBAMRCellData*>& a_phi);

      /*!
	@brief Compute gradients at domain faces
      */
      virtual void computeGradientsAtDomainFaces(Vector<EBAMRIFData*>&         a_grad,
						 const phase::which_phase&     a_phase,
						 const Vector<EBAMRCellData*>& a_phi);

      /*!
	@brief Compute cdr velocities
      */
      virtual void computeCdrDriftVelocities();

      /*!
	@brief Compute the cell-centered cdr velocities. Calls the level versions. 
      */
      virtual void computeCdrDriftVelocities(Vector<EBAMRCellData*>&       a_velocities,
					     const Vector<EBAMRCellData*>& a_cdr_densities,
					     const EBAMRCellData&          a_E,
					     const Real&                   a_time);

      /*!
	@brief Compute the cell-centered cdr velocities. Level version. 
      */
      virtual void computeCdrDriftVelocities(Vector<LevelData<EBCellFAB> *>&       a_velocities,
					     const Vector<LevelData<EBCellFAB> *>& a_cdr_densities,
					     const LevelData<EBCellFAB> &          a_E,
					     const int                             a_lvl,
					     const Real&                           a_time);

      virtual void computeCdrDriftVelocitiesRegular(Vector<EBCellFAB*>&       a_velocities,
						    const Vector<EBCellFAB*>& a_cdr_densities,
						    const EBCellFAB&          a_E,
						    const Box&                a_box,
						    const Real&               a_time,
						    const Real&               a_dx);

      virtual void computeCdrDriftVelocitiesRegularFast(Vector<EBCellFAB*>&       a_velocities,
							const Vector<EBCellFAB*>& a_cdr_densities,
							const EBCellFAB&          a_E,
							const Box&                a_box,
							const Real&               a_time,
							const Real&               a_dx);

      virtual void computeCdrDriftVelocitiesRegularFast2D(Vector<EBCellFAB*>&       a_velocities,
							  const Vector<EBCellFAB*>& a_cdr_densities,
							  const EBCellFAB&          a_E,
							  const Box&                a_box,
							  const Real&               a_time,
							  const Real&               a_dx);
  
      virtual void computeCdrDriftVelocitiesRegularFast3D(Vector<EBCellFAB*>&       a_velocities,
							  const Vector<EBCellFAB*>& a_cdr_densities,
							  const EBCellFAB&          a_E,
							  const Box&                a_box,
							  const Real&               a_time,
							  const Real&               a_dx);

      virtual void computeCdrDriftVelocitiesIrregular(Vector<EBCellFAB*>&       a_velocities,
						      const Vector<EBCellFAB*>& a_cdr_densities,
						      const EBCellFAB&          a_E,
						      const Box&                a_box,
						      const Real&               a_time,
						      const Real&               a_dx,
						      const int                 a_lvl,
						      const DataIndex&          a_dit);

      /*!
	@brief Compute the cell-centered fluxes and extrapolate these to the EB centroid
      */
      virtual void computeExtrapolatedFluxes(Vector<EBAMRIVData*>&        a_fluxes,
					     const Vector<EBAMRCellData*> a_densities,
					     const Vector<EBAMRCellData*> a_velocities,
					     const phase::which_phase     a_phase);

      /*!
	@brief Compute the cell-centered fluxes and extrapolate these to the EB centroid
      */
      virtual void computeExtrapolatedVelocities(Vector<EBAMRIVData*>&        a_ebvel,
						 const Vector<EBAMRCellData*> a_velocities,
						 const phase::which_phase     a_phase);

      /*!
	@brief Compute the cell-centered fluxes and extrapolate these to the domain faces
      */
      virtual void computeExtrapolatedDomainFluxes(Vector<EBAMRIFData*>&        a_fluxes,
						   const Vector<EBAMRCellData*> a_densities,
						   const Vector<EBAMRCellData*> a_velocities,
						   const phase::which_phase     a_phase);


      /*!
	@brief Compute flux
      */
      virtual void computeFlux(EBAMRCellData&       a_flux,
			       const EBAMRCellData& a_density,
			       const EBAMRCellData& a_velocity);

      virtual void computeFlux(LevelData<EBCellFAB>&       a_flux,
			       const LevelData<EBCellFAB>& a_density,
			       const LevelData<EBCellFAB>& a_velocity,
			       const int                   a_lvl);

      /*!
	@brief Compute charge flux. 
      */
      virtual void computeChargeFlux(EBAMRIVData& a_flux, Vector<EBAMRIVData*>& a_cdr_flxues);


      /*!
	@brief Compute the cell-centered space charge density by using data inside CdrLayout. 
	@details This assumes that the data contain in CdrLayout is cell-centered
      */
      virtual void computeSpaceChargeDensity();

      /*!
	@brief Compute cell-centered space charge density on a specific phase using whatever is available in CdrLayout
      */
      virtual void computeSpaceChargeDensity(EBAMRCellData& a_rho, const phase::which_phase a_phase);

      /*!
	@brief Compute the centroid-centered space charge density
	@details You have the option of specifying the centering for this computation. If the densities are cell-centered, we first 
	compute the cell-centered space charge density and the interpolate to the cell centroid. 
      */
      virtual void computeSpaceChargeDensity(MFAMRCellData&                 a_rho,
					     const Vector<EBAMRCellData*>&  a_densities,
					     const Centering                a_centering = Centering::CellCenter);


  
      /*!
	@brief Deallocate internal solver storages
      */
      virtual void deallocateSolverInternals();

      virtual void extrapolateToEb(LevelData<BaseIVFAB<Real> >& a_extrap,
				   const phase::which_phase     a_phase,
				   const LevelData<EBCellFAB>&  a_data,
				   const int                    a_lvl);

      /*!
	@brief Extrapolate data to EB centroid
      */
      virtual void extrapolateToEb(EBAMRIVData& a_extrap, const phase::which_phase a_phase, const EBAMRCellData& a_data);

      /*!
	@brief Extrapolate data to EB centroid
      */
      virtual void extrapolateToEb(Vector<EBAMRIVData*>&         a_extrap,
				   const phase::which_phase      a_phase,
				   const Vector<EBAMRCellData*>& a_data);

      /*!
	@brief Extrapolate to domain faces. Calls level version. 
      */
      virtual void extrapolateToDomainFaces(EBAMRIFData& a_extrap, const phase::which_phase a_phase, const EBAMRCellData& a_data);

      /*!
	@brief Extrapolate to domain faces
      */
      virtual void extrapolateToDomainFaces(LevelData<DomainFluxIFFAB>& a_extrap,
					    const phase::which_phase    a_phase,
					    const LevelData<EBCellFAB>& a_data,
					    const int                   a_lvl);

      /*!
	@brief Extrapolate to domain faces
      */
      virtual void extrapolateToDomainFaces(Vector<EBAMRIFData*>&         a_extrap,
					    const phase::which_phase      a_phase,
					    const Vector<EBAMRCellData*>& a_data);
  

      virtual void extrapolateToVectorDomainFaces(EBAMRIFData&             a_extrap,
						  const phase::which_phase a_phase,
						  const EBAMRCellData&     a_data);

      virtual void extrapolateToVectorDomainFaces(Vector<EBAMRIFData*>&         a_extrap,
						  const phase::which_phase      a_phase,
						  const Vector<EBAMRCellData*>& a_data);

      virtual void extrapolateVelocitiesVectorDomainFaces(Vector<EBAMRIFData*>&         a_extrap,
							  const phase::which_phase      a_phase,
							  const Vector<EBAMRCellData*>& a_data);

      /*!
	@brief Get maximum value of the cdr solver stuff
      */
      virtual void getCdrMax(Real& a_cdr_max, std::string& a_solver_name);

      /*!
	@brief Instantiate solvers
      */
      virtual void setupSolvers();

      /*!
	@brief Fill all solvers with initial data
      */
      virtual void initialData();

      /*!
	@brief Initialize sigma
      */
      virtual void initialSigma();

      /*!
	@brief Project flux
      */
      virtual void projectFlux(LevelData<BaseIVFAB<Real> >& a_proj_flux, const LevelData<BaseIVFAB<Real> >& a_flux, const int a_lvl);

      /*!
	@brief Extrapolated flux to EB and project it onto the normal
      */
      virtual void projectFlux(EBAMRIVData& a_projected_flux, const EBAMRIVData& a_flux);

      /*!
	@brief Project flux on the boundary normal
      */
      virtual void projectDomain(EBAMRIFData& a_projected_flux, const EBAMRIFData& a_flux);

      /*!
	@brief Regrid
      */
      virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

      /*!
	@brief Regrid internal storage. This MUST be overwritten.
      */
      virtual void regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) = 0;

      /*!
	@brief Regrid individual solvers
      */
      virtual void regridSolvers(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

      /*!
	@brief Set data defined over dielectric cells to zero
      */
      virtual void resetDielectricCells(EBAMRIVData& a_data);

      /*!
	@brief Do a sanity check
      */
      virtual void sanityCheck();

      /*!
	@brief Set fast rte
      */
      virtual void setFastRadiativeTransfer(const int a_fast_rte);

      /*!
	@brief Set fast rte
      */
      virtual void setFastPoisson(const int a_fast_poisson);

      /*!
	@brief Set hardcap on the time step
      */
      virtual void setMinDt(const Real a_min_dt);

      /*!
	@brief Set hardcap on the time step
      */
      virtual void setMaxDt(const Real a_min_dt);

      /*!
	@brief Set CFL condition
      */
      virtual void setCFL(const Real a_cfl);

      /*!
	@brief Set fudge factor for relaxation time
      */
      virtual void setRelaxTime(const Real a_relax_time);

      /*!
	@brief Set source growth
      */
      virtual void setSourceGrowth(const Real a_src_growth);

      /*!
	@brief Set source growth
      */
      virtual void setSourceGrowthTolerance(const Real a_src_tolerance);

      /*!
	@brief Set plasma kinetics
      */
      virtual void setCdrPlasmaPhysics(const RefCountedPtr<CdrPlasmaPhysics>& a_physics);

      /*!
	@brief Set solver verbosity
      */
      virtual void setSolverVerbosity();

      /*!
	@brief Setup the cdr solvers
      */
      virtual void setupCdr();

      /*!
	@brief Setup the poisson solver
      */
      virtual void setupPoisson();

      /*!
	@brief Setup the rte solvers
      */
      virtual void setupRadiativeTransfer();

      /*!
	@brief Instantiate the sigma solver
      */
      virtual void setupSigma();

      /*!
	@brief Solver dumps data files
      */
      virtual void solverDump();

      /*!
	@brief Set potential
      */
      virtual void setVoltage(std::function<Real(const Real a_time)> a_potential);

      /*!
	@brief RTE advance using data existing in other solvers
      */
      virtual void solveRadiativeTransfer(const Real a_dt);

      /*!
	@brief Set the cdr layout
      */
      virtual void setCdrSolvers(RefCountedPtr<CdrLayout<CdrSolver> >& a_cdr);

      /*!
	@brief Set the Poisson solver
      */
      virtual void setFieldSolver(RefCountedPtr<FieldSolver>& a_poisson);

      /*!
	@brief Set the RTE layout
      */
      virtual void setRadiativeTransferSolvers(RefCountedPtr<RtLayout<RtSolver> >& a_rte);

      /*!
	@brief General RTE advance routine
	@details This computes the source terms at the specified time and then advances the solution
      */
      virtual void solveRadiativeTransfer(Vector<EBAMRCellData*>&       a_phis,
					  Vector<EBAMRCellData*>&       a_rhs,
					  const Vector<EBAMRCellData*>& a_cdr_states,
					  const EBAMRCellData&          a_E,
					  const Real                    a_time,
					  const Real                    a_dt,
					  const Centering               a_centering = Centering::CellCenter);

      /*!
	@brief Synchronize solver times
      */
      virtual void synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt);

      /*!
	@brief Query about necessary number of ghost cells. This collects that from the individual solvers. 
      */
      virtual int queryGhost();

      /*!
	@brief Compute the ohmic current on electrodes
      */
      virtual Real computeElectrodeCurrent();

      /*!
	@brief Compute the ohmic current on dielectrics
      */
      virtual Real computeDielectricCurrent();

      /*!
	@brief Compute the ohmic current through the domain faces
      */
      virtual Real computeDomainCurrent();

      /*!
	@brief Compute induced current in external circuit due to ohmic conduction
      */
      virtual Real computeOhmicInductionCurrent();
			     
      /*!
	@brief Compute the relaxation time
      */
      virtual Real computeRelaxationTime();

      /*!
	@brief Restrict dt
      */
      virtual Real restrictDt() = 0;

      /*!
	@brief Return time
      */
      virtual Real getTime();

      /*!
	@brief Get dt
      */
      virtual Real getDt();

      /*!
	@brief Get the CFL time step
      */
      virtual Real getCflDt();

      /*!
	@brief Get cdr solvers
      */
      virtual RefCountedPtr<CdrLayout<CdrSolver> >& getCdrSolvers();

      /*!
	@brief Get cdr solvers
      */
      virtual RefCountedPtr<FieldSolver>& getFieldSolver();
  
      /*!
	@brief Get cdr solvers
      */
      virtual RefCountedPtr<RtLayout<RtSolver>>& getRadiativeTransferSolvers();

      /*!
	@brief Get sigma solver
      */
      virtual RefCountedPtr<SigmaSolver>& getSigmaSolver();
  
    protected:

      /*!
	@brief Realm
      */
      std::string m_realm;



      // Various parsing functions that are common to many derived classes
      virtual void parseVerbosity();
      virtual void parseSolverVerbosity();
      virtual void parseCFL();
      virtual void parseRelaxationTime();
      virtual void parseSourceGrowth();
      virtual void parseSourceTolerance();
      virtual void parseMinDt();
      virtual void parseMaxDt();
      virtual void parseFastPoisson();
      virtual void parseFastRadiativeTransfer();
      virtual void parseSourceComputation();

      /*!
	@brief Plasma phase
      */
      phase::which_phase m_phase;

      /*!
	@brief Time stepper class name
      */
      std::string m_className;

      /*!
	@brief Index space
      */
      RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

      /*!
	@brief Plasma kinetics
      */
      RefCountedPtr<CdrPlasmaPhysics> m_physics;

      /*!
	@brief CDR solvers
      */
      RefCountedPtr<CdrLayout<CdrSolver> > m_cdr;

      /*!
	@brief Radiative transfer solvers
      */
      RefCountedPtr<RtLayout<RtSolver>> m_rte;

      /*!
	@brief Poisson solver
      */
      RefCountedPtr<FieldSolver> m_fieldSolver;

      /*!
	@brief Surface charge solver
      */
      RefCountedPtr<SigmaSolver> m_sigma;

      /*!
	@brief Set the potential
      */
      std::function<Real(const Real a_time) > m_potential;
  
      /*!
	@brief Minimum allowed time step
      */
      Real m_min_dt;

      /*!
	@brief Maximum allowed time step
      */
      Real m_max_dt;

      /*!
	@brief CFL
      */
      Real m_cfl;

      /*!
	@brief Fudge factor for relaxation time
      */
      Real m_relax_time;

      /*!
	@brief Source growth
      */
      Real m_src_growth;

      /*!
	@brief TIme
      */
      Real m_time;

      /*!
	@brief Previous time step size
      */
      Real m_dt;

      /*!
	@brief Computed CFL time step
      */
      Real m_dt_cfl;

      /*!
	@brief Tolerance for source term growth
      */
      Real m_src_tolerance;

      /*!
	@brief Finest level on which to impose relaxation time constraints
      */
      int m_relax_level;

      /*!
	@brief Verbosity for solvers
      */
      int m_solver_verbosity;

      /*!
	@brief Fast RTE solver
      */
      int m_fast_rte;

      /*!
	@brief Fast RTE solver
      */
      int m_fast_poisson;

      /*!
	@brief Do source term time step restriction only for electrons
	@details If this one is true, computeDt passes information down the toolchain to only compute the source term
	restriction time for the first species defined in the plasma kinetics
      */
      bool m_src_elec_only;

      /*!
	@brief Which source term computation
      */
      SourceTermComputation m_whichSourceTermComputation;

      /*!
	@brief Upwind factor
      */
      int m_upwindFactor;


      /*!
	@brief Check if solver does subcyling. If it does, we need to tell AmrMesh to make some extra storage for us
      */
      bool m_subcycle;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
