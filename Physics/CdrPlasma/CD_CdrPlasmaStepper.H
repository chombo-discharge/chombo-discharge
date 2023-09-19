/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaStepper.H
  @brief  Declaration of a CdrPlasmaStepper, an abstract class for evolving cdr plasma equations in time.
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
#include <CD_SurfaceODESolver.H>
#include <CD_TimeStepper.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Abstract class for evolving the minimal plasma model in time.
      @details This class implements all the required coupling functions for more easily constructing algorithms. Subclasses must
      implement the advance method. This class handles all I/O.
    */
    class CdrPlasmaStepper : public TimeStepper
    {
    public:
      /*!
	@brief Empty constructor.
      */
      CdrPlasmaStepper();

      /*!
	@brief Full constructor. 
	@param[in] a_physics CDR plasma physics. 
      */
      CdrPlasmaStepper(RefCountedPtr<CdrPlasmaPhysics>& a_physics);

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaStepper();

      /*!
	@brief Instantiate solvers
      */
      virtual void
      setupSolvers() override;

      /*!
	@brief Allocate data for the time stepper and solvers. 
      */
      virtual void
      allocate() override;

      /*!
	@brief Fill all solvers with initial data
      */
      virtual void
      initialData() override;

      /*!
	@brief Post-initialize operations to be performed after filling solvers with initialData.
      */
      virtual void
      postInitialize() override;

      /*!
	@brief Post-initialize operations to be performed after filling solvers with data read from checkpoint files. 
      */
      virtual void
      postCheckpointSetup() override;

      /*!
	@brief Register realms to be used for the simulation.
      */
      virtual void
      registerRealms() override;

      /*!
	@brief Register operators to be used for the simulation
      */
      virtual void
      registerOperators() override;

      /*!
	@brief Parse runtime options. Subclasses must implement this.
      */
      virtual void
      parseRuntimeOptions() = 0;

      /*!
	@brief Compute some thing that go into plot files
      */
      virtual void
      prePlot() override;

#ifdef CH_USE_HDF5
      /*!  
	@brief Write checkpoint data to file.
	@param[inout] a_handle HDF5 fil
	@param[in]    a_lvl    Grid level
      */
      virtual void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Read checkpoint data from file. Used for restarts. 
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      virtual void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /*!
	@brief Get the number of plot variables for this time stepper. 
	@details This is necessary because Driver, not TimeStepper, is responsible for allocating the necessary memory. 
	@return Returns number of plot variables that will be written during writePlotData
      */
      virtual int
      getNumberOfPlotVariables() const override;

      /*!
	@brief Get plot variable names
      */
      virtual Vector<std::string>
      getPlotVariableNames() const override;

      /*!
	@brief Write plot data to output holder. 
	@param[inout] a_output      Output data holder.
	@param[inout] a_icomp       Starting component in a_output to begin at. 
	@param[in]    a_outputRealm Realm where a_output belongs
	@param[in]    a_level       Grid level
      */
      virtual void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string     a_outputRealm,
                    const int             a_level) const override;

      /*!
	@brief Compute the time step.
	@details Subclasses must implement this one. 
      */
      virtual Real
      computeDt() = 0;

      /*!
	@brief Advance method, advances equations.
	@param[in] a_dt Time step
	@returns Actual time step used (can be different from a_dt for adaptive methods)
      */
      virtual Real
      advance(const Real a_dt) = 0;

      /*!
	@brief Synchronize solver times
	@param[in] a_step Time step
	@param[in] a_time Time (in seconds)
	@param[in] a_dt   Time step that was used. 
      */
      virtual void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /*!
	@brief Print a step report. 
	@details This is called by Driver after time step. The routine can be used to display use information about the simulation progress. 
      */
      virtual void
      printStepReport() override;

      /*!
	@brief Pre-regrid function -- this is the one called by Driver. 
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      virtual void
      preRegrid(const int a_lmin, const int a_oldFinestLevel) override;

      /*!
	@brief Do a preRegrid operation for internal storage (if needed). 
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      virtual void
      preRegridInternals(const int a_lmin, const int a_oldFinestLevel);

      /*!
	@brief Regrid method.
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      virtual void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Perform post-regrid operations. 
	@details This includes all operations to be done AFTER interpolating data to new grids. 
      */
      virtual void
      postRegrid() override;

      /*!
	@brief Compute the current density
	@param[out] a_J Current density
      */
      virtual void
      computeJ(EBAMRCellData& a_J) const;

      /*!
	@brief Compute the cell-centered electric field on both phases
	@param[out] a_electricFieldCell Cell-centered electric field
	@param[in]  a_potential         Cell-centered electrostatic potential
      */
      virtual void
      computeElectricField(MFAMRCellData& a_electricFieldCell, const MFAMRCellData& a_potential) const;

      /*!
	@brief Compute the cell-centered electric field on a specific phase using whatever is available in m_fieldSolver
	@details This will use the potential from m_fieldSolver.
	@param[out] a_electricFieldCell Cell-centered electric field
	@param[in]  a_phase             Phase
      */
      virtual void
      computeElectricField(EBAMRCellData& a_electricFieldCell, const phase::which_phase a_phase) const;

      /*!
	@brief Compute the cell-centered electric field on a specific phase
	@param[out] a_electricFieldCell Cell-centered electric field
	@param[in]  a_phase             Phase
	@param[in]  a_potential         Electrostatic potential (cell-centered)
      */
      virtual void
      computeElectricField(EBAMRCellData&           a_E,
                           const phase::which_phase a_phase,
                           const MFAMRCellData&     a_potential) const;

      /*!
	@brief Compute the face-centered electric field. This will only do the normal component. 
	@param[out] a_electricFieldFace  Face-centered electric field
	@param[in]  a_phase              Phase
	@param[in]  a_electricFieldCell  Cell-centered electric field (all components). 
      */
      virtual void
      computeElectricField(EBAMRFluxData&           a_electricFieldFace,
                           const phase::which_phase a_phase,
                           const EBAMRCellData&     a_electricFieldCell) const;

      /*!
	@brief Compute the EB-centered electric field. 
	@param[out] a_electricFieldEB    EB-centered electric field. 
	@param[in]  a_phase              Phase
	@param[in]  a_electricFieldCell  Cell-centered electric field (all components). 
      */
      virtual void
      computeElectricField(EBAMRIVData&             a_electricFieldEB,
                           const phase::which_phase a_phase,
                           const EBAMRCellData&     a_electricFieldCell) const;

      /*!
	@brief Compute the maximum of the electric field
	@param[out] a_maxElectricField Maximum electric field magnitude evaluated over the input phase
	@param[in]  a_phase            Input phase
      */
      virtual void
      computeMaxElectricField(Real& a_maxElectricField, const phase::which_phase a_phase);

      /*!
	@brief Return time
	@return m_time.
      */
      virtual Real
      getTime() const;

      /*!
	@brief Check if RTE solvers are stationary
	@return True if the RTE solvers are stationary. 
      */
      virtual bool
      stationaryRTE();

      /*!
	@brief Allocate internal storage. Must be implemented by subclasses. 
      */
      virtual void
      allocateInternals() = 0;

      /*!
	@brief Deallocate internals. Must be implemented by subclasses. 
      */
      virtual void
      deallocateInternals() = 0;

      /*!
	@brief Parse TimeStepper class options. Must be implemented by subclasses.
      */
      virtual void
      parseOptions() = 0;

      /*!
	@brief Deallocation function. This will deallocate internal storage in the subclasses as well as the solvers. 
      */
      virtual void
      deallocate();

      /*!
	@brief Compute the centroid-centered space charge density by using the data inside the CDR solvers.
	@details This will call the other version with the storage holders being the space charge density data holder in m_fieldSolver and the 
	cell-centered densities in the CDR solvers. 
      */
      virtual void
      computeSpaceChargeDensity();

      /*!
	@brief Compute the centroid-centered space charge density.
	@details This will first compute the cell-centered space charge density and then interpolate that to the centroid. 
	@param[out] a_rho          Centroid centered space charge density.
	@param[in]  a_cdrDensities Cell-centered plasma species densities. Must be in the same order as the densities in m_cdr.
      */
      virtual void
      computeSpaceChargeDensity(MFAMRCellData& a_rho, const Vector<EBAMRCellData*>& a_cdrDensities);

      /*!
	@brief Compute the cell-centered conductivity.
	@param[out] a_cellConductivity Cell-centered conductivity. 
      */
      virtual void
      computeCellConductivity(EBAMRCellData& a_cellConductivity) const;

      /*!
	@brief Compute the face-centered conductivity from a cell-centered conductivity.
	@param[out] a_conductivityFace Face-centered conductivity.
	@param[out] a_conductivityEB   EB-centered conductivity. 
	@param[in]  a_conductivityCell Cell-centered conductivity. 
      */
      virtual void
      computeFaceConductivity(EBAMRFluxData&       a_conductivityFace,
                              EBAMRIVData&         a_conductivityEB,
                              const EBAMRCellData& a_conductivityCell) const;

      /*!
	@brief Set up a semi-implicit Poisson solver.
	@details This will compute the required conductivities first, using data already present in the CDR solvers. After that, it
	will set up the Poisson equation with spatially varying conductivities.
	@param[in] a_dt Time step. 
      */
      virtual void
      setupSemiImplicitPoisson(const Real a_dt);

      /*!
	@brief Set up a semi-implicit Poisson solver.
	@details This is the version that uses pre-existing conductivities in a_conductivity. The factor dictates the multiplication factor f in the 
	expression div((epsilonr + f*sigma) grad(phi)) where sigma is the conductivity and f is a scaling factor. 
	@param[in] a_conductivityFace Face-centered conductivities (the Poisson solver will interpolate them to face centroids in the multigrid solve).
	@param[in] a_conductivityEB   EB-centered conductivities.
	@param[in] a_factor Factor to multiply conductivities with (they could be pre-scaled by dt/eps0 in which case a_factor = 1).
      */
      virtual void
      setupSemiImplicitPoisson(const EBAMRFluxData& a_conductivityFace,
                               const EBAMRIVData&   a_conductivityEB,
                               const Real           a_factor);

      /*!
	@brief Solve the electrostatic Poisson equation. 
	@details This will first compute the centroid-centered space charge density (using the cell-centered data in the CDR solvers) and then
	ask the field solver to solve.
      */
      virtual bool
      solvePoisson();

      /*!
	@brief General electrostatic Poisson solver routine. 
	@details This will compute the space charge density onto a_rho and then solve for the potential. 
	@param[inout] a_potential    Electrostatic potential
	@param[out]   a_rho          Centroid-centered space charge density. 
	@param[in]    a_cdrDensities Plasma species densities. 
	@param[in]    a_sigma        Surface charge density on dielectric EBs.
      */
      virtual bool
      solvePoisson(MFAMRCellData&               a_potential,
                   MFAMRCellData&               a_rho,
                   const Vector<EBAMRCellData*> a_cdrDensities,
                   const EBAMRIVData&           a_sigma);

      /*!
	@brief Advance the reaction network. This will compute the electric field (on the appropriate phase) and call the other version.
	@param[in] a_time Time
	@param[in] a_dt   Time step
      */
      virtual void
      advanceReactionNetwork(const Real a_time, const Real a_dt);

      /*!
	@brief Compute reaction network source terms. 
	@details This will call the other version but will first compute the gradients of the plasma densities. 
	@param[out] a_cdrSources   Source terms for the CDR solvers. 
	@param[out] a_rteSources   Source terms for the radiative transfer solvers. 
	@param[in]  a_cdrDensities Cell-centered densities for the plasma species (i.e., CDR solvers).
	@param[in]  a_rteDensities Cell-centered densities for the radiative transfer equations (i.e., RTE solvers).
	@param[in]  a_E            Cell-centered electric field. 
	@param[in]  a_time         Time
	@param[in]  a_dt           Time step.
      */
      virtual void
      advanceReactionNetwork(Vector<EBAMRCellData*>&       a_cdrSources,
                             Vector<EBAMRCellData*>&       a_rteSources,
                             const Vector<EBAMRCellData*>& a_cdrDensities,
                             const Vector<EBAMRCellData*>& a_rteDensities,
                             const EBAMRCellData&          a_E,
                             const Real&                   a_time,
                             const Real&                   a_dt);

      /*!
	@brief Compute reaction network source terms. This is the AMR version that does all levels. 
	@param[out] a_cdrSources   Source terms for the CDR solvers. 
	@param[out] a_rteSources   Source terms for the radiative transfer solvers. 
	@param[in]  a_cdrDensities Cell-centered densities for the plasma species (i.e., CDR solvers).
	@param[in]  a_cdrGradients Cell-centered gradients for the plasma species (i.e., CDR solvers).
	@param[in]  a_rteDensities Cell-centered densities for the radiative transfer equations (i.e., RTE solvers).
	@param[in]  a_E            Cell-centered electric field. 
	@param[in]  a_time         Time
	@param[in]  a_dt           Time step.
      */
      virtual void
      advanceReactionNetwork(Vector<EBAMRCellData*>&       a_cdrSources,
                             Vector<EBAMRCellData*>&       a_rteSources,
                             const Vector<EBAMRCellData*>& a_cdrDensities,
                             const Vector<EBAMRCellData*>& a_cdrGradients,
                             const Vector<EBAMRCellData*>& a_rteDensities,
                             const EBAMRCellData&          a_E,
                             const Real&                   a_time,
                             const Real&                   a_dt);

      /*!
	@brief Compute reaction network source terms. This is the level version. 
	@param[out] a_cdrSources   Source terms for the CDR solvers. 
	@param[out] a_rteSources   Source terms for the radiative transfer solvers. 
	@param[in]  a_cdrDensities Cell-centered densities for the plasma species (i.e., CDR solvers).
	@param[in]  a_cdrGradients Cell-centered gradients for the plasma species (i.e., CDR solvers).
	@param[in]  a_rteDensities Cell-centered densities for the radiative transfer equations (i.e., RTE solvers).
	@param[in]  a_E            Cell-centered electric field. 
	@param[in]  a_time         Time
	@param[in]  a_dt           Time step.
	@param[in]  a_level        Grid level. 
      */
      virtual void
      advanceReactionNetwork(Vector<LevelData<EBCellFAB>*>&       a_cdrSources,
                             Vector<LevelData<EBCellFAB>*>&       a_rteSources,
                             const Vector<LevelData<EBCellFAB>*>& a_cdrDensities,
                             const Vector<LevelData<EBCellFAB>*>& a_cdrGradients,
                             const Vector<LevelData<EBCellFAB>*>& a_rteDensities,
                             const LevelData<EBCellFAB>&          a_E,
                             const Real&                          a_time,
                             const Real&                          a_dt,
                             const int                            a_level);
      /*!
	@brief Compute reaction network sources. This is the regular version which does not do the cut-cells.
	@param[out] a_cdrSources   Source terms for the CDR solvers. 
	@param[out] a_rteSources   Source terms for the radiative transfer solvers. 
	@param[in]  a_cdrDensities Cell-centered densities for the plasma species (i.e., CDR solvers).
	@param[in]  a_cdrGradients Cell-centered gradients for the plasma species (i.e., CDR solvers).
	@param[in]  a_rteDensities Cell-centered densities for the radiative transfer equations (i.e., RTE solvers).
	@param[in]  a_E            Cell-centered electric field. 
	@param[in]  a_time         Time
	@param[in]  a_dt           Time step.
	@param[in]  a_dx           Grid resolution.
	@param[in]  a_cellBox      Cell-centered computational region. 
      */
      virtual void
      advanceReactionNetworkRegularCells(Vector<FArrayBox*>&       a_cdrSources,
                                         Vector<FArrayBox*>&       a_rteSources,
                                         const Vector<FArrayBox*>& a_cdrDensities,
                                         const Vector<FArrayBox*>& a_cdrGradients,
                                         const Vector<FArrayBox*>& a_rteDensities,
                                         const FArrayBox&          a_E,
                                         const Real&               a_time,
                                         const Real&               a_dt,
                                         const Real&               a_dx,
                                         const Box&                a_cellBox);

      /*!
	@brief Compute reaction network sources. This is the irregular version which only does the cut-cells. 
	@param[out] a_cdrSources     Source terms for the CDR solvers. 
	@param[out] a_rteSources     Source terms for the radiative transfer solvers. 
	@param[in]  a_cdrDensities   Cell-centered densities for the plasma species (i.e., CDR solvers).
	@param[in]  a_cdrGradients   Cell-centered gradients for the plasma species (i.e., CDR solvers).
	@param[in]  a_cdrVelocities  Cell-centered velocities for the CDR solvers. 
	@param[in]  a_rteDensities   Cell-centered densities for the radiative transfer equations (i.e., RTE solvers).
	@param[in]  a_interpStencils Interpolation stencils for putting cell-centered data on the centroid. 
	@param[in]  a_E              Cell-centered electric field. 
	@param[in]  a_time           Time
	@param[in]  a_dt             Time step.
	@param[in]  a_dx             Grid resolution.
	@param[in]  a_cellBox        Cell-centered computational region. 
	@param[in]  a_lvl            Grid level
	@param[in]  a_dit            Grid patch index
	@note This will call one of the other versions. 
      */
      virtual void
      advanceReactionNetworkIrreg(Vector<EBCellFAB*>&          a_cdrSources,
                                  Vector<EBCellFAB*>&          a_rteSources,
                                  const Vector<EBCellFAB*>&    a_cdrDensities,
                                  const Vector<EBCellFAB*>&    a_cdrGradients,
                                  const Vector<EBCellFAB*>&    a_cdrVelocities,
                                  const Vector<EBCellFAB*>&    a_rteDensities,
                                  const BaseIVFAB<VoFStencil>& a_interpStencils,
                                  const EBCellFAB&             a_E,
                                  const Real&                  a_time,
                                  const Real&                  a_dt,
                                  const Real&                  a_dx,
                                  const Box&                   a_cellBox,
                                  const int                    a_lvl,
                                  const DataIndex&             a_dit);

      /*!
	@brief Compute reaction network sources. This version interpolates everything to the centroid and runs with that data. 
	@param[out] a_cdrSources     Source terms for the CDR solvers. 
	@param[out] a_rteSources     Source terms for the radiative transfer solvers. 
	@param[in]  a_cdrDensities   Cell-centered densities for the plasma species (i.e., CDR solvers).
	@param[in]  a_cdrGradients   Cell-centered gradients for the plasma species (i.e., CDR solvers).
	@param[in]  a_cdrVelocities  Cell-centered velocities for the CDR solvers. 
	@param[in]  a_rteDensities   Cell-centered densities for the radiative transfer equations (i.e., RTE solvers).
	@param[in]  a_interpStencils Interpolation stencils for putting cell-centered data on the centroid. 
	@param[in]  a_E              Cell-centered electric field. 
	@param[in]  a_time           Time
	@param[in]  a_dt             Time step.
	@param[in]  a_dx             Grid resolution.
	@param[in]  a_cellBox        Cell-centered computational region. 
	@param[in]  a_lvl            Grid level
	@param[in]  a_dit            Grid patch index
      */
      virtual void
      advanceReactionNetworkIrregInterp(Vector<EBCellFAB*>&          a_cdrSources,
                                        Vector<EBCellFAB*>&          a_rteSources,
                                        const Vector<EBCellFAB*>&    a_cdrDensities,
                                        const Vector<EBCellFAB*>&    a_cdrGradients,
                                        const Vector<EBCellFAB*>&    a_cdrVelocities,
                                        const Vector<EBCellFAB*>&    a_rteDensities,
                                        const BaseIVFAB<VoFStencil>& a_interpStencils,
                                        const EBCellFAB&             a_E,
                                        const Real&                  a_time,
                                        const Real&                  a_dt,
                                        const Real&                  a_dx,
                                        const Box&                   a_box,
                                        const int                    a_lvl,
                                        const DataIndex&             a_dit);
      /*!
	@brief Compute reaction network sources. This version assumes that some of the inputs are on the centroid. 
	@param[out] a_cdrSources     Source terms for the CDR solvers. 
	@param[out] a_rteSources     Source terms for the radiative transfer solvers. 
	@param[in]  a_cdrDensities   Centroid-centered densities for the plasma species (i.e., CDR solvers).
	@param[in]  a_cdrGradients   Cell-centered gradients for the plasma species (i.e., CDR solvers).
	@param[in]  a_rteDensities   Centroid-centered densities for the radiative transfer equations (i.e., RTE solvers).
	@param[in]  a_interpStencils Interpolation stencils for putting cell-centered data on the centroid. 
	@param[in]  a_E              Cell-centered electric field. 
	@param[in]  a_time           Time
	@param[in]  a_dt             Time step.
	@param[in]  a_dx             Grid resolution.
	@param[in]  a_cellBox        Cell-centered computational region. 
	@param[in]  a_lvl            Grid level
	@param[in]  a_dit            Grid patch index      
      */
      virtual void
      advanceReactionNetworkIrregKappa(Vector<EBCellFAB*>&          a_cdrSources,
                                       Vector<EBCellFAB*>&          a_rteSources,
                                       const Vector<EBCellFAB*>&    a_cdrDensities,
                                       const Vector<EBCellFAB*>&    a_cdrGradients,
                                       const Vector<EBCellFAB*>&    a_rteDensities,
                                       const BaseIVFAB<VoFStencil>& a_interpStencils,
                                       const EBCellFAB&             a_E,
                                       const Real&                  a_time,
                                       const Real&                  a_dt,
                                       const Real&                  a_dx,
                                       const Box&                   a_box,
                                       const int                    a_lvl,
                                       const DataIndex&             a_dit);

      /*!
	@brief Compute reaction network sources. 
	@details This version is called with the Villa et. al upwind correction. It can also be used with regularly interpolated data, in which 
	case it will turn off reactions in EB cells where the mass flows into the domain. This function assumes that all input is on the cell center. 
	@param[out] a_cdrSources     Source terms for the CDR solvers. 
	@param[out] a_rteSources     Source terms for the radiative transfer solvers. 
	@param[in]  a_cdrDensities   Cell-centered densities for the plasma species (i.e., CDR solvers).
	@param[in]  a_cdrGradients   Cell-centered gradients for the plasma species (i.e., CDR solvers).
	@param[in]  a_cdrVelocities  Cell-centered velocities for the CDR solvers. 
	@param[in]  a_rteDensities   Cell-centered densities for the radiative transfer equations (i.e., RTE solvers).
	@param[in]  a_interpStencils Interpolation stencils for putting cell-centered data on the centroid. 
	@param[in]  a_E              Cell-centered electric field. 
	@param[in]  a_time           Time
	@param[in]  a_dt             Time step.
	@param[in]  a_dx             Grid resolution.
	@param[in]  a_cellBox        Cell-centered computational region. 
	@param[in]  a_lvl            Grid level
	@param[in]  a_dit            Grid patch index
      */
      virtual void
      advanceReactionNetworkIrregUpwind(Vector<EBCellFAB*>&          a_cdrSources,
                                        Vector<EBCellFAB*>&          a_rteSources,
                                        const Vector<EBCellFAB*>&    a_cdrDensities,
                                        const Vector<EBCellFAB*>&    a_cdrGradients,
                                        const Vector<EBCellFAB*>&    a_cdrVelocities,
                                        const Vector<EBCellFAB*>&    a_rteDensities,
                                        const BaseIVFAB<VoFStencil>& a_interpStencils,
                                        const EBCellFAB&             a_E,
                                        const Real&                  a_time,
                                        const Real&                  a_dt,
                                        const Real&                  a_dx,
                                        const Box&                   a_box,
                                        const int                    a_lvl,
                                        const DataIndex&             a_dit);

      /*!
	@brief Compute diffusion things using whatever is available in the solvers. 
	@details This will compute the electric field on the cell center and EB, and then call the other version. 
      */
      virtual void
      computeCdrDiffusion();

      /*!
	@brief Compute diffusion things using whatever is available in the solvers. This version uses the input electric fields.
	@param[in] a_electricFieldCell Electric field on cell centers. 
	@param[in] a_electricFieldEB   Electric field on EB centroids
      */
      virtual void
      computeCdrDiffusion(const EBAMRCellData& a_electricFieldCell, const EBAMRIVData& a_electricFieldEB);

      /*!
	@brief Compute CDR diffusion coefficients on cell centers. This will call the level version. 
	@param[out] a_cdrDcoCell        CDR diffusion coefficients on cell centers. 
	@param[in]  a_cdrDensities      Cell-centered CDR densities
	@param[in]  a_electricFieldCell Cell-centered electric field
	@param[in]  a_time              Time
      */
      virtual void
      computeCdrDiffusionCell(Vector<EBAMRCellData>&        a_cdrDcoCell,
                              const Vector<EBAMRCellData*>& a_cdrDensities,
                              const EBAMRCellData&          a_electricFieldCell,
                              const Real&                   a_time);

      /*!
	@brief Compute diffusion coefficients. This is the level version and it will iterate through grid patches. 
	@param[out] a_cdrDcoCell        CDR diffusion coefficients on cell centers. 
	@param[in]  a_cdrDensities      Cell-centered CDR densities
	@param[in]  a_electricFieldCell Cell-centered electric field
	@param[in]  a_lvl               Grid level
	@param[in]  a_time              Time
      */
      virtual void
      computeCdrDiffusionCell(Vector<LevelData<EBCellFAB>*>&       a_cdrDcoCell,
                              const Vector<LevelData<EBCellFAB>*>& a_cdrDensities,
                              const LevelData<EBCellFAB>&          a_electricFieldCell,
                              const int                            a_lvl,
                              const Real&                          a_time);

      /*!
	@brief Compute diffusion coefficients on cell centers. This does all the regular cells. 
	@param[out] a_cdrDcoCell        CDR diffusion coefficients on cell centers. 
	@param[in]  a_cdrDensities      Cell-centered CDR densities
	@param[in]  a_electricFieldCell Cell-centered electric field
	@param[in]  a_cellBox           Kernel region
	@param[in]  a_dx                Grid resolution
	@param[in]  a_time              Time
      */
      virtual void
      computeCdrDiffusionCellRegular(Vector<FArrayBox*>&       a_cdrDcoCell,
                                     const Vector<FArrayBox*>& a_cdrDensities,
                                     const FArrayBox&          a_electricFieldCell,
                                     const Box                 a_cellBox,
                                     const Real                a_dx,
                                     const Real                a_time);

      /*!
	@brief Compute diffusion coefficients on cell centers. This version does the irregular cells. 
	@param[out] a_cdrDcoCell        CDR diffusion coefficients on cell centers. 
	@param[in]  a_cdrDensities      Cell-centered CDR densities
	@param[in]  a_electricFieldCell Cell-centered electric field
	@param[in]  a_dx                Grid resolution
	@param[in]  a_time              Time
	@param[in]  a_lvl               Grid level
	@param[in]  a_dit               Grid patch index
      */
      virtual void
      computeCdrDiffusionCellIrregular(Vector<EBCellFAB*>&       a_cdrDcoCell,
                                       const Vector<EBCellFAB*>& a_cdrDensities,
                                       const EBCellFAB&          a_electricFieldCell,
                                       const Real                a_dx,
                                       const Real&               a_time,
                                       const int                 a_lvl,
                                       const DataIndex&          a_dit);

      /*!
	@brief Compute diffusion coefficients on face centers. 
	@details This will first compute the diffusion coefficients on cell centers, and then interpolate the result to face centers.
	@param[out] a_cdrDcoCell        CDR diffusion coefficients on cell centers. 
	@param[in]  a_cdrDensities      Cell-centered CDR densities
	@param[in]  a_electricFieldCell Cell-centered electric field
	@param[in]  a_time              Time
      */
      virtual void
      computeCdrDiffusionFace(Vector<EBAMRFluxData*>&       a_cdrDcoFace,
                              const Vector<EBAMRCellData*>& a_cdrDensities,
                              const EBAMRCellData&          a_electricFieldCell,
                              const Real&                   a_time);

      /*!
	@brief Compute diffusion coefficients on the EB. This is the AMR version -- it will call the level version. 
	@param[out] a_cdrDcoEB        EB-centered_CDR diffusion coefficients
	@param[in]  a_cdrDensitiesEB  EB-centered CDR densities
	@param[in]  a_electricFieldEB EB-centered electric field
	@param[in]  a_time            Time
      */
      virtual void
      computeCdrDiffusionEb(Vector<EBAMRIVData*>&       a_cdrDcoEB,
                            const Vector<EBAMRIVData*>& a_cdrDensitiesEB,
                            const EBAMRIVData&          a_electricFieldEB,
                            const Real&                 a_time);

      /*!
	@brief Compute diffusion coefficients on the EB. This is a level version but it will go through all patches. 
	@param[out] a_cdrDcoEB        EB-centered_CDR diffusion coefficients
	@param[in]  a_cdrDensitiesEB  EB-centered CDR densities
	@param[in]  a_electricFieldEB EB-centered electric field
	@param[in]  a_time            Time
	@param[in]  a_lvl             Grid level. 
      */
      virtual void
      computeCdrDiffusionEb(Vector<LevelData<BaseIVFAB<Real>>*>&       a_cdrDcoEB,
                            const Vector<LevelData<BaseIVFAB<Real>>*>& a_cdrDensitiesEB,
                            const LevelData<BaseIVFAB<Real>>&          a_electricFieldEB,
                            const Real&                                a_time,
                            const int                                  a_lvl);

      /*!
	@brief Compute the CDR drift velocities.
	@details This will compute the cell-centered electric field and then call the other version.
      */
      virtual void
      computeCdrDriftVelocities();

      /*!
	@brief Compute the cell-centered CDR velocities. This will call the level version.
	@param[out] a_cdrVelocities Cell-centered CDR drift velocities
	@param[out] a_cdrDensities  Cell-centered CDR densities
	@param[in]  a_electricField The cell-centered electric field
	@param[in]  a_time          Time
      */
      virtual void
      computeCdrDriftVelocities(Vector<EBAMRCellData*>&       a_velocities,
                                const Vector<EBAMRCellData*>& a_cdrDensities,
                                const EBAMRCellData&          a_E,
                                const Real&                   a_time);

      /*!
	@brief Compute the cell-centered CDR velocities. 
	@details This is the level version -- it will call the patch version for each patch in the grid. 
	@param[out] a_cdrVelocities Cell-centered CDR drift velocities
	@param[out] a_cdrDensities  Cell-centered CDR densities
	@param[in]  a_electricField The cell-centered electric field
	@param[in]  a_lvl           Grid level
	@param[in]  a_time          Time
      */
      virtual void
      computeCdrDriftVelocities(Vector<LevelData<EBCellFAB>*>&       a_cdrVelocities,
                                const Vector<LevelData<EBCellFAB>*>& a_cdrDensities,
                                const LevelData<EBCellFAB>&          a_electricField,
                                const int                            a_lvl,
                                const Real&                          a_time);

      /*!
	@brief Compute the cell-centered CDR velocities. 
	@details This is patch version for regular grid cells. 
	@param[out] a_cdrVelocities Cell-centered CDR drift velocities
	@param[out] a_cdrDensities  Cell-centered CDR densities
	@param[in]  a_electricField The cell-centered electric field
	@param[in]  a_cellBox       Computational region. 
	@param[in]  a_time          Time
	@param[in]  a_dx            Grid resolution
      */
      virtual void
      computeCdrDriftVelocitiesRegular(Vector<FArrayBox*>&       a_cdrVelocities,
                                       const Vector<FArrayBox*>& a_cdrDensities,
                                       const FArrayBox&          a_electricField,
                                       const Box&                a_cellBox,
                                       const Real&               a_time,
                                       const Real&               a_dx);

      /*!
	@brief Compute the cell-centered CDR velocities. 
	@details This is patch version for irregular grid cells. 
	@param[out] a_cdrVelocities Cell-centered CDR drift velocities
	@param[out] a_cdrDensities  Cell-centered CDR densities
	@param[in]  a_electricField The cell-centered electric field
	@param[in]  a_time          Time
	@param[in]  a_dx            Grid resolution
	@param[in]  a_lvl           Grid level
	@param[in]  a_dit           Grid index
      */
      virtual void
      computeCdrDriftVelocitiesIrregular(Vector<EBCellFAB*>&       a_velocities,
                                         const Vector<EBCellFAB*>& a_cdrDensities,
                                         const EBCellFAB&          a_electricField,
                                         const Real&               a_time,
                                         const Real&               a_dx,
                                         const int                 a_lvl,
                                         const DataIndex&          a_dit);

      /*!
	@brief Compute CDR fluxes on the EB. This is the AMR version -- it will call the level version. 
	@param[out] a_cdrFluxesEB         CDR fluxes to be put in CDR solvers as a BC on the EBs.
	@param[in]  a_extrapCdrFluxes     CDR fluxes extrapolated from the inside. Can be used to construct inflow/outflow BCs.
	@param[in]  a_extrapCdrDensities  CDR densities extrapolated from the inside. 
	@param[in]  a_extrapCdrVelocities CDR velocities on the EB, extrapolated from the inside. 
	@param[in]  a_extrapCdrGradients  CDR gradients on the EB, extrapolated from the inside. 
	@param[in]  a_electricField       Electric field on the EB. 
	@param[in]  a_time                Time. 
      */
      virtual void
      computeCdrFluxes(Vector<EBAMRIVData*>&       a_cdrFluxesEB,
                       const Vector<EBAMRIVData*>& a_extrapCdrFluxes,
                       const Vector<EBAMRIVData*>& a_extrapCdrDensities,
                       const Vector<EBAMRIVData*>& a_extrapCdrVelocities,
                       const Vector<EBAMRIVData*>& a_extrapCdrGradients,
                       const Vector<EBAMRIVData*>& a_extrapRteFluxes,
                       const EBAMRIVData&          a_electricField,
                       const Real&                 a_time);

      /*!
	@brief Compute CDR fluxes on the EB. This is the level version which will fill EB fluxes on both electrode and dielectric cells. 
	@param[out] a_cdrFluxesEB         CDR fluxes to be put in CDR solvers as a BC on the EBs.
	@param[in]  a_extrapCdrFluxes     CDR fluxes extrapolated from the inside. Can be used to construct inflow/outflow BCs.
	@param[in]  a_extrapCdrDensities  CDR densities extrapolated from the inside. 
	@param[in]  a_extrapCdrVelocities CDR velocities on the EB, extrapolated from the inside. 
	@param[in]  a_extrapCdrGradients  CDR gradients on the EB, extrapolated from the inside. 
	@param[in]  a_electricField       Electric field on the EB. 
	@param[in]  a_time                Time. 
	@param[in]  a_lvl                 Grid level. 
      */
      virtual void
      computeCdrFluxes(Vector<LevelData<BaseIVFAB<Real>>*>&       a_cdrFluxesEB,
                       const Vector<LevelData<BaseIVFAB<Real>>*>& a_extrapCdrFluxes,
                       const Vector<LevelData<BaseIVFAB<Real>>*>& a_extrapCdrDensities,
                       const Vector<LevelData<BaseIVFAB<Real>>*>& a_extrapCdrVelocities,
                       const Vector<LevelData<BaseIVFAB<Real>>*>& a_extrapCdrGradients,
                       const Vector<LevelData<BaseIVFAB<Real>>*>& a_extrapRteFluxes,
                       const LevelData<BaseIVFAB<Real>>&          a_electricField,
                       const Real&                                a_time,
                       const int                                  a_lvl);

      /*!
	@brief Compute CDR fluxes on domain faces for usage in boundary conditions. This is the AMR version. 
	@param[out] a_cdrFluxes           CDR domain bounadry fluxes for solvers.
	@param[in]  a_extrapCdrFluxes     CDR fluxes extrapolated to the domain faces. 
	@param[in]  a_extrapCdrDensities  CDR densities extrapolated to the domain faces. 
	@param[in]  a_extrapCdrVelocities CDR velocities extrapolated to the domain faces. 
	@param[in]  a_extrapCdrGradients  CDR gradients extrapolated to the domain faces. 
	@param[in]  a_extrapRteFluxes     RTE fluxes on the domain faces. 
	@param[in]  a_electricField       Electric field on the domain faces. 
	@param[in]  a_time                Time.
      */
      virtual void
      computeCdrDomainFluxes(Vector<EBAMRIFData*>&       a_cdrFluxes,
                             const Vector<EBAMRIFData*>& a_extrapCdrFluxes,
                             const Vector<EBAMRIFData*>& a_extrapCdrDensities,
                             const Vector<EBAMRIFData*>& a_extrapCdrVelocities,
                             const Vector<EBAMRIFData*>& a_extrapCdrGradients,
                             const Vector<EBAMRIFData*>& a_extrapRteFluxes,
                             const EBAMRIFData&          a_field,
                             const Real&                 a_time);

      /*!
	@brief Compute CDR fluxes on domain faces for usage in boundary conditions. This is the level version. 
	@details This will iterate over the patches on each level and compute the domain fluxes. 
	@param[out] a_cdrFluxes           CDR domain bounadry fluxes for solvers.
	@param[in]  a_extrapCdrFluxes     CDR fluxes extrapolated to the domain faces. 
	@param[in]  a_extrapCdrDensities  CDR densities extrapolated to the domain faces. 
	@param[in]  a_extrapCdrVelocities CDR velocities extrapolated to the domain faces. 
	@param[in]  a_extrapCdrGradients  CDR gradients extrapolated to the domain faces. 
	@param[in]  a_extrapRteFluxes     RTE fluxes on the domain faces. 
	@param[in]  a_electricField       Electric field on the domain faces. 
	@param[in]  a_time                Time.
	@param[in]  a_lvl                 Grid level
      */
      virtual void
      computeCdrDomainFluxes(Vector<LevelData<DomainFluxIFFAB>*>        a_cdrFluxes,
                             const Vector<LevelData<DomainFluxIFFAB>*>& a_extrapCdrFluxes,
                             const Vector<LevelData<DomainFluxIFFAB>*>& a_extrapCdrDensities,
                             const Vector<LevelData<DomainFluxIFFAB>*>& a_extrapCdrVelocities,
                             const Vector<LevelData<DomainFluxIFFAB>*>& a_extrapCdrGradients,
                             const Vector<LevelData<DomainFluxIFFAB>*>& a_extrapRteFluxes,
                             const LevelData<DomainFluxIFFAB>&          a_electricField,
                             const Real&                                a_time,
                             const int                                  a_lvl);

      /*!
	@brief Extrapolate cell-centered fluxes to the EB centroid.
	@param[out] a_extrapCdrFluxesEB CDR fluxes extrapolated to the EB.
	@param[in]  a_cdrDensities      CDR densities on the cell center
	@param[in]  a_cdrVelocities     CDR velocities on the cell center
	@param[in]  a_phase             Phase
      */
      virtual void
      computeExtrapolatedFluxes(Vector<EBAMRIVData*>&        a_extrapCdrFluxesEB,
                                const Vector<EBAMRCellData*> a_cdrDensities,
                                const Vector<EBAMRCellData*> a_cdrVelocities,
                                const phase::which_phase     a_phase);

      /*!
	@brief Extrapolate drift fluxes to domain faces. 
	@details Note that this projects the flux on the domain face vectors. 
	@param[out] a_cdrDomainFluxes CDR fluxes extrapolated and projected onto the domain face. 
	@param[in]  a_cdrDensities    Cell-centered CDR solver densities. 
	@param[in]  a_cdrVelocities   Cell-centered CDR solver velocities. 
	@param[in]  a_phase           Phase
      */
      virtual void
      computeExtrapolatedDomainFluxes(Vector<EBAMRIFData*>&        a_cdrDomainFluxes,
                                      const Vector<EBAMRCellData*> a_cdrDensities,
                                      const Vector<EBAMRCellData*> a_cdrVelocities,
                                      const phase::which_phase     a_phase);

      /*!
	@brief Extrapolated cell-centered velocities to the EB. 
	@details This only computes the normal velocity on the EB, i.e. the output is the extrapolated velocity projected along the EB normal. 
	@param[out] a_cdrVelocitiesEB   EB-centered velocities
	@param[in]  a_cdrVelocitiesCell Cell-centered velocities
	@param[in]  a_phase             Phase
      */
      virtual void
      computeExtrapolatedVelocities(Vector<EBAMRIVData*>&        a_cdrVelocitiesEB,
                                    const Vector<EBAMRCellData*> a_cdrVelocitiesCell,
                                    const phase::which_phase     a_phase);

      /*!
	@brief Deallocate internal solver storages
      */
      virtual void
      deallocateSolverInternals();

      /*!
	@brief Extrapolated cell-centered data to the EB. This is the AMR version for all species. 
	@param[out] a_ebData   EB-centered data
	@param[in]  a_phase    Phase 
	@param[in]  a_cellData Cell-centered data
      */
      virtual void
      extrapolateToEb(Vector<EBAMRIVData*>&         a_ebData,
                      const phase::which_phase      a_phase,
                      const Vector<EBAMRCellData*>& a_cellData);

      /*!
	@brief Extrapolated cell-centered data to the EB. This is the AMR version.
	@param[out] a_ebData   EB-centered data
	@param[in]  a_phase    Phase 
	@param[in]  a_cellData Cell-centered data
      */
      virtual void
      extrapolateToEb(EBAMRIVData& a_ebData, const phase::which_phase a_phase, const EBAMRCellData& a_cellData);

      /*!
	@brief Extrapolated cell-centered data to the EB. This is the level version. 
	@param[out] a_ebData   EB-centered data
	@param[in]  a_phase    Phase 
	@param[in]  a_cellData Cell-centered data
	@param[in]  a_lvl      Grid level
      */
      virtual void
      extrapolateToEb(LevelData<BaseIVFAB<Real>>& a_ebData,
                      const phase::which_phase    a_phase,
                      const LevelData<EBCellFAB>& a_cellData,
                      const int                   a_lvl);

      /*!
	@brief Extrapolate cell-centered data to domain faces. This is the AMR version that does all species. 
	@param[out] a_domainData Domain-centered data
	@param[in]  a_phase      Phase
	@param[in]  a_cellData   Cell-centered data. 
      */
      virtual void
      extrapolateToDomainFaces(Vector<EBAMRIFData*>&         a_domainData,
                               const phase::which_phase      a_phase,
                               const Vector<EBAMRCellData*>& a_cellData);

      /*!
	@brief Extrapolate cell-centered data to domain faces. This is the AMR version.
	@param[out] a_domainData Domain-centered data
	@param[in]  a_phase      Phase
	@param[in]  a_cellData   Cell-centered data. 
      */
      virtual void
      extrapolateToDomainFaces(EBAMRIFData&             a_domainData,
                               const phase::which_phase a_phase,
                               const EBAMRCellData&     a_cellData);

      /*!
	@brief Extrapolate cell-centered data to domain faces. This is the AMR version.
	@param[out] a_domainData Domain-centered data
	@param[in]  a_phase      Phase
	@param[in]  a_cellData   Cell-centered data. 
	@param[in]  a_lvl        Grid level
      */
      virtual void
      extrapolateToDomainFaces(LevelData<DomainFluxIFFAB>& a_domainData,
                               const phase::which_phase    a_phase,
                               const LevelData<EBCellFAB>& a_cellData,
                               const int                   a_lvl);

      /*!
	@brief Extrapolate vector data to domain faces. This is the AMR version that does all species. 
	@details This will extrapolate the vector data to domain faces and project it along the domain normal. 
	@param[out] a_domainData Domain-centered data. Must have one component. 
	@param[in]  a_phase      Phase
	@param[in]  a_cellData   Cell-centered data. Must have SpaceDim components. 
      */
      virtual void
      extrapolateVectorToDomainFaces(Vector<EBAMRIFData*>&         a_domainData,
                                     const phase::which_phase      a_phase,
                                     const Vector<EBAMRCellData*>& a_cellData);

      /*!
	@brief Extrapolate vector data to domain faces. This is the AMR version.
	@details This will extrapolate the vector data to domain faces and project it along the domain normal. 
	@param[out] a_domainData Domain-centered data. Must have one component. 
	@param[in]  a_phase      Phase
	@param[in]  a_cellData   Cell-centered data. Must have SpaceDim components. 
      */
      virtual void
      extrapolateVectorToDomainFaces(EBAMRIFData&             a_domainData,
                                     const phase::which_phase a_phase,
                                     const EBAMRCellData&     a_cellData);

      /*!
	@brief Extrapolate velocities data to domain faces. This is the AMR version that does all species. 
	@details This will extrapolate the velocities data to domain faces.
	@param[out] a_domainVelocities Domain-centered velocities. Must have SpaceDim components. 
	@param[in]  a_phase            Phase
	@param[in]  a_cellData         Cell-centered velocities. Must have SpaceDim components. 
      */
      virtual void
      extrapolateVelocitiesToDomainFaces(Vector<EBAMRIFData*>&         a_domainVelocities,
                                         const phase::which_phase      a_phase,
                                         const Vector<EBAMRCellData*>& a_cellVelocities);

      /*!
	@brief Get maximum density in the CDR species. This will fetch both the value and the solver name. 
	@param[out] a_cdrMax     Maximum density found in one of the CDR solvers
	@param[out] a_solverName The solver name
      */
      virtual void
      getCdrMax(Real& a_cdrMax, std::string& a_solverName) const;

      /*!
	@brief Initialize the surface charge. 
      */
      virtual void
      initialSigma();

      /*!
	@brief Project a flux onto the EB normal. This is the AMR version. 
	@param[out] a_projectedFlux Projected flux -- must have one component. 
	@param[in]  a_flux          Vector flux -- must have SpaceDim components.
      */
      virtual void
      projectFlux(EBAMRIVData& a_projectedFlux, const EBAMRIVData& a_flux);

      /*!
	@brief Project a flux onto the EB normal. This is the level version. 
	@param[out] a_projectedFlux Projected flux -- must have one component. 
	@param[in]  a_flux          Vector flux -- must have SpaceDim components.
      */
      virtual void
      projectFlux(LevelData<BaseIVFAB<Real>>&       a_projectedFlux,
                  const LevelData<BaseIVFAB<Real>>& a_flux,
                  const int                         a_lvl);

      /*!
	@brief Project flux onto domain edges/faces.
	@param[out] a_projectedFlux Projected flux -- must have one component. 
	@param[in]  a_flux          Vector flux -- must have SpaceDim components.
      */
      virtual void
      projectDomain(EBAMRIFData& a_projectedFlux, const EBAMRIFData& a_flux);

      /*!
	@brief Regrid internal storage. This MUST be overwritten by implementations. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      virtual void
      regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) = 0;

      /*!
	@brief Regrid all of our solvers. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      virtual void
      regridSolvers(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

      /*!
	@brief Set data defined over dielectric cells to zero.
	@param[in] a_data EB-centered data. 
      */
      virtual void
      resetDielectricCells(EBAMRIVData& a_data) const;

      /*!
	@brief Do a sanity check. Only works when DEBUG=TRUE.
      */
      virtual void
      sanityCheck() const;

      /*!
	@brief Setup the CDR solvers. 
      */
      virtual void
      setupCdr();

      /*!
	@brief Setup the field solver
      */
      virtual void
      setupPoisson();

      /*!
	@brief Setup the RTE solvers
      */
      virtual void
      setupRadiativeTransfer();

      /*!
	@brief Set up the surface charge solver. 
      */
      virtual void
      setupSigma();

      /*!
	@brief Set solver verbosity.
	@details This will set the solver verbosity to m_solverVerbosity. 
      */
      virtual void
      setSolverVerbosity();

      /*!
	@brief All solvers writes plot files. These are *not* the same as the plot files written by Driver. 
      */
      virtual void
      solverDump();

      /*!
	@brief Set voltage curve.
	@param[in] a_voltage Voltage curve
      */
      virtual void
      setVoltage(std::function<Real(const Real a_time)> a_voltage);

      /*!
	@brief Solve the radiative transfer problem. 
	@details This will call the other version with states taken from the solvers. 
	@param[in] a_dt Time step.
      */
      virtual void
      solveRadiativeTransfer(const Real a_dt);

      /*!
	@brief Set the CDR layout. 
	@param[in] a_cdr CDR solvers.
      */
      virtual void
      setCdrSolvers(RefCountedPtr<CdrLayout<CdrSolver>>& a_cdr);

      /*!
	@brief Set the field solver
	@param[in] a_fieldSolver Field solver
      */
      virtual void
      setFieldSolver(RefCountedPtr<FieldSolver>& a_fieldSolver);

      /*!
	@brief Set the RTE solvers
	@param[in] a_rte RTE solvers
      */
      virtual void
      setRadiativeTransferSolvers(RefCountedPtr<RtLayout<RtSolver>>& a_rte);

      /*!
	@brief Compute the ohmic current on electrodes.
	@details This will compute the total charge flux through electrode surfaces and return it. 
	@note In 2D this will yield the current per unit length into the domain. 
      */
      virtual Real
      computeElectrodeCurrent();

      /*!
	@brief Compute the ohmic current on dielectrics
	@details This will compute the total charge flux through dielectric surfaces and return it. 
      */
      virtual Real
      computeDielectricCurrent();

      /*!
	@brief Compute the ohmic current through the domain faces
	@details This will compute the total current through domain faces and return it. 
      */
      virtual Real
      computeDomainCurrent();

      /*!
	@brief Compute induced current in external circuit due to ohmic conduction.
	@details This will return int(E*J dV) over the entire domain. 
      */
      virtual Real
      computeOhmicInductionCurrent();

      /*!
	@brief Compute the relaxation time as dt = eps0/conductivity. 
      */
      virtual Real
      computeRelaxationTime();

      /*!
	@brief Get dt
      */
      virtual Real
      getDt();

    protected:
      /*!
	@brief Enum for switching between source term computations. 'Interpolated' interpolates the cell-centered data
	to cell centroids. InterpolatedStable is the same as Interpolated but with a stability fix near inflow EBs.
	CellAverage assumes centroid-centered data, and Upwind uses the upwind formulation in Villa et. al. 
      */
      enum class SourceTermComputation
      {
        Interpolated,
        InterpolatedStable,
        CellAverage,
        Upwind
      };

      /*!
	@brief An enum for encapsulating how time steps were restricted. 
      */
      enum class TimeCode
      {
        Advection,
        Diffusion,
        AdvectionDiffusion,
        Source,
        RelaxationTime,
        Restricted,
        Hardcap,
        Error,
        Physics
      };

      /*!
	@brief Realm where the fluid and field solvers live
      */
      std::string m_fluidRealm;

      /*!
	@brief If using dual-grid, RTE solvers live here.
      */
      std::string m_particleRealm;

      /*!
	@brief Time stepper class name
      */
      std::string m_className;

      /*!
	@brief Plasma phase
      */
      phase::which_phase m_phase;

      /*!
	@brief Which source term computation
      */
      SourceTermComputation m_whichSourceTermComputation;

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
      RefCountedPtr<CdrLayout<CdrSolver>> m_cdr;

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
      RefCountedPtr<SurfaceODESolver<1>> m_sigma;

      /*!
	@brief Storage for the current
      */
      EBAMRCellData m_currentDensity;

      /*!
	@brief Storage for physics-observables
      */
      EBAMRCellData m_physicsPlotVars;

      /*!
	@brief Set the potential
      */
      std::function<Real(const Real a_time)> m_voltage;

      /*!
	@brief Minimum allowed time step
      */
      Real m_minDt;

      /*!
	@brief Maximum allowed time step
      */
      Real m_maxDt;

      /*!
	@brief CFL number.
      */
      Real m_cfl;

      /*!
	@brief Fudge factor for relaxation time
      */
      Real m_relaxTime;

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
      Real m_dtCFL;

      /*!
	@brief Time code for step restriction
      */
      TimeCode m_timeCode;

      /*!
	@brief Verbosity for solvers
      */
      int m_solverVerbosity;

      /*!
	@brief Fast RTE solver
      */
      int m_fastRTE;

      /*!
	@brief Fast RTE solver
      */
      int m_fastPoisson;

      /*!
	@brief Upwind factor
      */
      int m_upwindFactor;

      /*!
	@brief Parse class verbosity
      */
      virtual void
      parseVerbosity();

      /*!
	@brief Parse solver verbosities. 
      */
      virtual void
      parseSolverVerbosity();

      /*!
	@brief Parse the CFL number. 
      */
      virtual void
      parseCFL();

      /*!
	@brief Parse the relaxation time restriction. 
      */
      virtual void
      parseRelaxationTime();

      /*!
	@brief Parse the minimum allowed time step
      */
      virtual void
      parseMinDt();

      /*!
	@brief Parse the maximum allowed time step
      */
      virtual void
      parseMaxDt();

      /*!
	@brief Parse the "fast poisson" method, i.e. how often we solve the Poisson equation. 
      */
      virtual void
      parseFastPoisson();

      /*!
	@brief Parse the "fast rte" method, i.e. how often we solve the Poisson equation. 
      */
      virtual void
      parseFastRadiativeTransfer();

      /*!
	@brief Parse how we compute the source terms. 
	@details This is for modifications near the EB.
      */
      virtual void
      parseSourceComputation();

      /*!
	@brief Compute and put the current density in output data holder
	@param[inout] a_output Output data holder.
	@param[inout] a_icomp  Starting component in a_output to begin at. On output this is incremented by SpaceDim.
	@param[in]    a_outputRealm Realm where a_output lives    
	@param[in]    a_level  Grid level
      */
      virtual void
      writeJ(LevelData<EBCellFAB>& a_output, int& a_comp, const std::string a_outputRealm, const int a_level) const;

      /*!
	@brief Calculate the user-supplied plot variables
      */
      virtual void
      computePhysicsPlotVars(EBAMRCellData& a_physicsPlotVars) const noexcept;

      /*!
	@brief Write data to output. Convenience function. 
	@param[inout] a_output Output data holder.
	@param[inout] a_icomp  Starting component where this solver begins writing the output. 
	@param[in]    a_data   Data to write. 
	@param[in]    a_outputRealm Realm where a_output lives
	@param[in]    a_level  Grid level
	@param[in]    a_interpToCentroids If true, a_data will be interpolated to cell centroids before writing to a_output. 
	@param[in]    a_interpGhost If true, interpolate ghost cells
      */
      virtual void
      writeData(LevelData<EBCellFAB>& a_output,
                int&                  a_comp,
                const EBAMRCellData&  a_data,
                const std::string     a_outputRealm,
                const int             a_level,
                const bool            a_interpToCentroids,
                const bool            a_interpGhost) const noexcept;
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
