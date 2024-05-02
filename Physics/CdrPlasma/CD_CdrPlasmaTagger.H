/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaTagger.H
  @brief  Declearation of an abstract class for tagging cells in a Physics/CdrPlasma context.
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaTagger_H
#define CD_CdrPlasmaTagger_H

// Our includes
#include <CD_CellTagger.H>
#include <CD_CdrPlasmaStepper.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Abstract Class which implements CellTagger for the CDR plasma physics module. 
      @details This class exposes a number of mesh data holders which can be used for flagging cells for refinement. It is only
      meant to be used as an interface, providing the user access to solver data and physics through the constructor. From that, users
      can fetch the necessary data from the solvers (or the physics class) and more easily use the data to flag cells for refinement
      and coarsening. 

      This class uses 'tracer fields' for flagging cells for refinement and coarsening. It is up to the user to compute these fields by implementing
      the routine CdrPlasmaTagger::computeTracerFields. These fields can then be used in the other routines CdrPlasmaTagger::refineCell
      and CdrPlasmaTagger::coarsenCell
    */
    class CdrPlasmaTagger : public CellTagger
    {
    public:
      /*!
	@brief Weak constructor. Must subsequently call define
      */
      CdrPlasmaTagger();

      /*!
	@brief Full constructor. Leaves object in usable state. 
	@param[in] a_physics               Reference to CdrPlasmaPhysics. 
	@param[in] a_timeStepper           Reference to the TimeStepper being used. 
	@param[in] a_amr                   Reference to AmrMesh. 
	@param[in] a_computationalGeometry Computational geometry
      */
      CdrPlasmaTagger(const RefCountedPtr<CdrPlasmaPhysics>&      a_physics,
                      const RefCountedPtr<CdrPlasmaStepper>&      a_timeStepper,
                      const RefCountedPtr<AmrMesh>&               a_amr,
                      const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaTagger();

      /*!
	@brief Define function constructor. Leaves object in usable state. 
	@param[in] a_physics               Reference to CdrPlasmaPhysics. 
	@param[in] a_timeStepper           Reference to the TimeStepper being used. 
	@param[in] a_amr                   Reference to AmrMesh. 
	@param[in] a_computationalGeometry Computational geometry
      */
      void
      define(const RefCountedPtr<CdrPlasmaPhysics>&      a_physics,
             const RefCountedPtr<CdrPlasmaStepper>&      a_timeStepper,
             const RefCountedPtr<AmrMesh>&               a_amr,
             const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

      /*!
	@brief Compute tracer fields.
	@details This should fill m_tracers -- user must implement it. 
      */
      virtual void
      computeTracers() const = 0;

      /*!
	@brief Parse class options -- must be implemented by user. 
      */
      virtual void
      parseOptions() override = 0;

      /*!
	@brief Pre-plot magic for Driver
      */
      virtual void
      prePlot() const noexcept override;

      /*!
	@brief Perform pre-regrid operations
      */
      virtual void
      preRegrid() noexcept override;

      /*!
	@brief Regrid method -- used for reallocating storage for the tracer fields. 
      */
      virtual void
      regrid() override final;

      /*!
	@brief Tag cells for refinement. 
	@param[inout] a_tags Tags on grid levels
	@details EBAMRTags is a data-type Vector<RefCountedPtr<LayoutData<DenseIntVectSet> > >. The vector indicates the grid level, the LayoutData indicates data ownership (in
	much the same way as LevelData. The DenseIntVectSet is essentially an IntVectSet restricted to the patch (i.e. one cannot add IntVects that are outside the patch). 
	This function cells for refinement or coarsening. The user's responsibility is to add (or remove) tags from a_tags. 
      */
      virtual bool
      tagCells(EBAMRTags& a_tags) override final;

      /*!
	@brief Get number of plot variables. This is fetched by Driver and can be used. 
      */
      virtual int
      getNumberOfPlotVariables() const override final;

      /*!
	@brief Get plot variable names. 
      */
      virtual Vector<std::string>
      getPlotVariableNames() const override final;

      /*!
	@brief Write plot data.
	@param[inout] a_output Output data holder
	@param[inout] a_icomp Starting variable in a_output where we begin appending data.
	@pπaram[in] a_outputRealm Realm where a_output lives
	@param[in] a_level Grid level
      */
      virtual void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string     a_outputRealm,
                    const int             a_level) const override final;

    protected:
      /*!
	@brief Realm where this cell tagger allocates its data. 
      */
      std::string m_realm;

      /*!
	@brief Phase where this cell tagger allocates its data. 
      */
      phase::which_phase m_phase;

      /*!
	@brief Plasma physics. 
      */
      RefCountedPtr<CdrPlasmaPhysics> m_physics;

      /*!
	@brief Time stepper
      */
      RefCountedPtr<CdrPlasmaStepper> m_timeStepper;

      /*!
	@brief AmrMesh. Needed for computing tracer stuff
      */
      RefCountedPtr<AmrMesh> m_amr;

      /*!
	@brief Computational geometry. Not needed but useful to have. 
      */
      RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

      /*!
	@brief Number of tracer fields
      */
      int m_numTracers;

      /*!
	@brief Tracer fields
      */
      mutable Vector<EBAMRCellData> m_tracers;

      /*!
	@brief Gradient of tracer fields
      */
      mutable Vector<EBAMRCellData> m_gradTracers;

      /*!
	@brief Per-box refinement method. 
	@details Iterates through all cells and check if they need refinement. 
	@param[out] a_refinedCells Cells flagged for refinement
	@param[in]  a_tracers      Tracer fields on this grid patch. 
	@param[in]  a_gradTracers  Gradient of tracer fields on this grid patch. 
	@param[in]  a_lvl          Grid level
	@param[in]  a_dit          Grid index
	@param[in]  a_box          Cell centered box
	@param[in]  a_ebisbox      EBIS box
	@param[in]  a_time         Current time
	@param[in]  a_dx           Grid resolution
	@param[in]  a_probLo       Lower-left corner of simulation domain. 
      */
      virtual void
      refineCellsBox(DenseIntVectSet&          a_refinedCells,
                     const Vector<EBCellFAB*>& a_tracers,
                     const Vector<EBCellFAB*>& a_gradTracers,
                     const int                 a_lvl,
                     const DataIndex           a_dit,
                     const Box                 a_box,
                     const EBISBox&            a_ebisbox,
                     const Real                a_time,
                     const Real                a_dx,
                     const RealVect            a_probLo);

      /*!
	@brief Per-box coarsening method. 
	@details Iterates through all cells and check if they need coarsening. 
	@param[out] a_coarsenedCells Cells flagged for coarsening.
	@param[in]  a_tracers        Tracer fields on this grid patch. 
	@param[in]  a_gradTracers    Gradient of tracer fields on this grid patch. 
	@param[in]  a_lvl            Grid level
	@param[in]  a_dit          Grid index
	@param[in]  a_box            Cell centered box
	@param[in]  a_ebisbox        EBIS box
	@param[in]  a_time           Current time
	@param[in]  a_dx             Grid resolution
	@param[in]  a_probLo         Lower-left corner of simulation domain. 
      */
      virtual void
      coarsenCellsBox(DenseIntVectSet&          a_coarsenedCells,
                      const Vector<EBCellFAB*>& a_tracers,
                      const Vector<EBCellFAB*>& a_gradTracers,
                      const int                 a_lvl,
                      const DataIndex           a_dit,
                      const Box                 a_box,
                      const EBISBox&            a_ebisbox,
                      const Real                a_time,
                      const Real                a_dx,
                      const RealVect            a_probLo);

      /*!
	@brief Coarsen a cell based on a tracer field
	@param[in] a_pos         Physical position
	@param[in] a_time        Current time
	@param[in] a_dx          Grid resolution
	@param[in] a_lvl         Grid level
	@param[in] a_tracers     Tracer fields
	@param[in] a_gradTracers Gradient of tracer fields
	@return Return true if the cell should be coarsened and false otherwise. 
      */
      virtual bool
      coarsenCell(const RealVect         a_pos,
                  const Real             a_time,
                  const Real             a_dx,
                  const int              a_lvl,
                  const Vector<Real>     a_tracers,
                  const Vector<RealVect> a_gradTracers) const = 0;

      /*!
	@brief Refine a cell based on a tracer field
	@param[in] a_pos         Physical position
	@param[in] a_time        Current time
	@param[in] a_dx          Grid resolution
	@param[in] a_lvl         Grid level
	@param[in] a_tracers     Tracer fields
	@param[in] a_gradTracers Gradient of tracer fields
	@return True if the cell should be refined and false otherwise. 
      */
      virtual bool
      refineCell(const RealVect         a_pos,
                 const Real             a_time,
                 const Real             a_dx,
                 const int              a_lvl,
                 const Vector<Real>     a_tracers,
                 const Vector<RealVect> a_gradTracers) const = 0;
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
