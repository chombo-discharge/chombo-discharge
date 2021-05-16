/*!
  @file   ito_plasma_tagger.H
  @brief  Shameless ripoff of cdr_plasma_tagger.H
  @author Robert Marskar
  @date   June 2020. 
  @todo   Rename to ito_plasma_tagger. 
*/

#ifndef _ITO_PLASMA_TAGGER_
#define _ITO_PLASMA_TAGGER_

#include "CD_CellTagger.H"
#include "ito_plasma_stepper.H"
#include "ito_plasma_physics.H"

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace ito_plasma {
  
    /*!
      @brief Abstract class for tagging cells across an AMR hierarchy
    */
    class ito_plasma_tagger : public CellTagger{
    public:

      /*!
	@brief Weak constructor, must subsequently call define
      */
      ito_plasma_tagger();

      /*!
	@brief Full constructor
      */
      ito_plasma_tagger(const RefCountedPtr<ito_plasma_physics>&     a_physics,
			const RefCountedPtr<ito_plasma_stepper>&     a_timeStepper,
			const RefCountedPtr<AmrMesh>&               a_amr,
			const RefCountedPtr<computational_geometry>& a_computationalGeometry);

      /*!
	@brief Destructor
      */
      virtual ~ito_plasma_tagger();

      /*!
	@brief Define function
	@param[in] a_physics     Plasma kinetics
	@param[in] a_timeStepper Time stepper
	@param[in] a_amr         AmrMesh instance
	@param[in] a_computationalGeometry    Computational geometry
      */
      virtual void define(const RefCountedPtr<ito_plasma_physics>&     a_physics,
			  const RefCountedPtr<ito_plasma_stepper>&     a_timeStepper,
			  const RefCountedPtr<AmrMesh>&               a_amr,
			  const RefCountedPtr<computational_geometry>& a_computationalGeometry);


      /*!
	@brief Compute tracer fields
      */
      virtual void compute_tracers() = 0;

      /*!
	@brief Parse class options
      */
      virtual void parseOptions() = 0;

      /*!
	@brief Parse class options
      */
      virtual void parseRuntimeOptions() = 0;

      /*!
	@brief
      */
      virtual int getNumberOfPlotVariables();

      /*!
	@brief Regrid storage for tracer field
      */
      virtual void regrid();

      /*!
	@brief Tag cells
      */
      virtual bool tagCells(EBAMRTags& a_tags);

      /*!
	@brief Set phase
      */
      virtual void setPhase(const phase::which_phase a_phase);

      /*!
	@brief Tracer fields
      */
      virtual Vector<EBAMRCellData>& get_tracer_fields();

      /*!
	@brief Write plot data
      */
      virtual void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp);

    protected:

      /*!
	@brief Realm
      */
      std::string m_Realm;

      /*!
	@brief Phase
      */
      phase::which_phase m_phase;
  
      /*!
	@brief Plasma kinetics. 
      */
      RefCountedPtr<ito_plasma_physics> m_physics;

      /*!
	@brief Time stepper
      */
      RefCountedPtr<ito_plasma_stepper> m_timeStepper;

      /*!
	@brief AmrMesh. Needed for computing tracer stuff
      */
      RefCountedPtr<AmrMesh> m_amr;

      /*!
	@brief Computational geometry. Needed for what exactly?
      */
      RefCountedPtr<computational_geometry> m_computationalGeometry;

      /*!
	@brief Number of tracer fields
      */
      int m_num_tracers;

      /*!
	@brief Tracer fields
      */
      Vector<EBAMRCellData> m_tracer;

      /*!
	@brief Gradient of tracer fields
      */
      Vector<EBAMRCellData> m_grad_tracer;

      /*!
	@brief Per-box refinement tags
      */
      virtual void refine_cells_box(DenseIntVectSet&          a_refined_tags,
				    const Vector<EBCellFAB*>& a_tracers,
				    const Vector<EBCellFAB*>& a_grad_tracers,
				    const int                 a_lvl,
				    const Box                 a_box,
				    const EBISBox&            a_ebisbox,
				    const Real                a_time,
				    const Real                a_dx,
				    const RealVect            a_origin);

      /*!
	@brief Per-box coarsen
      */
      virtual void coarsen_cells_box(DenseIntVectSet&          a_refined_tags,
				     const Vector<EBCellFAB*>& a_tracers,
				     const Vector<EBCellFAB*>& a_grad_tracers,
				     const int                 a_lvl,
				     const Box                 a_box,
				     const EBISBox&            a_ebisbox,
				     const Real                a_time,
				     const Real                a_dx,
				     const RealVect            a_origin);



      /*!
	@brief Coarsen a cell based on a tracer field
      */
      virtual bool coarsen_cell(const RealVect         a_pos,
				const Real             a_time,
				const Real             a_dx,
				const int              a_lvl,
				const Vector<Real>     a_tracer,
				const Vector<RealVect> a_grad_tracer) = 0;

      /*!
	@brief Refine a cell based on a tracer field
      */
      virtual bool refine_cell(const RealVect         a_pos,
			       const Real             a_time,
			       const Real             a_dx,
			       const int              a_lvl,
			       const Vector<Real>     a_tracer,
			       const Vector<RealVect> a_grad_tracer) = 0;
    };
  }
}
#include "CD_NamespaceFooter.H"
#endif
