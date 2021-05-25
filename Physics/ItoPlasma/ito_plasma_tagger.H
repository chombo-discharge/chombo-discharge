/*!
  @file   ito_plasma_tagger.H
  @brief  Shameless ripoff of CdrPlasmaTagger.H
  @author Robert Marskar
  @date   June 2020. 
  @todo   Rename to ito_plasma_tagger. 
*/

#ifndef _ITO_PLASMA_TAGGER_
#define _ITO_PLASMA_TAGGER_

#include "CD_CellTagger.H"
#include <CD_ItoPlasmaStepper.H>
#include <CD_ItoPlasmaPhysics.H>

#include "CD_NamespaceHeader.H"
namespace Physics {
  namespace ItoPlasma {
  
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
      ito_plasma_tagger(const RefCountedPtr<ItoPlasmaPhysics>&     a_physics,
			const RefCountedPtr<ItoPlasmaStepper>&     a_timeStepper,
			const RefCountedPtr<AmrMesh>&               a_amr,
			const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

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
      virtual void define(const RefCountedPtr<ItoPlasmaPhysics>&     a_physics,
			  const RefCountedPtr<ItoPlasmaStepper>&     a_timeStepper,
			  const RefCountedPtr<AmrMesh>&               a_amr,
			  const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);


      /*!
	@brief Compute tracer fields
      */
      virtual void computeTracers() = 0;

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
      virtual Vector<EBAMRCellData>& getTracerFields();

      /*!
	@brief Write plot data
      */
      virtual void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp);

    protected:

      /*!
	@brief Realm
      */
      std::string m_realm;

      /*!
	@brief Phase
      */
      phase::which_phase m_phase;
  
      /*!
	@brief Plasma kinetics. 
      */
      RefCountedPtr<ItoPlasmaPhysics> m_physics;

      /*!
	@brief Time stepper
      */
      RefCountedPtr<ItoPlasmaStepper> m_timeStepper;

      /*!
	@brief AmrMesh. Needed for computing tracer stuff
      */
      RefCountedPtr<AmrMesh> m_amr;

      /*!
	@brief Computational geometry. Needed for what exactly?
      */
      RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

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
      virtual void refineCellsBox(DenseIntVectSet&          a_refined_tags,
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
      virtual void coarsenCellsBox(DenseIntVectSet&          a_refined_tags,
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
      virtual bool coarsenCell(const RealVect         a_pos,
				const Real             a_time,
				const Real             a_dx,
				const int              a_lvl,
				const Vector<Real>     a_tracer,
				const Vector<RealVect> a_grad_tracer) = 0;

      /*!
	@brief Refine a cell based on a tracer field
      */
      virtual bool refineCell(const RealVect         a_pos,
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
