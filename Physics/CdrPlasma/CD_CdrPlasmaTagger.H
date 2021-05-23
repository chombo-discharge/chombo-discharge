/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaTagger.H
  @brief  Class for tagging cells in a Physics/CdrPlasma context
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
      @brief Abstract class for tagging cells across an AMR hierarchy
    */
    class CdrPlasmaTagger : public CellTagger{
    public:

      /*!
	@brief Weak constructor, must subsequently call define
      */
      CdrPlasmaTagger();

      /*!
	@brief Full constructor
      */
      CdrPlasmaTagger(const RefCountedPtr<CdrPlasmaPhysics>&     a_physics,
		      const RefCountedPtr<CdrPlasmaStepper>&     a_timeStepper,
		      const RefCountedPtr<AmrMesh>&               a_amr,
		      const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaTagger();

      /*!
	@brief Define function
	@param[in] a_physics     Plasma kinetics
	@param[in] a_timeStepper Time stepper
	@param[in] a_amr         AmrMesh instance
	@param[in] a_computationalGeometry    Computational geometry
      */
      virtual void define(const RefCountedPtr<CdrPlasmaPhysics>&     a_physics,
			  const RefCountedPtr<CdrPlasmaStepper>&     a_timeStepper,
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
	@brief Get number of plot variables. 
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

#include <CD_NamespaceFooter.H>

#endif
