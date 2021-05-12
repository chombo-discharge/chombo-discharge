/*!
  @file   cdr_plasma_tagger.H
  @brief  Class for tagging cells across an AMR hierarchy
  @author Robert Marskar
  @date   Nov. 2017
  @todo   Rename to cdr_plasma_tagger. 
*/

#ifndef _CDR_PLASMA_TAGGER_
#define _CDR_PLASMA_TAGGER_

#include "cell_tagger.H"
#include "cdr_plasma_stepper.H"

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace cdr_plasma {
  
    /*!
      @brief Abstract class for tagging cells across an AMR hierarchy
    */
    class cdr_plasma_tagger : public cell_tagger{
    public:

      /*!
	@brief Weak constructor, must subsequently call define
      */
      cdr_plasma_tagger();

      /*!
	@brief Full constructor
      */
      cdr_plasma_tagger(const RefCountedPtr<cdr_plasma_physics>&     a_physics,
			const RefCountedPtr<cdr_plasma_stepper>&     a_timestepper,
			const RefCountedPtr<AmrMesh>&               a_amr,
			const RefCountedPtr<computational_geometry>& a_compgeom);

      /*!
	@brief Destructor
      */
      virtual ~cdr_plasma_tagger();

      /*!
	@brief Define function
	@param[in] a_physics     Plasma kinetics
	@param[in] a_timestepper Time stepper
	@param[in] a_amr         AmrMesh instance
	@param[in] a_compgeom    Computational geometry
      */
      virtual void define(const RefCountedPtr<cdr_plasma_physics>&     a_physics,
			  const RefCountedPtr<cdr_plasma_stepper>&     a_timestepper,
			  const RefCountedPtr<AmrMesh>&               a_amr,
			  const RefCountedPtr<computational_geometry>& a_compgeom);


      /*!
	@brief Compute tracer fields
      */
      virtual void compute_tracers() = 0;

      /*!
	@brief Parse class options
      */
      virtual void parseOptions() = 0;

      /*!
	@brief
      */
      virtual int get_num_plot_vars();

      /*!
	@brief Regrid storage for tracer field
      */
      virtual void regrid();

      /*!
	@brief Tag cells
      */
      virtual bool tag_cells(EBAMRTags& a_tags);

      /*!
	@brief Set phase
      */
      virtual void set_phase(const phase::which_phase a_phase);

      /*!
	@brief Tracer fields
      */
      virtual Vector<EBAMRCellData>& get_tracer_fields();

      /*!
	@brief Write plot data
      */
      virtual void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp);

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
      RefCountedPtr<cdr_plasma_physics> m_physics;

      /*!
	@brief Time stepper
      */
      RefCountedPtr<cdr_plasma_stepper> m_timestepper;

      /*!
	@brief AmrMesh. Needed for computing tracer stuff
      */
      RefCountedPtr<AmrMesh> m_amr;

      /*!
	@brief Computational geometry. Needed for what exactly?
      */
      RefCountedPtr<computational_geometry> m_compgeom;

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
