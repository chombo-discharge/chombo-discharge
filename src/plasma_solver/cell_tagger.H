/*!
  @file   cell_tagger.H
  @brief  Class for tagging cells across an AMR hierarchy
  @author Robert Marskar
  @date   Nov. 2017
  @todo   Rename to cell_tagger. 
*/

#ifndef _CELL_TAGGER_
#define _CELL_TAGGER_

#include "type_definitions.H"
#include "computational_geometry.H"
#include "physical_domain.H"
#include "time_stepper.H"
#include "plasma_kinetics.H"
#include "amr_mesh.H"

/*!
  @brief Abstract class for tagging cells across an AMR hierarchy
*/
class cell_tagger {
public:

  /*!
    @brief Weak constructor, must subsequently call define
  */
  cell_tagger(const int a_num_tracers = 0);

  /*!
    @brief Destructor
  */
  virtual ~cell_tagger();

  /*!
    @brief Define function
    @param[in] a_plasKin     Plasma kinetics
    @param[in] a_amrParams   AMR parameters
    @param[in] a_timeStepper Time stepper
  */
  virtual void define(const RefCountedPtr<plasma_kinetics>&        a_plaskin,
		      const RefCountedPtr<time_stepper>&           a_timestepper,
		      const RefCountedPtr<amr_mesh>&               a_amr,
		      const RefCountedPtr<computational_geometry>& a_compgeom,
		      const RefCountedPtr<physical_domain>&        a_physdom);


  /*!
    @brief Compute tracer fields
  */
  virtual void compute_tracers() = 0;

  /*!
    @brief Regrid storage for tracer field
  */
  void regrid();

  /*!
    @brief Tag cells
  */
  void tag_cells(EBAMRTags& a_tags);

  /*!
    @brief Set phase
  */
  void set_phase(const phase::which_phase a_phase);

  /*!
    @brief Set verbosity
  */
  void set_verbosity(const int a_verbosity);

  /*!
    @brief Tracer fields
  */
  Vector<EBAMRCellData>& get_tracer_fields();

protected:

    /*!
    @brief Plasma kinetics. 
  */
  RefCountedPtr<plasma_kinetics> m_plaskin;

  /*!
    @brief Time stepper
  */
  RefCountedPtr<time_stepper> m_timestepper;

  /*!
    @brief amr_mesh. Needed for computing tracer stuff
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Computational geometry. Needed for what exactly?
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief Physical domain. 
  */
  RefCountedPtr<physical_domain> m_physdom;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Set name
  */
  std::string m_name;

  /*!
    @brief Number of tracer fields
  */
  int m_num_tracers;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Tracer fields
  */
  Vector<EBAMRCellData> m_tracer;

  /*!
    @brief Gradient of tracer fields
  */
  Vector<EBAMRCellData> m_grad_tracer;

  /*!
    @brief Coarsen a cell based on a tracer field
  */
  virtual bool coarsen_cell(const RealVect&         a_pos,
			    const Real&             a_time,
			    const Real&             a_dx,
			    const int&              a_lvl,
			    const Vector<Real>&     a_tracer,
			    const Vector<RealVect>& a_grad_tracer) = 0;

  /*!
    @brief Refine a cell based on a tracer field
  */
  virtual bool refine_cell(const RealVect&         a_pos,
			   const Real&             a_time,
			   const Real&             a_dx,
			   const int&              a_lvl,
			   const Vector<Real>&     a_tracer,
			   const Vector<RealVect>& a_grad_tracer) = 0;
};

#endif
