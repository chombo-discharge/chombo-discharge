/*!
  @file   air7_zheleznyak_tagger.H
  @brief  Implementation of CellTagger that tags cells based on the gradient of the electric field
  @author Robert Marskar
  @date   Feb. 2018
*/

#ifndef _AIR7_ZHELEZNYAK_TAGGER_
#define _AIR7_ZHELEZNYAK_TAGGER_

#include "ultralw_tagger.H"

/*!
  @brief Refinement class for morrow_zheleznyak plasma model only
*/
class air7_zheleznyak_tagger : public ultralw_tagger {
public:
  
  /*!
    @brief Constructor
  */
  air7_zheleznyak_tagger();

  /*!
    @brief Destructor
  */
  virtual ~air7_zheleznyak_tagger();

  /*!
    @brief Parse class options
  */
  virtual void parseOptions();

  /*!
    @brief Return tracer field
  */
  virtual Vector<Real> tracer(const RealVect         a_pos,
			      const Real             a_time,
			      const Real             a_dx,
			      const RealVect         a_E,
			      const Real             a_min_E,
			      const Real             a_max_E,
			      const RealVect         a_grad_E,
			      const Real             a_min_grad_E,
			      const Real             a_max_grad_E);

  /*!
    @brief Return coarsening criterion. 
  */
  virtual bool coarsen_cell(const RealVect         a_pos,
			    const Real             a_time,
			    const Real             a_dx,
			    const int              a_lvl,
			    const Vector<Real>     a_tracer,
			    const Vector<RealVect> a_grad_tracer);

  /*!
    @brief Return refinement criterion. 
  */
  virtual bool refine_cell(const RealVect         a_pos,
			   const Real             a_time,
			   const Real             a_dx,
			   const int              a_lvl,
			   const Vector<Real>     a_tracer,
			   const Vector<RealVect> a_grad_tracer);
protected:

  /*!
    @brief Threshold for coarsening based on curvature
  */
  Real m_coar_curv;

  /*!
    @brief Threshold for refinement based on curvature
  */
  Real m_refi_curv;
  
  /*!
    @brief Threshold for ionization refinement
  */
  Real m_refi_alpha;

  /*!
    @brief Threshold for ionization coarsening
  */
  Real m_coar_alpha;

  /*!
    @brief Set coarsening level
  */
  int m_max_coarsen_level;
#include "CD_NamespaceFooter.H"
#endif
