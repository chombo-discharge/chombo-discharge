/*!
  @file   sisdc_tagger.H
  @brief  Declaration of a  class for tagging streamer heads
  @author Robert Marskar
  @date   Oct. 2018
*/

#ifndef _SISDC_TAGGER_H
#define _SISDC_TAGGER_H

#include "cell_tagger.H"

/*!
  @brief Class that tags streamer heads based on the electron source term. 
  @note  The electrons must be the first class that is defined in the cdr solver layout
*/
class sisdc_tagger : public cell_tagger {
public:
  
  sisdc_tagger();
  ~sisdc_tagger();

  void compute_tracers();
  
protected:

  int m_cdr_idx;
  Real m_err_thresh;
  Real m_curv_thresh;
  Real m_mag_thresh;

  bool coarsen_cell(const RealVect&         a_pos,
		    const Real&             a_time,
		    const Real&             a_dx,
		    const int&              a_lvl,
		    const Vector<Real>&     a_tracer,
		    const Vector<RealVect>& a_grad_tracer);

  bool refine_cell(const RealVect&         a_pos,
		   const Real&             a_time,
		   const Real&             a_dx,
		   const int&              a_lvl,
		   const Vector<Real>&     a_tracer,
		   const Vector<RealVect>& a_grad_tracer);
#include "CD_NamespaceFooter.H"

#endif
