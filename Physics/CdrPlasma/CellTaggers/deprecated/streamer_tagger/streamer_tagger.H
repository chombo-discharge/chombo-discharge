/*!
  @file   CdrPlasmaStreamerTagger.H
  @brief  Declaration of a  class for tagging streamer heads
  @author Robert Marskar
  @date   Oct. 2018
*/

#ifndef _STREAMER_TAGGER_H
#define _STREAMER_TAGGER_H

#include "CD_CellTagger.H"

/*!
  @brief Class that tags streamer heads based on the electron source term. 
  @note  The electrons must be the first class that is defined in the cdr solver layout
*/
class CdrPlasmaStreamerTagger : public CellTagger {
public:
  
  CdrPlasmaStreamerTagger();
  ~CdrPlasmaStreamerTagger();

  void computeTracers();
  
protected:

  Real m_fudge;
  int m_electron_idx;
  Vector<Real> m_thresh1;
  Vector<Real> m_thresh2;



  bool coarsenCell(const RealVect&         a_pos,
		    const Real&             a_time,
		    const Real&             a_dx,
		    const int&              a_lvl,
		    const Vector<Real>&     a_tracer,
		    const Vector<RealVect>& a_grad_tracer);

  bool refineCell(const RealVect&         a_pos,
		   const Real&             a_time,
		   const Real&             a_dx,
		   const int&              a_lvl,
		   const Vector<Real>&     a_tracer,
		   const Vector<RealVect>& a_grad_tracer);
#include "CD_NamespaceFooter.H"

#endif
