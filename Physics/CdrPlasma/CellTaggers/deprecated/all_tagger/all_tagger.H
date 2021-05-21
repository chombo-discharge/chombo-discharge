/*!
  @file   all_tagger
  @brief  Declaration of a class that tags all cells and coarsens none
  @author Robert Marskar
  @date   Aug. 2018
*/

#ifndef _ALL_TAGGER_
#define _ALL_TAGGER_

#include "CD_CellTagger.H"

/*!
  @brief Class that tags all cells
*/
class all_tagger : public CellTagger {
public:

  /*!
    @brief Constructor
  */
  all_tagger(){
    m_num_tracers = 1;
  }
  
  void computeTracers(){
    return;
  }

  bool coarsenCell(const RealVect&         a_pos,
		    const Real&             a_time,
		    const Real&             a_dx,
		    const int&              a_lvl,
		    const Vector<Real>&     a_tracer,
		    const Vector<RealVect>& a_grad_tracer){
    return false;
  }

  bool refineCell(const RealVect&         a_pos,
		   const Real&             a_time,
		   const Real&             a_dx,
		   const int&              a_lvl,
		   const Vector<Real>&     a_tracer,
		   const Vector<RealVect>& a_grad_tracer){
    return true;
  }
#include "CD_NamespaceFooter.H"

#endif
