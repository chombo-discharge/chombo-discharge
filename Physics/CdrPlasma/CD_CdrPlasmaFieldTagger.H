/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaFieldTagger.H
  @brief  Cell tagger class for Physics/CdrPlasma. 
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaFieldTagger_H
#define CD_CdrPlasmaFieldTagger_H

// Our includes
#include <CD_CdrPlasmaTagger.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {
  
    /*!
      @brief Lightweight tagging class. Only tags cells based on electric field. 
      @details This class fetches the point-wise electric field. The class is pure, so users must implement 
      the actual tagging criterion themselves. 
    */
    class CdrPlasmaFieldTagger : public CdrPlasmaTagger {
    public:

      /*!
	@brief Weak constructor, must subsequently call define
      */
      CdrPlasmaFieldTagger();

      /*!
	@brief Destructor
      */
      ~CdrPlasmaFieldTagger();

      /*!
	@brief Compute tracer fields
      */
      void computeTracers() const override final;

      /*!
	@brief Parse options
      */
      virtual void parseOptions() = 0;

    protected:

      /*!
	@brief Scratch storage
      */
      mutable EBAMRCellData m_scratch;
  
      /*!
	@brief Storage for electric field
      */
      mutable EBAMRCellData m_E;

      /*!
	@brief Storage for electric field
      */
      mutable EBAMRCellData m_grad_E;

      /*!
	@brief Allocate memory
      */
      void allocateStorage() const;

      /*!
	@brief Deallocate memory
      */
      void deallocateStorage() const;

      /*!
	@brief Compute electric field
      */
      void computeElectricField(EBAMRCellData& a_E, EBAMRCellData& a_grad_E) const;
    
      /*!
	@brief Compute tracer field
      */
      virtual Vector<Real> tracer(const RealVect          a_pos,
				  const Real              a_time,
				  const Real              a_dx,
				  const RealVect          a_E,
				  const Real              a_min_E,
				  const Real              a_max_E,
				  const RealVect          a_grad_E,
				  const Real              a_min_grad_E,
				  const Real              a_max_grad_E) const = 0;


      /*!
	@brief Coarsen a cell based on a tracer field
      */
      virtual bool coarsenCell(const RealVect         a_pos,
			       const Real             a_time,
			       const Real             a_dx,
			       const int              a_lvl,
			       const Vector<Real>     a_tracer,
			       const Vector<RealVect> a_grad_tracer) const = 0;

      /*!
	@brief Refine a cell based on a tracer field
      */
      virtual bool refineCell(const RealVect         a_pos,
			      const Real             a_time,
			      const Real             a_dx,
			      const int              a_lvl,
			      const Vector<Real>     a_tracer,
			      const Vector<RealVect> a_grad_tracer) const = 0;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
