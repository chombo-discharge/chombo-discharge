/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaFieldTagger.H
  @brief  Declaration of an abstract field-only tagging class for ito plasmas
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaFieldTagger_H
#define CD_ItoPlasmaFieldTagger_H

// Our includes
#include <CD_ItoPlasmaTagger.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    /*!
      @brief Lightweight tagging class. Only tags cells based on electric field. 
      @details This class fetches the point-wise electric field. The class is pure, so users must implement 
      the actual tagging criterion themselves. 
    */
    template <typename S>
    class ItoPlasmaFieldTagger : public ItoPlasmaTagger<S>
    {
    public:
      /*!
	@brief Weak constructor, must subsequently call define
      */
      ItoPlasmaFieldTagger();

      /*!
	@brief Destructor
      */
      virtual ~ItoPlasmaFieldTagger();

      /*!
	@brief Compute tracer fields
      */
      virtual void
      computeTracers();

      /*!
	@brief Parse options
      */
      virtual void
      parseOptions() = 0;

      /*!
	@brief Parse options
      */
      virtual void
      parseRuntimeOptions() = 0;

    protected:
      /*!
	@brief Scratch storage
      */
      EBAMRCellData m_scratch;

      /*!
	@brief Storage for electric field
      */
      EBAMRCellData m_E;

      /*!
	@brief Storage for electric field
      */
      EBAMRCellData m_grad_E;

      /*!
	@brief Allocate memory
      */
      virtual void
      allocateStorage();

      /*!
	@brief Deallocate memory
      */
      virtual void
      deallocateStorage();

      /*!
	@brief Compute electric field
      */
      virtual void
      computeElectricField(EBAMRCellData& a_E, EBAMRCellData& a_grad_E);

      /*!
	@brief Compute tracer field
      */
      virtual Vector<Real>
      tracer(const RealVect a_pos,
             const Real     a_time,
             const Real     a_dx,
             const RealVect a_E,
             const Real     a_min_E,
             const Real     a_max_E,
             const RealVect a_grad_E,
             const Real     a_min_grad_E,
             const Real     a_max_grad_E) = 0;

      /*!
	@brief Coarsen a cell based on a tracer field
      */
      virtual bool
      coarsenCell(const RealVect         a_pos,
                  const Real             a_time,
                  const Real             a_dx,
                  const int              a_lvl,
                  const Vector<Real>     a_tracer,
                  const Vector<RealVect> a_grad_tracer) = 0;

      /*!
	@brief Refine a cell based on a tracer field
      */
      virtual bool
      refineCell(const RealVect         a_pos,
                 const Real             a_time,
                 const Real             a_dx,
                 const int              a_lvl,
                 const Vector<Real>     a_tracer,
                 const Vector<RealVect> a_grad_tracer) = 0;
    };
  } // namespace ItoPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoPlasmaFieldTaggerImplem.H>

#endif
