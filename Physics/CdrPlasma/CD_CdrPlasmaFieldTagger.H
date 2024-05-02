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
      @brief Lightweight tagging class which flags cells for refinement/coarsening using the electric field only. 
      @details This class fetches the point-wise electric field. The class is pure, so users must implement 
      the actual tagging criterion themselves. 
    */
    class CdrPlasmaFieldTagger : public CdrPlasmaTagger
    {
    public:
      /*!
	@brief Weak constructor, must subsequently call define
      */
      CdrPlasmaFieldTagger();

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaFieldTagger();

      /*!
	@brief Compute tracer fields. 
	@details Implements parent class method but introduces a cell-wise method CdrPlasmaFieldTagger::tracer which is used on per-cell basis. 
      */
      void
      computeTracers() const override final;

      /*!
	@brief Parse options. Must be implemented by users. 
      */
      virtual void
      parseOptions() override = 0;

    protected:
      /*!
	@brief Scratch storage
      */
      mutable EBAMRCellData m_scratch;

      /*!
	@brief Storage for electric field
      */
      mutable EBAMRCellData m_electricField;

      /*!
	@brief Storage for electric field
      */
      mutable EBAMRCellData m_gradElectricField;

      /*!
	@brief Allocate memory for electric field and scratch. 
      */
      virtual void
      allocateStorage() const;

      /*!
	@brief Deallocate memory for electric field and scratch. 
      */
      virtual void
      deallocateStorage() const;

      /*!
	@brief Compute electric field onto the input data holders. 
	@param[out] a_electricField         Cell-centered electric field
	@param[out] a_gradientElectricField Gradient of |E| where E is the electric field. 
      */
      virtual void
      computeElectricField(EBAMRCellData& a_electricField, EBAMRCellData& a_gradientElectricField) const;

      /*!
	@brief Compute tracer field. 
	@param[in] a_pos                  Physical coordinates
	@param[in] a_time                 Time
	@param[in] a_electricField        Electric field (in SI units)
	@param[in] a_minElectricField     Minimum value of electric field
	@param[in] a_maxElectricField     Maximum value of electric field
	@param[in] a_gradElectricField    Gradient of |E|, i.e. grad(|E|)
	@param[in] a_minGradElectricField Maximum value of |grad(|E|)|
	@param[in] a_maxGradElectricField Maximum value of |grad(|E|)|
	@return This should return a Vector of values whose length is equal to the number of tracer fields the subclasses uses. 
      */
      virtual Vector<Real>
      tracer(const RealVect a_pos,
             const Real     a_time,
             const Real     a_dx,
             const RealVect a_electricField,
             const Real     a_minElectricField,
             const Real     a_maxElectricField,
             const RealVect a_gradElectricField,
             const Real     a_minGradElectricField,
             const Real     a_maxGradElectricField) const = 0;

      /*!
	@brief Coarsen a cell based on a tracer field
	@param[in] a_pos         Physical position
	@param[in] a_time        Current time
	@param[in] a_dx          Grid resolution
	@param[in] a_lvl         Grid level
	@param[in] a_tracers     Tracer fields
	@param[in] a_gradTracers Gradient of tracer fields
	@return Return true if the cell should be coarsened and false otherwise. 
      */
      virtual bool
      coarsenCell(const RealVect         a_pos,
                  const Real             a_time,
                  const Real             a_dx,
                  const int              a_lvl,
                  const Vector<Real>     a_tracers,
                  const Vector<RealVect> a_gradTracers) const override = 0;

      /*!
	@brief Refine a cell based on a tracer field
	@param[in] a_pos         Physical position
	@param[in] a_time        Current time
	@param[in] a_dx          Grid resolution
	@param[in] a_lvl         Grid level
	@param[in] a_tracers     Tracer fields
	@param[in] a_gradTracers Gradient of tracer fields
	@return True if the cell should be refined and false otherwise. 
      */
      virtual bool
      refineCell(const RealVect         a_pos,
                 const Real             a_time,
                 const Real             a_dx,
                 const int              a_lvl,
                 const Vector<Real>     a_tracers,
                 const Vector<RealVect> a_gradTracers) const override = 0;
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
