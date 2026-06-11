/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ItoKMCFieldTagger.H
  @brief  Declaration of the Physics::ItoKMC::ItoKMCFieldTagger abstract CellTagger
  @author Robert Marskar
*/

#ifndef CD_ITOKMCFIELDTAGGER_H
#define CD_ITOKMCFIELDTAGGER_H

// Our includes
#include <CD_ItoKMCTagger.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /**
      @brief Lightweight tagging class. Only tags cells based on the electric field. 
      @details This class fetches the point-wise electric field. The class is pure, so users must still implement 
      the actual tagging criterion themselves. However, this class at the very least exposes the per-cell refinement/coarsening
      in terms of the electric field only.
    */
    template <typename S>
    class ItoKMCFieldTagger : public ItoKMCTagger<S>
    {
    public:
      /**
	@brief Weak constructor. User MUST subsequently call define.
      */
      ItoKMCFieldTagger() noexcept;

      /**
	@brief Destructor
      */
      virtual ~ItoKMCFieldTagger() noexcept;

      /**
	@brief Compute tagging fields. 
      */
      virtual void
      computeTagFields() const noexcept override;

      /**
	@brief Parse options
      */
      virtual void
      parseOptions() override = 0;

      /**
	@brief Parse run-time configurable options
      */
      virtual void
      parseRuntimeOptions() override = 0;

    protected:
      /**
	@brief Scratch storage
      */
      mutable EBAMRCellData m_scratch;

      /**
	@brief Storage for the electric field E
      */
      mutable EBAMRCellData m_E;

      /**
	@brief Storage for grad(|E|)
      */
      mutable EBAMRCellData m_gradE;

      /**
	@brief Allocate memory for scratch, electric field, and gradient of electric field. 
      */
      virtual void
      allocateStorage() const noexcept;

      /**
	@brief Deallocate memory
      */
      virtual void
      deallocateStorage() const noexcept;

      /**
	@brief Compute the electric field and its gradient magnitude.
	@param[out] a_E     Cell-centred electric field vector.
	@param[out] a_gradE Cell-centred gradient of |E|.
      */
      virtual void
      computeElectricField(EBAMRCellData& a_E, EBAMRCellData& a_gradE) const noexcept;

      /**
	@brief Compute per-cell tagging field values from the electric field and its gradient.
	@param[in] a_pos      Physical coordinates.
	@param[in] a_time     Current simulation time.
	@param[in] a_dx       Grid resolution.
	@param[in] a_E        Electric field vector at this cell.
	@param[in] a_minE     Minimum |E| over the domain.
	@param[in] a_maxE     Maximum |E| over the domain.
	@param[in] a_gradE    grad(|E|) at this cell.
	@param[in] a_minGradE Minimum |grad(|E|)| over the domain.
	@param[in] a_maxGradE Maximum |grad(|E|)| over the domain.
	@return Per-cell tagging field values.
      */
      virtual Vector<Real>
      computeTagFields(const RealVect& a_pos,
                       const Real      a_time,
                       const Real      a_dx,
                       const RealVect& a_E,
                       const Real      a_minE,
                       const Real      a_maxE,
                       const RealVect& a_gradE,
                       const Real      a_minGradE,
                       const Real      a_maxGradE) const noexcept = 0;

      /**
	@brief Determine if a particular cell should be coarsened or not. 
	@param[in] a_pos Physical position. 
	@param[in] a_time Current time
	@param[in] a_dx Grid resolution
	@param[in] a_lvl Grid level 
	@param[in] a_tagFields Cell tagging fields
	@param[in] a_gradTagFields Gradient of cell tagging fields
        @return Return value
      */
      virtual bool
      coarsenCell(const RealVect&         a_pos,
                  const Real              a_time,
                  const Real              a_dx,
                  const int               a_lvl,
                  const Vector<Real>&     a_tagFields,
                  const Vector<RealVect>& a_gradTagFields) const noexcept override = 0;

      /**
	@brief Determine if a particular cell should be refined or not. 
	@param[in] a_pos Physical position. 
	@param[in] a_time Current time
	@param[in] a_dx Grid resolution
	@param[in] a_lvl Grid level 
	@param[in] a_tagFields Cell tagging fields
	@param[in] a_gradTagFields Gradient of cell tagging fields
        @return Return value
      */
      virtual bool
      refineCell(const RealVect&         a_pos,
                 const Real              a_time,
                 const Real              a_dx,
                 const int               a_lvl,
                 const Vector<Real>&     a_tagFields,
                 const Vector<RealVect>& a_gradTagFields) const noexcept override = 0;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCFieldTaggerImplem.H>

#endif
