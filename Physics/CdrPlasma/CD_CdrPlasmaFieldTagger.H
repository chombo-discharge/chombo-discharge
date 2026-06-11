/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_CdrPlasmaFieldTagger.H
   @brief  Declaration of the Physics::CdrPlasma::CdrPlasmaFieldTagger class
   @author Robert Marskar
*/

#ifndef CD_CDRPLASMAFIELDTAGGER_H
#define CD_CDRPLASMAFIELDTAGGER_H

// Our includes
#include <CD_CdrPlasmaTagger.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /**
       @brief Intermediate CdrPlasmaTagger subclass that computes tracer fields from the electric field only.
       @details This class implements computeTracers() by fetching the electric field and its gradient
       magnitude from the time stepper and calling the pure virtual per-cell tracer() method on every
       grid cell. Users subclass this and implement tracer(), refineCell(), and coarsenCell() to define
       their specific tagging criteria based on the local electric field.
    */
    class CdrPlasmaFieldTagger : public CdrPlasmaTagger
    {
    public:
      /**
         @brief Default constructor. Sets the tagger name.
      */
      CdrPlasmaFieldTagger();

      /**
         @brief Destructor.
      */
      virtual ~CdrPlasmaFieldTagger();

      /**
         @brief Compute tracer fields from the electric field on the AMR hierarchy.
         @details Fetches the cell-centered electric field and its gradient, then calls
         the per-cell tracer() method on every regular and irregular cell. Tracer fields
         and their gradients are stored in m_tracers and m_gradTracers.
      */
      void
      computeTracers() const override final;

      /**
         @brief Parse user-defined options. Must be implemented by subclasses.
      */
      virtual void
      parseOptions() override = 0;

    protected:
      /**
         @brief Scratch storage for scalar intermediate quantities.
      */
      mutable EBAMRCellData m_scratch;

      /**
         @brief Cell-centered electric field vector (SpaceDim components).
      */
      mutable EBAMRCellData m_electricField;

      /**
         @brief Cell-centered gradient of the electric field magnitude, i.e., grad(|E|).
      */
      mutable EBAMRCellData m_gradElectricField;

      /**
         @brief Allocate transient storage for the electric field, its gradient, and scratch.
      */
      virtual void
      allocateStorage() const;

      /**
         @brief Release transient storage allocated by allocateStorage().
      */
      virtual void
      deallocateStorage() const;

      /**
         @brief Compute the centroid-centered electric field and its gradient magnitude.
         @param[out] a_electricField         Cell-centered electric field vector.
         @param[out] a_gradientElectricField Cell-centered gradient of |E|.
      */
      virtual void
      computeElectricField(EBAMRCellData& a_electricField, EBAMRCellData& a_gradElectricField) const;

      /**
         @brief Per-cell tracer field evaluation. Must be implemented by subclasses.
         @param[in] a_pos                  Cell centroid position.
         @param[in] a_time                 Current simulation time.
         @param[in] a_dx                   Grid cell spacing.
         @param[in] a_electricField        Electric field vector at this cell.
         @param[in] a_minElectricField     Global minimum of |E| over the AMR hierarchy.
         @param[in] a_maxElectricField     Global maximum of |E| over the AMR hierarchy.
         @param[in] a_gradElectricField    Gradient of |E| at this cell.
         @param[in] a_minGradElectricField Global minimum of |grad(|E|)| over the AMR hierarchy.
         @param[in] a_maxGradElectricField Global maximum of |grad(|E|)| over the AMR hierarchy.
         @return Vector of tracer field values at this cell. Length must equal m_numTracers.
      */
      virtual Vector<Real>
      tracer(const RealVect& a_pos,
             const Real      a_time,
             const Real      a_dx,
             const RealVect& a_electricField,
             const Real      a_minElectricField,
             const Real      a_maxElectricField,
             const RealVect& a_gradElectricField,
             const Real      a_minGradElectricField,
             const Real      a_maxGradElectricField) const = 0;

      /**
         @brief Decide whether a cell should be coarsened based on tracer fields.
         @param[in] a_pos         Cell centroid position.
         @param[in] a_time        Current simulation time.
         @param[in] a_dx          Grid cell spacing.
         @param[in] a_lvl         AMR level of this cell.
         @param[in] a_tracers     Tracer field values at this cell.
         @param[in] a_gradTracers Gradient of tracer fields at this cell.
         @return True if the cell should be coarsened.
      */
      virtual bool
      coarsenCell(const RealVect&         a_pos,
                  const Real              a_time,
                  const Real              a_dx,
                  const int               a_lvl,
                  const Vector<Real>&     a_tracers,
                  const Vector<RealVect>& a_gradTracers) const override = 0;

      /**
         @brief Decide whether a cell should be refined based on tracer fields.
         @param[in] a_pos         Cell centroid position.
         @param[in] a_time        Current simulation time.
         @param[in] a_dx          Grid cell spacing.
         @param[in] a_lvl         AMR level of this cell.
         @param[in] a_tracers     Tracer field values at this cell.
         @param[in] a_gradTracers Gradient of tracer fields at this cell.
         @return True if the cell should be refined.
      */
      virtual bool
      refineCell(const RealVect&         a_pos,
                 const Real              a_time,
                 const Real              a_dx,
                 const int               a_lvl,
                 const Vector<Real>&     a_tracers,
                 const Vector<RealVect>& a_gradTracers) const override = 0;
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
