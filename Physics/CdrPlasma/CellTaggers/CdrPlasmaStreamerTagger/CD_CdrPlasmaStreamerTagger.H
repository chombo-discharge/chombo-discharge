/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_CdrPlasmaStreamerTagger.H
  @brief  Declaration of the CdrPlasmaStreamerTagger class for AMR cell tagging based on streamer criteria
  @author Robert Marskar
*/

#ifndef CD_CDRPLASMASTREAMERTAGGER_H
#define CD_CDRPLASMASTREAMERTAGGER_H

// Our includes
#include <CD_CdrPlasmaFieldTagger.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /**
      @brief Concrete CdrPlasmaFieldTagger that tags cells based on streamer inception criteria.
      @details Cells are flagged for refinement if either the ionization criterion
      alpha(E)*dx exceeds a user-specified threshold m_refiAlpha, or the electric field
      curvature |grad(|E|)|*dx/|E| exceeds m_refiCurv. Coarsening is the reverse of both
      criteria. Both thresholds are runtime-configurable via ParmParse.
    */
    class CdrPlasmaStreamerTagger : public CdrPlasmaFieldTagger
    {
    public:
      /**
	@brief Disallowed weak constructor. 
      */
      CdrPlasmaStreamerTagger() = delete;

      /**
	@brief Full constructor
	@param[in] a_physics               Reference to CdrPlasmaPhysics. 
	@param[in] a_timeStepper           Reference to the TimeStepper being used. 
	@param[in] a_amr                   Reference to AmrMesh. 
	@param[in] a_computationalGeometry Computational geometry
      */
      CdrPlasmaStreamerTagger(const RefCountedPtr<CdrPlasmaPhysics>&      a_physics,
                              const RefCountedPtr<CdrPlasmaStepper>&      a_timeStepper,
                              const RefCountedPtr<AmrMesh>&               a_amr,
                              const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

      /**
	@brief Destructor
      */
      virtual ~CdrPlasmaStreamerTagger();

      /**
	@brief Parse class options
      */
      virtual void
      parseOptions() override;

      /**
	@brief Parse runtime options. 
      */
      virtual void
      parseRuntimeOptions() override;

      /**
	@brief Compute tracer fields. 
	@param[in] a_pos                  Physical coordinates
	@param[in] a_time                 Time
	@param[in] a_electricField        Electric field (in SI units)
	@param[in] a_minElectricField     Minimum value of electric field
	@param[in] a_maxElectricField     Maximum value of electric field
	@param[in] a_gradElectricField    Gradient of |E|, i.e. grad(|E|)
	@param[in] a_minGradElectricField Minimum value of |grad(|E|)| over the AMR hierarchy
	@param[in] a_maxGradElectricField Maximum value of |grad(|E|)|
	@return This should return a Vector of values whose length is equal to the number of tracer fields the subclasses uses. 
        @param[in] a_dx Description
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
             const Real      a_maxGradElectricField) const override;

      /**
	@brief Cell-coarsening method.
	@param[in] a_pos         Physical position
	@param[in] a_time        Current time
	@param[in] a_dx          Grid resolution
	@param[in] a_lvl         Grid level
	@param[in] a_tracers     Tracer fields
	@param[in] a_gradTracers Gradient of tracer fields
	@return Return true if the cell should be coarsened and false otherwise. 
      */
      virtual bool
      coarsenCell(const RealVect&         a_pos,
                  const Real              a_time,
                  const Real              a_dx,
                  const int               a_lvl,
                  const Vector<Real>&     a_tracers,
                  const Vector<RealVect>& a_gradTracers) const override;

      /**
	@brief Cell-refinement method.
	@param[in] a_pos         Physical position
	@param[in] a_time        Current time
	@param[in] a_dx          Grid resolution
	@param[in] a_lvl         Grid level
	@param[in] a_tracers     Tracer fields
	@param[in] a_gradTracers Gradient of tracer fields
	@return True if the cell should be refined and false otherwise. 
      */
      virtual bool
      refineCell(const RealVect&         a_pos,
                 const Real              a_time,
                 const Real              a_dx,
                 const int               a_lvl,
                 const Vector<Real>&     a_tracers,
                 const Vector<RealVect>& a_gradTracers) const override;

    protected:
      /**
	@brief Threshold for coarsening based on curvature
      */
      Real m_coarCurv;

      /**
	@brief Threshold for refinement based on curvature
      */
      Real m_refiCurv;

      /**
	@brief Threshold for refinement based on Townsend ionization coefficient. 
      */
      Real m_refiAlpha;

      /**
	@brief Threshold for coarsening based on the Townsend ionization coefficient.
      */
      Real m_coarAlpha;

      /**
	@brief Maximum AMR level allowed to be coarsened.
      */
      int m_maxCoarsenLevel;
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
