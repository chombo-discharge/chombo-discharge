/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_AdvectionDiffusionTagger.H
  @brief  Declaration of a cell tagger class for AdvectionDiffusion physics.
  @author Robert Marskar
*/

#ifndef CD_AdvectionDiffusionTagger_H
#define CD_AdvectionDiffusionTagger_H

#include <CD_CellTagger.H>
#include <CD_CdrSolver.H>
#include <CD_AmrMesh.H>

#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace AdvectionDiffusion {

    /*!
      @brief Class for tagging cells in AdvectionDiffusion physics
      @details This class refines cells if grad(phi)*dx/phi > r1 or if phi > r2.
    */
    class AdvectionDiffusionTagger : public CellTagger {
    public:

      /*!
	@brief Disallowed weak construction.
      */
      AdvectionDiffusionTagger() = delete;

      /*!
	@brief Constructor 
	@param[in] a_solver Cdr solver reference -- needed in order to fetch data used for tagging
	@param[in] a_amr    AmrMesh solver reference
      */
      AdvectionDiffusionTagger(RefCountedPtr<CdrSolver>& a_solver,
			       RefCountedPtr<AmrMesh>&   a_amr);

      /*!
	@brief Destructor (does nothing)
      */
      virtual ~AdvectionDiffusionTagger();

      /*!
	@brief Regrid method.
	@details Does nothing, but is overriden because parent method is pure.
      */
      void regrid() override;

      /*!
	@brief Parse options
	@details Parses m_refCurv and m_refMagn.
      */
      void parseOptions() override;

      /*!
	@brief Flag cells for refinement.
	@param[inout] a_tags Cell tags across AMR hierarchy. See CellTagger for details.
      */
      bool tagCells(EBAMRTags& a_tags) override;

    protected:

      /*!
	@brief Realm where solver is registered.
      */
      std::string m_realm;

      /*!
	@brief Reference to CdrSolver
      */
      RefCountedPtr<CdrSolver> m_solver;

      /*!
	@brief Reference to AmrMesh
      */      
      RefCountedPtr<AmrMesh> m_amr;

      /*!
	@brief Refinement threshold for curvature-based refinement.
      */
      Real m_refCurv;
      
      /*!
	@brief Refinement threshold for magnitude-based refinement.
      */      
      Real m_refMagn;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
