/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BrownianWalkerTagger.H
  @brief  Class for tagging cells across an AMR hierarchy
  @author Robert Marskar
*/

#ifndef CD_BrownianWalkerTagger_H
#define CD_BrownianWalkerTagger_H

// Our includes
#include <CD_CellTagger.H>
#include <CD_ItoSolver.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace BrownianWalker {

    /*!
      @brief Implementation of CellTagger for tagging cells for Brownian walker physics. 
      @details This class just implements the interface but does not flag any cells for refinement. However, it returns true from the tagging routine
      so that TimeStepper and Driver can load balance the application. 
    */
    class BrownianWalkerTagger : public CellTagger
    {
    public:
      /*!
	@brief Disallowed weak construction -- use the full constructor
      */
      BrownianWalkerTagger() = delete;

      /*!
	@brief Full constructor.
	@param[in] a_solver Solver to be tracked
	@param[in] a_amr    Reference to AmrMesh -- needed because we need access to the grids. 
      */
      BrownianWalkerTagger(RefCountedPtr<ItoSolver>& a_solver, RefCountedPtr<AmrMesh>& a_amr);

      /*!
	@brief Destructor (does nothing)
      */
      ~BrownianWalkerTagger();

      /*!
	@brief Regrid method (does nothing)
      */
      void
      regrid() override;

      /*!
	@brief Parse options
      */
      void
      parseOptions() override;

      /*!
	@brief Tag cells (does nothing, returns false)
	@param[inout] a_tags Cell tags
      */
      bool
      tagCells(EBAMRTags& a_tags) override;

    protected:
      /*!
	@brief Reference to solver
      */
      RefCountedPtr<ItoSolver> m_solver;

      /*!
	@brief Reference to AmrMesh
      */
      RefCountedPtr<AmrMesh> m_amr;

      /*!
	@brief Magnitude refinement threshold (not currently used)
      */
      Real m_refMagn;
    };
  } // namespace BrownianWalker
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
