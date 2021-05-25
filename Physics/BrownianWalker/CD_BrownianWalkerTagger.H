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
      @brief Abstract class for tagging cells across an AMR hierarchy
    */
    class BrownianWalkerTagger : public CellTagger {
    public:

      BrownianWalkerTagger() = delete;

      BrownianWalkerTagger(RefCountedPtr<ItoSolver>& a_solver, RefCountedPtr<AmrMesh>& a_amr);
      ~BrownianWalkerTagger();

      void regrid() override;
      void parseOptions() override;
      bool tagCells(EBAMRTags& a_tags) override;

    protected:
      RefCountedPtr<ItoSolver> m_solver;
      RefCountedPtr<AmrMesh>   m_amr;

      Real m_refi_curv;
      Real m_refi_magn;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
