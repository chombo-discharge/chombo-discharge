/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_AdvectionDiffusionTagger.H
  @brief  Tagging class for advection-diffusion physics
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
      @brief Abstract class for tagging cells across an AMR hierarchy
    */
    class AdvectionDiffusionTagger : public CellTagger {
    public:

      AdvectionDiffusionTagger(RefCountedPtr<CdrSolver>& a_solver,
				 RefCountedPtr<AmrMesh>&   a_amr);
      ~AdvectionDiffusionTagger();

      void regrid() override;
      void parseOptions() override;
      bool tagCells(EBAMRTags& a_tags) override;

    protected:
      std::string m_realm;
      
      RefCountedPtr<CdrSolver> m_solver;
      RefCountedPtr<AmrMesh>   m_amr;

    private:
      AdvectionDiffusionTagger() = delete;

      Real m_refi_curv;
      Real m_refi_magn;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
