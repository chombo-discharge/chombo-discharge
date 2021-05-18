/*!
  @file   AdvectionDiffusionTagger.H
  @brief  Class for tagging cells across an AMR hierarchy
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _ADVECTION_DIFFUSION_TAGGER_
#define _ADVECTION_DIFFUSION_TAGGER_

#include "CD_CellTagger.H"
#include "CD_CdrSolver.H"
#include "CD_AmrMesh.H"

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace advection_diffusion {

    /*!
      @brief Abstract class for tagging cells across an AMR hierarchy
    */
    class AdvectionDiffusionTagger : public CellTagger {
    public:

      AdvectionDiffusionTagger(RefCountedPtr<CdrSolver>& a_solver,
				 RefCountedPtr<AmrMesh>&   a_amr);
      ~AdvectionDiffusionTagger();

      void regrid();
      void parseOptions();
      bool tagCells(EBAMRTags& a_tags);

    protected:
      std::string m_realm;
      
      RefCountedPtr<CdrSolver> m_solver;
      RefCountedPtr<AmrMesh>   m_amr;

    private:
      AdvectionDiffusionTagger(){
	MayDay::Abort("AdvectionDiffusionTagger::AdvectionDiffusionTagger() - disabled constructor");
      }

      Real m_refi_curv;
      Real m_refi_magn;
    };
  }
}
#include "CD_NamespaceFooter.H"
#endif
