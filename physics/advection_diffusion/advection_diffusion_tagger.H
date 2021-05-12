/*!
  @file   advection_diffusion_tagger.H
  @brief  Class for tagging cells across an AMR hierarchy
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _ADVECTION_DIFFUSION_TAGGER_
#define _ADVECTION_DIFFUSION_TAGGER_

#include "cell_tagger.H"
#include "cdr_solver.H"
#include "CD_AmrMesh.H"

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace advection_diffusion {

    /*!
      @brief Abstract class for tagging cells across an AMR hierarchy
    */
    class advection_diffusion_tagger : public cell_tagger {
    public:

      advection_diffusion_tagger(RefCountedPtr<cdr_solver>& a_solver,
				 RefCountedPtr<AmrMesh>&   a_amr);
      ~advection_diffusion_tagger();

      void regrid();
      void parse_options();
      bool tag_cells(EBAMRTags& a_tags);

    protected:
      std::string m_realm;
      
      RefCountedPtr<cdr_solver> m_solver;
      RefCountedPtr<AmrMesh>   m_amr;

    private:
      advection_diffusion_tagger(){
	MayDay::Abort("advection_diffusion_tagger::advection_diffusion_tagger() - disabled constructor");
      }

      Real m_refi_curv;
      Real m_refi_magn;
    };
  }
}
#include "CD_NamespaceFooter.H"
#endif
