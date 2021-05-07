/*!
  @file   brownian_walker_tagger.H
  @brief  Class for tagging cells across an AMR hierarchy
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _BROWNIAN_WALKER_TAGGER_
#define _BROWNIAN_WALKER_TAGGER_

#include "cell_tagger.H"
#include "ito_solver.H"
#include "amr_mesh.H"

namespace ChomboDischarge {
  namespace physics {
    namespace brownian_walker {

      /*!
	@brief Abstract class for tagging cells across an AMR hierarchy
      */
      class brownian_walker_tagger : public cell_tagger {
      public:

	brownian_walker_tagger(RefCountedPtr<ito_solver>& a_solver, RefCountedPtr<amr_mesh>& a_amr);
	~brownian_walker_tagger();

	void regrid() override;
	void parse_options() override;
	bool tag_cells(EBAMRTags& a_tags) override;

      protected:
	RefCountedPtr<ito_solver> m_solver;
	RefCountedPtr<amr_mesh>   m_amr;

      private:
	brownian_walker_tagger(){
	  MayDay::Abort("brownian_walker_tagger::brownian_walker_tagger() - disabled constructor");
	}

	Real m_refi_curv;
	Real m_refi_magn;
      };
    }
  }
}
#endif
