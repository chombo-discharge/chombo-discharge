/*!
  @file   photon_reaction.H
  @brief  Declaration of a class for holding photoionization reaction types in ito_plasma physics
  @author Robert Marskar
  @date   Aug. 2020
*/

#ifndef _PHOTO_REACTION_H_
#define _PHOTO_REACTION_H_

#include <Vector.H>
#include <RealVect.H>

namespace physics {
  namespace ito_plasma {
    
    /*!
      @brief Reaction class for ito_plasma_physics. 
      @details To use this class you should initialize a photoreaction e.g. photo_reaction({0}, {0,1,3})
      to the reaction Y0 => X_0 + X_1 + X_3.
      @note Technically, this class can also do Y0 + Y1 => X_0 + X1 + X3, but I don't see where that would be used...
    */
    class photo_reaction {
    public:

      /*!
	@brief Allow, but you need full ctor to do anything. 
      */
      photo_reaction() = default;

      /*!
	@brief Full Ctor, includes possible photons. 
      */
      photo_reaction(const int a_source,
		     const std::initializer_list<int> a_targets) : m_source(a_source),
								   m_targets(a_targets),
								   m_energies(0) {

      }

      /*!
	@brief Full Ctor, includes energy of generated species. 
      */
      photo_reaction(const int a_source,
		     const std::initializer_list<int> a_targets,
		     const std::initializer_list<std::pair<int, Real> > a_energies) : m_source(a_source),
										      m_targets(a_targets),
										      m_energies(a_energies){

      }

      /*!
	@brief Dtor. Not much to see here.
      */
      ~photo_reaction(){

      }

      /*!
	@brief Get the photon source
      */
      inline
      const int& source() const;

      /*!
	@brief Get the photon targets
      */
      inline
      const std::deque<int>& targets() const;

      /*!
	@brief Get the energies of the produced targets
      */
      inline
      const std::deque<std::pair<int, Real> >& get_target_energies() const;
      

    protected:

      int m_source;  // Source
      std::deque<int> m_targets; // Targets
      std::deque<std::pair<int, Real> > m_energies;
    };
  }
}

#include "photo_reactionI.H"

#endif
