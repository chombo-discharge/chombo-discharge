/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BrownianWalkerSpecies.H
  @brief  Declaration of a Brownian walker species for the ItoSolver
  @author Robert Marskar
*/

#ifndef CD_BrownianWalkerSpecies_H
#define CD_BrownianWalkerSpecies_H

// Std includes
#include <random>

// Our includes
#include <CD_ItoSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace BrownianWalker {

    /*!
      @brief Species which is advected/diffused in Brownian walker code
      @details The initial particles are drawn from a Gaussian distribution. The user specifies how many, and the parameters for the distribution.
    */
    class BrownianWalkerSpecies : public ItoSpecies {
    public:

      /*!
	@brief Constructor -- reads in input parameters
      */
      BrownianWalkerSpecies();

      /*!
	@brief Destructor (does nothing)
      */
      ~BrownianWalkerSpecies();
      
    protected:

      /*!
	@brief Number of initial particles
      */
      int m_numParticles;      

      /*!
	@brief Radius for the initial Gaussian distribution of particles
      */
      Real m_blobRadius;

      /*!
	@brief Center-point for the Gaussian particle distribution
      */
      RealVect m_blobCenter;

      /*!
	@brief Draw initial particles
      */
      void drawInitParticles();
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
