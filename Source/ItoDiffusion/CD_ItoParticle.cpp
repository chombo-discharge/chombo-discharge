/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoParticleImplem.cpp
  @brief  Implementation of CD_ItoParticle.H
  @author Robert Marskar
*/

// Our includes
#include <CD_ItoParticle.H>
#include <CD_NamespaceHeader.H>

std::vector<std::string> ItoParticle::s_realVariables = {"weight", "mobility", "diffusion", "energy", "tmpReal"};
std::vector<std::string> ItoParticle::s_vectVariables = {"oldPos", "velocity", "tmpVect"};

#include <CD_NamespaceFooter.H>
