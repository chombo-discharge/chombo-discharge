/*!
  @file   type_definitions.H
  @brief  Typedefs, useful shortcuts
  @author Robert Marskar
  @date   Nov. 2017
*/

#ifndef _TYPE_DEFINITIONS_
#define _TYPE_DEFINITIONS_

#include <LevelData.H>
#include <RefCountedPtr.H>
#include <EBCellFAB.H>
#include <EBFluxFAB.H>
#include <EBFaceFAB.H>
#include <BaseIVFAB.H>
#include <BaseIFFAB.H>
#include <MFCellFAB.H>
#include <MFFluxFAB.H>
#include <Particle.H>
#include <ParticleData.H>
#include <ParticleValidRegion.H>
#include <BinItem.H>

#include "DomainFluxIFFAB.H"
#include "MFBaseIVFAB.H"
#include "photon.H"

/*!
  @brief Types of wall boundary conditions for Poisson equation
*/
namespace wallbc {
  enum which_bc {
    dirichlet     = 0,
    dirichlet_map = 1,
    neumann       = 2,
    neumann_map   = 3,
    robin         = 4,
    robin_map     = 5,
    outflow       = 6,
    symmetry      = 7,
    wall          = 8,
  };
};

#endif
