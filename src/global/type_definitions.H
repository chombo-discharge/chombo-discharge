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

#include "EBAMRData.H"

/*!
  @brief phase names
*/
namespace phase {
  const int num_phases = 2;
  
  enum which_phase {
    gas   = 0,
    solid = 1
  };
};

/*!
  @brief Relaxation types for multigrid
*/
namespace relax {
  enum which_relax {
    jacobi       = 0,
    gauss_seidel = 1,
    gsrb_fast    = 2,
  };
};

/*!
  @brief Multigrid cycle types
*/
namespace amrmg {
  enum which_mg {
    full   = 0,
    vcycle = 1,
    fcycle = 2,
  };
};

/*!
  @brief Stencil types
*/
namespace stencil_type {
  enum which_type {
    linear = 0,
    taylor = 1,
    ray    = 2,
    lsq    = 3,
    pwl    = 4,
  };
};

/*!
  @brief Ghost cell interpolation methods
*/
namespace ghost_interpolation {
  enum which_type {
    pwl  = 0,
    quad = 1,
  };
};

/*!
  @brief Namespace for various centerings
*/
namespace centering {
  enum which_center {
    cell_center   = 0,
    cell_centroid = 1,
    eb_centroid   = 2,
  };
};



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

/*!
  @brief Species iterator. So that we can iterate over all cdr equations, or a subset of them
*/
namespace species_iteration {
  enum which_species {
    all      = 0,
    charged  = 1,
    negative = 2,
    positive = 3,
    neutral  = 4,
  };
};
    


//typedef Vector<RefCountedPtr<LayoutData<DenseIntVectSet> > > EBAMRTags;


#endif
