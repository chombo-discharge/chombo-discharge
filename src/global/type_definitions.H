/*!
  @file type_definitions.H
  @brief Typedefs, useful shortcuts
  @author Robert Marskar
  @date Nov. 2017
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
#include "tags.H"
#include "MFBaseIVFAB.H"
#include "photon.H"

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
  @brief Time code for output
*/
namespace time_code {
  enum which_code {
    cfl             = 0,
    diffusion       = 1,
    source          = 2,
    relaxation_time = 3,
    restricted      = 4,
    hardcap         = 5,
    error           = 6,
  };
};

/*!
  @brief Output mode for plasma_engine
*/
namespace output_mode{
  enum which_mode{
    ultra_light = 0,
    light       = 1,
    medium      = 2,
    full        = 3,
  };
};

/*!
  @brief Potential types on domain edges (faces). 
*/
namespace potential {
  enum ground_live {
    ground = 0,
    live   = 1,
    mapped = 2, // Not currently supported. But it will. 
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
  };
};

/*!
  @brief Boundary condition types for the CDR solvers
*/
namespace cdr_bc {
  enum which_bc {
    external = 0,
    wall     = 1,
    outflow  = 2,
    extrap   = 3,
  };
};
    

/*!
  @brief Restart modes
*/
namespace restart_mode {
  enum which_mode {
    full                = 0,
    surface_charge_only = 1,
    volume_charge_only  = 2,
  };
};

namespace memory_report_mode {
  enum which_mode {
    overall   = 0,
    unfreed   = 1,
    allocated = 2,
  };
};

typedef Vector<RefCountedPtr<LevelData<MFCellFAB> > >        MFAMRCellData;
typedef Vector<RefCountedPtr<LevelData<MFFluxFAB> > >        MFAMRFluxData;
typedef Vector<RefCountedPtr<LevelData<MFBaseIVFAB> > >      MFAMRIVData;
typedef Vector<RefCountedPtr<LevelData<EBCellFAB> > >        EBAMRCellData;
typedef Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        EBAMRFluxData;
typedef Vector<RefCountedPtr<LevelData<EBFaceFAB> > >        EBAMRFaceData;
typedef Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > EBAMRIVData;
typedef Vector<RefCountedPtr<LevelData<DomainFluxIFFAB > > > EBAMRIFData;
typedef Vector<RefCountedPtr<LevelData<tags> > >             EBAMRTags;
typedef Vector<RefCountedPtr<ParticleData<Particle> > >      EBAMRParticles;
typedef Vector<RefCountedPtr<ParticleData<photon> > >        EBAMRPhotons;
typedef Vector<RefCountedPtr<ParticleValidRegion> >          EBAMRPVR;

#endif
