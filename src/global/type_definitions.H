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
#include <MFCellFAB.H>
#include <MFFluxFAB.H>

#include "MFBaseIVFAB.H"

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

namespace ghost_interpolation {
  enum which_type {
    pwl  = 0,
    quad = 1,
  };
};

namespace centering {
  enum which_center {
    cell_center   = 0,
    cell_centroid = 1,
    eb_centroid   = 2,
  };
};

namespace time_code {
  enum which_code {
    cfl             = 0,
    diffusion       = 1,
    source          = 2,
    relaxation_time = 3,
    restricted      = 4,
    hardcap         = 5,
  };
};

namespace output_mode{
  enum which_mode{
    ultra_light = 0,
    light       = 1,
    full        = 2,
  };
};

typedef Vector<RefCountedPtr<LevelData<MFCellFAB> > >        MFAMRCellData;
typedef Vector<RefCountedPtr<LevelData<MFFluxFAB> > >        MFAMRFluxData;
typedef Vector<RefCountedPtr<LevelData<MFBaseIVFAB> > >      MFAMRIVData;
typedef Vector<RefCountedPtr<LevelData<EBCellFAB> > >        EBAMRCellData;
typedef Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        EBAMRFluxData;
typedef Vector<RefCountedPtr<LevelData<EBFaceFAB> > >        EBAMRFaceData;
typedef Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > EBAMRIVData;

#endif
