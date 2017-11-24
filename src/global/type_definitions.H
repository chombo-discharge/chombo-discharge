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

namespace relax {
  enum which_relax {
    jacobi       = 0,
    gauss_seidel = 1,
    gsrb_fast    = 2,
  };
};

namespace amrmg {
  enum which_mg {
    full   = 0,
    vcycle = 1,
    fcycle = 2,
  };
};

typedef Vector<RefCountedPtr<LevelData<MFCellFAB> > >        MFAMRCellData;
typedef Vector<RefCountedPtr<LevelData<EBCellFAB> > >        EBAMRCellData;
typedef Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        EBAMRFluxData;
typedef Vector<RefCountedPtr<LevelData<EBFaceFAB> > >        EBAMRFaceData;
typedef Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > EBAMRIVData;

#endif
