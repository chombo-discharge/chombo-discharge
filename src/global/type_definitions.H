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

typedef Vector<RefCountedPtr<LevelData<MFCellFAB> > >        MFAMRCellData;
typedef Vector<RefCountedPtr<LevelData<EBCellFAB> > >        EBAMRCellData;
typedef Vector<RefCountedPtr<LevelData<EBFluxFAB> > >        EBAMRFluxData;
typedef Vector<RefCountedPtr<LevelData<EBFaceFAB> > >        EBAMRFaceData;
typedef Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > EBAMRIVData;

#endif
