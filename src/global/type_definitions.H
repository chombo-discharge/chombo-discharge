/*!
  @file type_definitions.H
  @brief Typedefs, useful shortcuts
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _TYPE_DEFINITIONS_
#define _TYPE_DEFINITIONS_

#include <LevelData.H>
#include <EBCellFAB.H>
#include <EBFluxFAB.H>
#include <EBFaceFAB.H>
#include <BaseIVFAB.H>
#include <MFCellFAB.H>

typedef Vector<LevelData<MFCellFAB>* >        MFAMRCellData;
typedef Vector<LevelData<EBCellFAB>* >        EBAMRCellData;
typedef Vector<LevelData<EBFluxFAB>* >        EBAMRFluxData;
typedef Vector<LevelData<EBFaceFAB>* >        EBAMRFaceData;
typedef Vector<LevelData<BaseIVFAB<Real> >* > EBAMRIVData;

#endif
