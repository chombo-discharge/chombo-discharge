/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_AmrMesh.cpp
  @brief  Implementation of CD_AmrMesh.H
  @author Robert Marskar
*/

// Chombo includes
#include <BRMeshRefine.H>
#include <EBEllipticLoadBalance.H>
#include <EBLevelDataOps.H>
#include <MFLevelDataOps.H>
#include <EBArith.H>
#include <ParmParse.H>
#include <BaseIFFactory.H>
#include <BaseIVFactory.H>

// Our includes
#include <CD_AmrMesh.H>
#include <CD_MultifluidAlias.H>
#include <CD_LoadBalancing.H>
#include <CD_Timer.H>
#include <CD_DomainFluxIFFABFactory.H>
#include <CD_TiledMeshRefine.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

AmrMesh::AmrMesh()
{

  // Default things
  this->parseOptions();

  m_finestLevel = 0;
  m_hasGrids    = false;

  // Some things might require a vector which is just a tiny bit longer.
  m_refinementRatios.resize(m_maxAmrDepth);
  m_refinementRatios.push_back(2);
}

AmrMesh::~AmrMesh() {}

EBAMRCellData
AmrMesh::alias(const phase::which_phase a_phase, const MFAMRCellData& a_mfdata) const
{
  CH_TIME("AmrMesh::alias(phase::which_phase, MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::alias(phase::which_phase, MFAMRCellData)" << endl;
  }

  EBAMRCellData ret;

  const int finestLevel = a_mfdata.size() - 1;

  this->allocatePointer(ret, finestLevel);

  this->alias(ret, a_phase, a_mfdata, finestLevel);

  ret.setRealm(a_mfdata.getRealm());

  return ret;
}

EBAMRFluxData
AmrMesh::alias(const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata) const
{
  CH_TIME("AmrMesh::alias(phase::which_phase, MFAMRFluxData)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::alias(phase::which_phase, MFAMRFluxData)" << endl;
  }

  EBAMRFluxData ret;

  const int finestLevel = a_mfdata.size() - 1;

  this->allocatePointer(ret, finestLevel);

  this->alias(ret, a_phase, a_mfdata, finestLevel);

  ret.setRealm(a_mfdata.getRealm());

  return ret;
}

EBAMRIVData
AmrMesh::alias(const phase::which_phase a_phase, const MFAMRIVData& a_mfdata) const
{
  CH_TIME("AmrMesh::alias(phase, mfamrivdata)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::alias(phase, mfamrivdata)" << endl;
  }

  EBAMRIVData ret;

  this->allocatePointer(ret);

  this->alias(ret, a_phase, a_mfdata);

  ret.setRealm(a_mfdata.getRealm());

  return ret;
}

void
AmrMesh::alias(EBAMRCellData&           a_data,
               const phase::which_phase a_phase,
               const MFAMRCellData&     a_mfdata,
               const int                a_finestLevel) const
{
  CH_TIME("AmrMesh::alias(EBAMRCellData, phase::which_phase, MFAMRCellData, int))");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::alias(EBAMRCellData, phase::which_phase, MFAMRCellData, int)" << endl;
  }

  for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
    MultifluidAlias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void
AmrMesh::alias(EBAMRFluxData&           a_data,
               const phase::which_phase a_phase,
               const MFAMRFluxData&     a_mfdata,
               const int                a_finestLevel) const
{
  CH_TIME("AmrMesh::alias(EBAMRFluxData, phase::which_phase, EBAMRFluxData, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::alias(hardcap)" << endl;
  }

  for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
    MultifluidAlias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void
AmrMesh::alias(EBAMRCellData& a_data, const phase::which_phase a_phase, const MFAMRCellData& a_mfdata) const
{
  CH_TIME("AmrMesh::alias(EBAMRCellData, phase::which_phase, MFAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::alias(EBAMRCellData, phase::which_phase, MFAMRCellData)" << endl;
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    MultifluidAlias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void
AmrMesh::alias(EBAMRFluxData& a_data, const phase::which_phase a_phase, const MFAMRFluxData& a_mfdata) const
{
  CH_TIME("AmrMesh::alias(EBAMRFluxData, phase::which_phase, MFAMRFluxData)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::alias(EBAMRFluxData, phase::which_phase, MFAMRFluxData)" << endl;
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    MultifluidAlias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void
AmrMesh::alias(EBAMRIVData& a_data, const phase::which_phase a_phase, const MFAMRIVData& a_mfdata) const
{
  CH_TIME("AmrMesh::alias(EBAMRIVData, phase::which_phase, MFAMRIVData)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::alias(EBAMRIVData, phase::which_phase, MFAMRIVData)" << endl;
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    MultifluidAlias::aliasMF(*a_data[lvl], a_phase, *a_mfdata[lvl]);
  }
}

void
AmrMesh::allocate(EBAMRCellData&           a_data,
                  const std::string        a_realm,
                  const phase::which_phase a_phase,
                  const int                a_nComp,
                  const int                a_ghost) const
{
  CH_TIME("AmrMesh::allocate(EBAMRCellData, string, phase::which_phase, int, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(EBAMRCellData, string, phase::which_phase, int, int)" << endl;
  }

  CH_assert(a_nComp > 0);

  // This allocates data on a specific realm and phase with specified number of components and ghost cells. If a_ghost < 0
  // we use the default number of ghost cells in AmrMesh.

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::allocate(EBAMRCellData, string, phase::which_phase, int, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost < 0) ? m_numGhostCells : a_ghost;

  a_data.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout&        ebisl = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];

    a_data[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(
      new LevelData<EBCellFAB>(dbl, a_nComp, ghost * IntVect::Unit, EBCellFactory(ebisl)));
  }

  a_data.setRealm(a_realm);
}

void
AmrMesh::allocate(EBAMRFluxData&           a_data,
                  const std::string        a_realm,
                  const phase::which_phase a_phase,
                  const int                a_nComp,
                  const int                a_ghost) const
{
  CH_TIME("AmrMesh::allocate(EBAMRFluxData, string, phase::which_phase, int, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(EBAMRFluxData, string, phase::which_phase, int, int)" << endl;
  }

  CH_assert(a_nComp > 0);

  // This allocates data on a specific realm and phase with specified number of components and ghost cells. If a_ghost < 0
  // we use the default number of ghost cells in AmrMesh.

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::allocate(EBAMRFluxData, string, phase::which_phase, int, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost < 0) ? m_numGhostCells : a_ghost;

  a_data.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout&        ebisl = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];

    a_data[lvl] = RefCountedPtr<LevelData<EBFluxFAB>>(
      new LevelData<EBFluxFAB>(dbl, a_nComp, ghost * IntVect::Unit, EBFluxFactory(ebisl)));
  }

  a_data.setRealm(a_realm);
}

void
AmrMesh::allocate(EBAMRIVData&             a_data,
                  const std::string        a_realm,
                  const phase::which_phase a_phase,
                  const int                a_nComp,
                  const int                a_ghost) const
{
  CH_TIME("AmrMesh::allocate(EBAMRIVData, string, phase::which_phase, int, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(EBAMRIVData, string, phase::which_phase, int, int)" << endl;
  }

  CH_assert(a_nComp > 0);

  // This allocates data on a specific realm and phase with specified number of components and ghost cells. If a_ghost < 0
  // we use the default number of ghost cells in AmrMesh.

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::allocate(EBAMRIVData, string, phase::which_phase, int, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost < 0) ? m_numGhostCells : a_ghost;

  a_data.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl    = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout&        ebisl  = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];
    const ProblemDomain&     domain = m_realms[a_realm]->getDomains()[lvl];

    LayoutData<IntVectSet> irregCells(dbl);
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      Box box = dbl.get(dit());
      box.grow(ghost);
      box &= domain;

      irregCells[dit()] = ebisl[dit()].getIrregIVS(box);
    }

    a_data[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
      new LevelData<BaseIVFAB<Real>>(dbl, a_nComp, ghost * IntVect::Unit, BaseIVFactory<Real>(ebisl, irregCells)));
  }

  a_data.setRealm(a_realm);
}

void
AmrMesh::allocate(EBAMRIFData&             a_data,
                  const std::string        a_realm,
                  const phase::which_phase a_phase,
                  const int                a_nComp,
                  const int                a_ghost) const
{
  CH_TIME("AmrMesh::allocate(EBAMRIFData, string, phase::which_phase, int, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(EBAMRIFData, string, phase::which_phase, int, int)" << endl;
  }

  CH_assert(a_nComp > 0);

  // This allocates data on a specific realm and phase with specified number of components and ghost cells. If a_ghost < 0
  // we use the default number of ghost cells in AmrMesh.

  if (!this->queryRealm(a_realm)) {
    std::string str =
      "AmrMesh::allocate(EBAMRIFData, string, phase::which_phase, int, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost < 0) ? m_numGhostCells : a_ghost;

  a_data.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl    = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout&        ebisl  = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];
    const ProblemDomain&     domain = m_realms[a_realm]->getDomains()[lvl];

    a_data[lvl] = RefCountedPtr<LevelData<DomainFluxIFFAB>>(
      new LevelData<DomainFluxIFFAB>(dbl, a_nComp, ghost * IntVect::Unit, DomainFluxIFFABFactory(ebisl, domain)));
  }

  a_data.setRealm(a_realm);
}

void
AmrMesh::allocate(EBAMRBool& a_data, const std::string a_realm, const int a_nComp, const int a_ghost) const
{
  CH_TIME("AmrMesh::allocate(EBAMRBool, string, int, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(EBAMRBool, string, int, int)" << endl;
  }

  CH_assert(a_nComp > 0);

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(EBAMRBool, string, int, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost = (a_ghost < 0) ? m_numGhostCells : a_ghost;

  a_data.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];

    a_data[lvl] =
      RefCountedPtr<LevelData<BaseFab<bool>>>(new LevelData<BaseFab<bool>>(dbl, a_nComp, ghost * IntVect::Unit));
  }

  a_data.setRealm(a_realm);
}

void
AmrMesh::allocate(MFAMRCellData& a_data, const std::string a_realm, const int a_nComp, const int a_ghost) const
{
  CH_TIME("AmrMesh::allocate(MFAMRCellData, string, int, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(MFAMRCellData, string, int, int)" << endl;
  }

  CH_assert(a_nComp > 0);

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(MFAMRCellData, string, int, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost   = (a_ghost < 0) ? m_numGhostCells : a_ghost;
  const int ignored = a_nComp; // Strange but true thing -- number of components come in through the factory.
  const int nphases = m_multifluidIndexSpace->numPhases();

  a_data.resize(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace> ebisGas = m_realms[a_realm]->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace> ebisSol = m_realms[a_realm]->getEBIndexSpace(phase::solid);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];

    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_nComp);

    if (!ebisGas.isNull())
      ebisl[phase::gas] = m_realms[a_realm]->getEBISLayout(phase::gas)[lvl];
    if (!ebisSol.isNull())
      ebisl[phase::solid] = m_realms[a_realm]->getEBISLayout(phase::solid)[lvl];

    MFCellFactory factory(ebisl, comps);

    a_data[lvl] =
      RefCountedPtr<LevelData<MFCellFAB>>(new LevelData<MFCellFAB>(dbl, ignored, ghost * IntVect::Unit, factory));
  }

  a_data.setRealm(a_realm);
}

void
AmrMesh::allocate(MFAMRFluxData& a_data, const std::string a_realm, const int a_nComp, const int a_ghost) const
{
  CH_TIME("AmrMesh::allocate(MFAMRFluxData, string, int, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(MFAMRFluxData, string, int, int)" << endl;
  }

  CH_assert(a_nComp > 0);

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(MFAMRFluxData, string, int, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost   = (a_ghost < 0) ? m_numGhostCells : a_ghost;
  const int ignored = a_nComp; // Strange but true thing -- number of components come in through the factory.
  const int nphases = m_multifluidIndexSpace->numPhases();

  a_data.resize(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_realms[a_realm]->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_realms[a_realm]->getEBIndexSpace(phase::solid);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];

    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_nComp);

    if (!ebisGas.isNull())
      ebisl[phase::gas] = m_realms[a_realm]->getEBISLayout(phase::gas)[lvl];
    if (!ebisSol.isNull())
      ebisl[phase::solid] = m_realms[a_realm]->getEBISLayout(phase::solid)[lvl];

    MFFluxFactory factory(ebisl, comps);

    a_data[lvl] =
      RefCountedPtr<LevelData<MFFluxFAB>>(new LevelData<MFFluxFAB>(dbl, ignored, ghost * IntVect::Unit, factory));
  }

  a_data.setRealm(a_realm);
}

void
AmrMesh::allocate(MFAMRIVData& a_data, const std::string a_realm, const int a_nComp, const int a_ghost) const
{
  CH_TIME("AmrMesh::allocate(MFAMRIVData, string, int, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(MFAMRIVData, string, int, int)" << endl;
  }

  CH_assert(a_nComp > 0);

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::allocate(MFAMRIVData, string, int, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const int ghost   = (a_ghost < 0) ? m_numGhostCells : a_ghost;
  const int ignored = a_nComp;
  const int nphases = m_multifluidIndexSpace->numPhases();

  a_data.resize(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_realms[a_realm]->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_realms[a_realm]->getEBIndexSpace(phase::solid);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];

    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, a_nComp);

    if (!ebisGas.isNull())
      ebisl[phase::gas] = m_realms[a_realm]->getEBISLayout(phase::gas)[lvl];
    if (!ebisSol.isNull())
      ebisl[phase::solid] = m_realms[a_realm]->getEBISLayout(phase::solid)[lvl];

    MFBaseIVFABFactory factory(ebisl, comps);

    a_data[lvl] =
      RefCountedPtr<LevelData<MFBaseIVFAB>>(new LevelData<MFBaseIVFAB>(dbl, ignored, ghost * IntVect::Unit, factory));
  }

  a_data.setRealm(a_realm);
}

void
AmrMesh::reallocate(EBAMRCellData& a_data, const phase::which_phase a_phase, const int a_lmin) const
{
  CH_TIME("AmrMesh::reallocate(EBAMRCellData, phase::which_phase, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::reallocate(EBAMRCellData, phase::which_phase, int)" << endl;
  }

  CH_assert(a_lmin >= 0);

  // TLDR: This reallocates data from a_lmin on the same realm as before. Phase needs to be specified because EBAMRCellData does not know about phases.
  //       Data under a_lmin is not touched.

  const std::string a_realm = a_data.getRealm();

  if (!this->queryRealm(a_realm)) {
    std::string str =
      "AmrMesh::reallocate(EBAMRCellData, phase::which_phase, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int     nComp = a_data[0]->nComp();

  a_data.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout&        ebisl = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];

    a_data[lvl] =
      RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>(dbl, nComp, ghost, EBCellFactory(ebisl)));
  }
}

void
AmrMesh::reallocate(EBAMRFluxData& a_data, const phase::which_phase a_phase, const int a_lmin) const
{
  CH_TIME("AmrMesh::reallocate(EBAMRFluxData, phase::which_phase, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::reallocate(EBAMRFluxData, phase::which_phase, int)" << endl;
  }

  CH_assert(a_lmin >= 0);

  const std::string a_realm = a_data.getRealm();

  if (!this->queryRealm(a_realm)) {
    std::string str =
      "AmrMesh::reallocate(EBAMRFluxData, phase::which_phase, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int     nComp = a_data[0]->nComp();

  a_data.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl   = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout&        ebisl = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];

    a_data[lvl] =
      RefCountedPtr<LevelData<EBFluxFAB>>(new LevelData<EBFluxFAB>(dbl, nComp, ghost, EBFluxFactory(ebisl)));
  }
}

void
AmrMesh::reallocate(EBAMRIVData& a_data, const phase::which_phase a_phase, const int a_lmin) const
{
  CH_TIME("AmrMesh::reallocate(EBAMRIVData, phase::which_phase, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::reallocate(EBAMRIVData, phase::which_phase, int)" << endl;
  }

  CH_assert(a_lmin >= 0);

  const std::string a_realm = a_data.getRealm();

  if (!this->queryRealm(a_realm)) {
    std::string str =
      "AmrMesh::reallocate(EBAMRIVData, phase::which_phase, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int     nComp = a_data[0]->nComp();

  a_data.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl    = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout&        ebisl  = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];
    const ProblemDomain&     domain = m_realms[a_realm]->getDomains()[lvl];

    LayoutData<IntVectSet> irregCells(dbl);
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      Box box = dbl.get(dit());
      box.grow(ghost);
      box &= domain;

      irregCells[dit()] = ebisl[dit()].getIrregIVS(box);
    }

    a_data[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
      new LevelData<BaseIVFAB<Real>>(dbl, nComp, ghost, BaseIVFactory<Real>(ebisl, irregCells)));
  }
}

void
AmrMesh::reallocate(EBAMRIFData& a_data, const phase::which_phase a_phase, const int a_lmin) const
{
  CH_TIME("AmrMesh::reallocate(EBAMRIFData, phase::which_phase, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::reallocate(EBAMRIFData, phase::which_phase, int)" << endl;
  }

  CH_assert(a_lmin >= 0);

  const std::string a_realm = a_data.getRealm();

  if (!this->queryRealm(a_realm)) {
    std::string str =
      "AmrMesh::reallocate(EBAMRIFData, phase::which_phase, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int     nComp = a_data[0]->nComp();

  a_data.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl    = m_realms[a_realm]->getGrids()[lvl];
    const EBISLayout&        ebisl  = m_realms[a_realm]->getEBISLayout(a_phase)[lvl];
    const ProblemDomain&     domain = m_realms[a_realm]->getDomains()[lvl];

    a_data[lvl] = RefCountedPtr<LevelData<DomainFluxIFFAB>>(
      new LevelData<DomainFluxIFFAB>(dbl, nComp, ghost, DomainFluxIFFABFactory(ebisl, domain)));
  }
}

void
AmrMesh::reallocate(EBAMRBool& a_data, const int a_lmin) const
{
  CH_TIME("AmrMesh::reallocate(EBAMRBool, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::reallocate(EBAMRBool, int)" << endl;
  }

  CH_assert(a_lmin >= 0);

  const std::string a_realm = a_data.getRealm();

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::reallocate(EBAMRBool, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int     nComp = a_data[0]->nComp();

  a_data.resize(1 + m_finestLevel);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];

    a_data[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(new LevelData<BaseFab<bool>>(dbl, nComp, ghost));
  }
}

void
AmrMesh::reallocate(MFAMRCellData& a_data, const int a_lmin) const
{
  CH_TIME("AmrMesh::reallocate(MFAMRCellData, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::reallocate(MFAMRCellData, int)" << endl;
  }

  CH_assert(a_lmin >= 0);

  const std::string a_realm = a_data.getRealm();

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::reallocate(MFAMRCellData,int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int     nComp = a_data[0]->nComp();
  const int     ignored =
    nComp; // A strange but true thing -- for multifluid data we pass in the number of components through the factory.
  const int nphases = m_multifluidIndexSpace->numPhases();

  a_data.resize(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_realms[a_realm]->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_realms[a_realm]->getEBIndexSpace(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, nComp);

    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];

    if (!ebisGas.isNull())
      ebisl[phase::gas] = m_realms[a_realm]->getEBISLayout(phase::gas)[lvl];
    if (!ebisSol.isNull())
      ebisl[phase::solid] = m_realms[a_realm]->getEBISLayout(phase::solid)[lvl];

    MFCellFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFCellFAB>>(new LevelData<MFCellFAB>(dbl, ignored, ghost, factory));
  }
}

void
AmrMesh::reallocate(MFAMRFluxData& a_data, const int a_lmin) const
{
  CH_TIME("AmrMesh::allocate(MFAMRFluxData, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(MFAMRFluxData, int)" << endl;
  }

  CH_assert(a_lmin >= 0);

  const std::string a_realm = a_data.getRealm();

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::reallocate(MFAMRFluxData, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int     nComp = a_data[0]->nComp();
  const int     ignored =
    nComp; // Strange but true thing, for multifluid data the number of components come in through the factory.
  const int nphases = m_multifluidIndexSpace->numPhases();

  a_data.resize(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_realms[a_realm]->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_realms[a_realm]->getEBIndexSpace(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, nComp);

    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];

    if (!ebisGas.isNull())
      ebisl[phase::gas] = m_realms[a_realm]->getEBISLayout(phase::gas)[lvl];
    if (!ebisSol.isNull())
      ebisl[phase::solid] = m_realms[a_realm]->getEBISLayout(phase::solid)[lvl];

    MFFluxFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFFluxFAB>>(new LevelData<MFFluxFAB>(dbl, ignored, ghost, factory));
  }
}

void
AmrMesh::reallocate(MFAMRIVData& a_data, const int a_lmin) const
{
  CH_TIME("AmrMesh::allocate(MFAMRIVData, int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::allocate(MFAMRIVData, int)" << endl;
  }

  CH_assert(a_lmin >= 0);

  const std::string a_realm = a_data.getRealm();

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::reallocate(MFAMRIVData, int) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IntVect ghost = a_data[0]->ghostVect();
  const int     nComp = a_data[0]->nComp();
  const int     ignored =
    nComp; // Strange but true thing, for multifluid data the number of components come in through the factory.
  const int nphases = m_multifluidIndexSpace->numPhases();

  a_data.resize(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_realms[a_realm]->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_realms[a_realm]->getEBIndexSpace(phase::solid);

  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
    Vector<EBISLayout> ebisl(nphases);
    Vector<int>        comps(nphases, nComp);

    const DisjointBoxLayout& dbl = m_realms[a_realm]->getGrids()[lvl];

    if (!ebisGas.isNull())
      ebisl[phase::gas] = m_realms[a_realm]->getEBISLayout(phase::gas)[lvl];
    if (!ebisSol.isNull())
      ebisl[phase::solid] = m_realms[a_realm]->getEBISLayout(phase::solid)[lvl];

    MFBaseIVFABFactory factory(ebisl, comps);

    a_data[lvl] = RefCountedPtr<LevelData<MFBaseIVFAB>>(new LevelData<MFBaseIVFAB>(dbl, ignored, ghost, factory));
  }
}

void
AmrMesh::setMultifluidIndexSpace(const RefCountedPtr<MultiFluidIndexSpace>& a_multiFluidIndexSpace)
{
  CH_TIME("AmrMesh::setMultifluidIndexSpace(RefCountedPtr<MultifluidIndexSpace>");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::setMultifluidIndexSpace(RefCountedPtr<MultifluidIndexSpace>)" << endl;
  }

  CH_assert(!a_multiFluidIndexSpace.isNull());

  m_multifluidIndexSpace = a_multiFluidIndexSpace;
}

void
AmrMesh::setBaseImplicitFunction(const phase::which_phase a_phase, const RefCountedPtr<BaseIF>& a_baseIF)
{
  CH_TIME("AmrMesh::setBaseImplicitFunction(phase::which_phase, RefCountedPtr<BaseIF>)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::setBaseImplicitFunction(phase::which_phase, RefCountedPtr<BaseIF>)" << endl;
  }

  m_baseif.emplace(a_phase, a_baseIF);
}

void
AmrMesh::parseOptions()
{
  CH_TIME("AmrMesh::parseOptions()");

  this->parseVerbosity();

  if (m_verbosity > 5) {
    pout() << "AmrMesh::parseOptions()" << endl;
  }

  this->parseCoarsestLevelNumCells();
  this->parseMaxAmrDepth();
  this->parseMaxSimulationDepth();
  this->parseRefinementRatios();
  this->parseBlockingFactor();
  this->parseMaxBoxSize();
  this->parseMaxEbisBoxSize();
  this->parseGridGeneration();
  this->parseBrBufferSize();
  this->parseBrFillRatio();
  this->parseRedistributionRadius();
  ;
  this->parseNumGhostCells();
  this->parseEbGhostCells();
  this->parseProbLoHiCorners();
  this->parseCentroidStencils();
  this->parseEbCentroidStencils();
  this->parseMultigridInterpolator();

  this->sanityCheck();
  this->buildDomains();
}

void
AmrMesh::parseRuntimeOptions()
{
  CH_TIME("AmrMesh::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::parseRuntimeOptions()" << endl;
  }

  this->parseMaxSimulationDepth();
  this->parseVerbosity();
  this->parseBlockingFactor();
  this->parseMaxBoxSize();
  this->parseGridGeneration();
  this->parseBrBufferSize();
  this->parseBrFillRatio();
  this->parseMultigridInterpolator();
}

void
AmrMesh::parseProbLoHiCorners()
{
  CH_TIME("AmrMesh::parseProbLoHiCorners()");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::parseProbLoHiCorners()" << endl;
  }

  ParmParse pp("AmrMesh");

  Vector<Real> v(SpaceDim);

  pp.getarr("lo_corner", v, 0, SpaceDim);
  m_probLo = RealVect(D_DECL(v[0], v[1], v[2]));
  pp.getarr("hi_corner", v, 0, SpaceDim);
  m_probHi = RealVect(D_DECL(v[0], v[1], v[2]));
}

void
AmrMesh::buildDomains()
{
  CH_TIME("AmrMesh::buildDomains()");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::buildDomains()" << endl;
  }

  const int numLevels = 1 + m_maxAmrDepth;

  m_domains.resize(numLevels);
  m_dx.resize(numLevels);
  m_grids.resize(numLevels);

  m_dx[0]      = (m_probHi[0] - m_probLo[0]) / m_numCells[0];
  m_domains[0] = ProblemDomain(IntVect::Zero, m_numCells - IntVect::Unit);

  for (int lvl = 1; lvl <= m_maxAmrDepth; lvl++) {
    m_dx[lvl]      = m_dx[lvl - 1] / m_refinementRatios[lvl - 1];
    m_domains[lvl] = m_domains[lvl - 1];

    m_domains[lvl].refine(m_refinementRatios[lvl - 1]);
  }
}

void
AmrMesh::regridAmr(const Vector<IntVectSet>& a_tags, const int a_lmin, const int a_hardcap)
{
  CH_TIME("AmrMesh::regridAmr(Vector<IntVectSet>, int, int)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::regridAmr(Vector<IntVectSet>, int, int)" << endl;
  }

  CH_assert(a_lmin >= 0);

  // TLDR: This is the version that reads boxes. AmrMesh makes the grids from the tags and load balances them
  //       by using the patch volume. Those grids are then sent to the various Realms.

  Vector<IntVectSet> tags = a_tags; // buildGrids destroys tags, so we actually have to copy them.

  this->buildGrids(tags, a_lmin, a_hardcap); // Build AMR grids -- the realm grids are defined with these grids.
  this->defineRealms();                      // Define Realms with the new grids and redo the Realm stuff

  for (auto& r : m_realms) {
    r.second->regridBase(a_lmin);
  }
}

void
AmrMesh::regridOperators(const int a_lmin)
{
  CH_TIME("AmrMesh::regridOperators(int)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::regridOperators(int)" << endl;
  }

  for (auto& r : m_realms) {
    r.second->regridOperators(a_lmin);
  }
}

void
AmrMesh::buildGrids(Vector<IntVectSet>& a_tags, const int a_lmin, const int a_hardcap)
{
  CH_TIME("AmrMesh::buildGrids(Vector<IntVectSet>, int, int)");
  if (m_verbosity > 2) {
    pout() << "AmrMesh::buildGrids(Vector<IntVectSet>, int, int)" << endl;
  }

  // TLDR: a_lmin is the coarsest level that changes and a_hardcap is a hardcap for the maximum grid level that can
  //       be generated. If a_hardcap < 0 the restriction is m_maxAmrDepth.

  // baseLevel is the coarsest level which does not change. topLevel is the finest level where we have tags. We should never
  // have tags on max_amr_depth, and we make that restriction here.
  const int baseLevel = std::max(0, a_lmin - 1);
  const int topLevel  = (m_finestLevel == m_maxAmrDepth) ? m_finestLevel - 1 : a_tags.size() - 1;

  // New and old grid boxes
  Vector<Vector<Box>> newBoxes(1 + topLevel);
  Vector<Vector<Box>> oldBoxes(1 + topLevel);

  // Enforce potential hardcap.
  const int hardcap = (a_hardcap < 0) ? m_maxAmrDepth : a_hardcap;

  // Inside this loop we make the boxes.
  if (m_maxAmrDepth > 0 && hardcap > 0) {
    domainSplit(m_domains[0], oldBoxes[0], m_maxBoxSize, m_blockingFactor);

    // If have old grids, we can use the old boxes as input to the grid generators (since they don't necessarily regrid all levels).
    if (!m_hasGrids) {
      for (int lvl = 1; lvl <= topLevel; lvl++) {
        oldBoxes[lvl].resize(0);
        oldBoxes[lvl].push_back(m_domains[lvl].domainBox());
      }
    }
    else {
      for (int lvl = 0; lvl <= topLevel; lvl++) {
        oldBoxes[lvl] = m_grids[lvl].boxArray();
      }
    }

    // Berger-Rigoutsos grid generation
    int newFinestLevel;

    switch (m_gridGenerationMethod) {
    case GridGenerationMethod::BergerRigoutsous: {
      BRMeshRefine meshRefine(m_domains[0],
                              m_refinementRatios,
                              m_fillRatioBR,
                              m_blockingFactor,
                              m_bufferSizeBR,
                              m_maxBoxSize);
      newFinestLevel = meshRefine.regrid(newBoxes, a_tags, baseLevel, topLevel, oldBoxes);
      break;
    }
    case GridGenerationMethod::Tiled: {
      TiledMeshRefine meshRefine(m_domains[0], m_refinementRatios, m_blockingFactor * IntVect::Unit);
      newFinestLevel = meshRefine.regrid(newBoxes, a_tags, baseLevel, topLevel, oldBoxes);
      break;
    }
    default:
      MayDay::Error(
        "AmrMesh::buildGrids(Vector<IntVectSet>, int, int) - logic bust, regridding with unknown regrid algorithm");
      break;
    }

    // Identify the new finest grid level.
    m_finestLevel = Min(newFinestLevel, m_maxAmrDepth);       // Don't exceed m_maxAmrDepth
    m_finestLevel = Min(m_finestLevel, m_maxSimulationDepth); // Don't exceed maximum simulation depth
    m_finestLevel = Min(m_finestLevel, hardcap);              // Don't exceed hardcap
  }
  else { // Only end up here if we have a single grid level, i.e. just single-level grid decomposition.
    newBoxes.resize(1);
    domainSplit(m_domains[0], newBoxes[0], m_maxBoxSize, m_blockingFactor);

    m_finestLevel = 0;
  }

  // Coarsest level also changes in this case, but that's not actually caught by the regridders. We have to do this because the blocking
  // factor could have been changed during runtime option parsing.
  if (a_lmin == 0) {
    domainSplit(m_domains[0], newBoxes[0], m_maxBoxSize, m_blockingFactor);
  }

  // Sort the boxes and then load balance them, using the patch volume as an initial proxy.
  Vector<Vector<int>> processorIDs(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    LoadBalancing::sort(newBoxes[lvl], m_boxSort);
    LoadBalancing::makeBalance(processorIDs[lvl], newBoxes[lvl]);
  }

  // Now we define the grids. If a_lmin=0 every grid is new, otherwise keep old grids up to but not including a_lmin
  if (a_lmin == 0) {
    m_grids.resize(1 + m_finestLevel);
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      m_grids[lvl] = DisjointBoxLayout();
      m_grids[lvl].define(newBoxes[lvl], processorIDs[lvl], m_domains[lvl]);
      m_grids[lvl].close();
    }
  }
  else {
    Vector<DisjointBoxLayout> old_grids = m_grids;
    m_grids.resize(1 + m_finestLevel);
    for (int lvl = 0; lvl < a_lmin; lvl++) { // Copy old grids
      m_grids[lvl] = old_grids[lvl];
    }
    for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) { // Create new ones from tags
      m_grids[lvl] = DisjointBoxLayout();
      m_grids[lvl].define(newBoxes[lvl], processorIDs[lvl], m_domains[lvl]);
      m_grids[lvl].close();
    }
  }

  m_hasGrids = true;
}

void
AmrMesh::computeGradient(LevelData<EBCellFAB>&       a_gradient,
                         const LevelData<EBCellFAB>& a_phi,
                         const std::string           a_realm,
                         const phase::which_phase    a_phase,
                         const int                   a_lvl) const
{
  CH_TIME("AmrMesh::computeGradient(LD<EBCellFAB>, LD<EBCellFAB>, string, phase::which_phase,int)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::computeGradient(LD<EBCellFAB>, LD<EBCellFAB>, string, phase::which_phase,int)" << endl;
  }

  const RefCountedPtr<EBGradient>& gradientOp = m_realms[a_realm]->getGradientOp(a_phase)[a_lvl];

  gradientOp->computeLevelGradient(a_gradient, a_phi);
}

void
AmrMesh::computeGradient(EBAMRCellData&           a_gradient,
                         const EBAMRCellData&     a_phi,
                         const std::string        a_realm,
                         const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::computeGradient(EBAMRCellData, EBAMRCellData, string, phase::which_phase)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::computeGradient(EBAMRCellData, EBAMRCellData, string, phase::which_phase)" << endl;
  }

  CH_assert(a_gradient[0]->nComp() == SpaceDim);
  CH_assert(a_phi[0]->nComp() == 1);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

    const RefCountedPtr<EBGradient>& gradientOp = m_realms[a_realm]->getGradientOp(a_phase)[lvl];

    const bool hasFine = lvl < m_finestLevel;

    if (hasFine) {
      gradientOp->computeAMRGradient(*a_gradient[lvl], *a_phi[lvl], *a_phi[lvl + 1]);
    }
    else {
      gradientOp->computeLevelGradient(*a_gradient[lvl], *a_phi[lvl]);
    }
  }
}

void
AmrMesh::computeGradient(EBAMRFluxData&           a_gradient,
                         const EBAMRCellData&     a_phi,
                         const std::string        a_realm,
                         const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::computeGradient(EBAMRFluxData, EBAMRCellData, string, phase::which_phase)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::computeGradient(EBAMRFluxData, EBAMRCellData, string, phase::which_phase)" << endl;
  }

  // TLDR: This routine first computes the cell-centered gradient and it averages that to faces. For regular cells this will yield
  //       an 8-point stencil in 2D. We want to reduce that to a six-point stencil in 2D so that the normal derivative does not reach
  //       over a differencing distance of 2*dx (which the tangential derivatives do). So, after averaging the cell-centered gradient to
  //       faces we replace the component normal to the face with its tighter stencil variant.
  //
  //       Finally, we set the coarse-face gradients to be the arithmetic average of the fine faces. We do this to avoid having stencils
  //       that reach under the embedded boundary (where bogus data could be found).

  CH_assert(a_gradient[0]->nComp() == SpaceDim);
  CH_assert(a_phi[0]->nComp() == 1);

  // Make some scratch data.
  EBAMRCellData scratch;
  this->allocate(scratch, a_realm, a_phase, SpaceDim);
  this->computeGradient(scratch, a_phi, a_realm, a_phase);

  this->conservativeAverage(scratch, a_realm, a_phase);
  this->interpGhost(scratch, a_realm, a_phase);

  // Average the cells to face and replace the normal derivative with a tighter stencil.
  DataOps::averageCellToFace(a_gradient, scratch, this->getDomains());
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const RefCountedPtr<EBGradient>& gradientOp = m_realms[a_realm]->getGradientOp(a_phase)[lvl];

    gradientOp->computeNormalDerivative(*a_gradient[lvl], *a_phi[lvl]);
  }

  // Coarsen the faces.
  this->arithmeticAverage(a_gradient, a_realm, a_phase);
}

void
AmrMesh::computeGradient(MFAMRCellData& a_gradient, const MFAMRCellData& a_phi, const std::string a_realm) const
{
  CH_TIME("AmrMesh::computeGradient(MFAMRCellData, MFAMRCellData, string)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::computeGradient(MFAMRCellData, MFAMRCellData, string)" << endl;
  }

  for (int iphase = 0; iphase < m_multifluidIndexSpace->numPhases(); iphase++) {
    EBAMRCellData aliasGrad(1 + m_finestLevel);
    EBAMRCellData aliasPhi(1 + m_finestLevel);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      aliasGrad[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());
      aliasPhi[lvl]  = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());

      MultifluidAlias::aliasMF(*aliasGrad[lvl], iphase, *a_gradient[lvl]);
      MultifluidAlias::aliasMF(*aliasPhi[lvl], iphase, *a_phi[lvl]);

      CH_assert(aliasGrad[lvl]->nComp() == SpaceDim);
      CH_assert(aliasPhi[lvl]->nComp() == 1);
    }

    if (iphase == 0) {
      this->computeGradient(aliasGrad, aliasPhi, a_realm, phase::gas);
    }
    else if (iphase == 1) {
      this->computeGradient(aliasGrad, aliasPhi, a_realm, phase::solid);
    }
  }
}

void
AmrMesh::computeGradient(MFAMRFluxData& a_gradient, const MFAMRCellData& a_phi, const std::string a_realm) const
{
  CH_TIME("AmrMesh::computeGradient(MFAMRFluxData, MFAMRCellData, string)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::computeGradient(MFAMRFluxData, MFAMRCellData, string)" << endl;
  }

  for (int iphase = 0; iphase < m_multifluidIndexSpace->numPhases(); iphase++) {
    EBAMRFluxData aliasGrad(1 + m_finestLevel);
    EBAMRCellData aliasPhi(1 + m_finestLevel);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      aliasGrad[lvl] = RefCountedPtr<LevelData<EBFluxFAB>>(new LevelData<EBFluxFAB>());
      aliasPhi[lvl]  = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());

      MultifluidAlias::aliasMF(*aliasGrad[lvl], iphase, *a_gradient[lvl]);
      MultifluidAlias::aliasMF(*aliasPhi[lvl], iphase, *a_phi[lvl]);

      CH_assert(aliasGrad[lvl]->nComp() == SpaceDim);
      CH_assert(aliasPhi[lvl]->nComp() == 1);
    }

    if (iphase == 0) {
      this->computeGradient(aliasGrad, aliasPhi, a_realm, phase::gas);
    }
    else if (iphase == 1) {
      this->computeGradient(aliasGrad, aliasPhi, a_realm, phase::solid);
    }
  }
}

void
AmrMesh::average(MFAMRCellData& a_data, const std::string a_realm, const Average& a_average) const
{
  CH_TIME("AmrMesh::average(MFAMRCellData, string, Average)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::average(MFAMRCellData, string,Average)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str = "AmrMesh::average(MFAMRCellData) - could not find realm '" + a_realm + "'";

    MayDay::Abort(str.c_str());
  }

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_realms[a_realm]->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_realms[a_realm]->getEBIndexSpace(phase::solid);

  EBAMRCellData aliasGas(1 + m_finestLevel);
  EBAMRCellData aliasSol(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    aliasGas[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());
    aliasSol[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());

    if (!ebisGas.isNull()) {
      MultifluidAlias::aliasMF(*aliasGas[lvl], phase::gas, *a_data[lvl]);
    }
    if (!ebisSol.isNull()) {
      MultifluidAlias::aliasMF(*aliasSol[lvl], phase::solid, *a_data[lvl]);
    }
  }

  if (!ebisGas.isNull()) {
    this->average(aliasGas, a_realm, phase::gas, a_average);
  }
  if (!ebisSol.isNull()) {
    this->average(aliasSol, a_realm, phase::solid, a_average);
  }
}

void
AmrMesh::arithmeticAverage(MFAMRCellData& a_data, const std::string a_realm) const
{
  CH_TIME("AmrMesh::arithmeticAverage(MFAMRCellData, string)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::arithmeticAverage(MFAMRCellData, string)" << endl;
  }

  this->average(a_data, a_realm, Average::Arithmetic);
}

void
AmrMesh::harmonicAverage(MFAMRCellData& a_data, const std::string a_realm) const
{
  CH_TIME("AmrMesh::harmonic(MFAMRCellData, string)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::harmonic(MFAMRCellData, string)" << endl;
  }

  this->average(a_data, a_realm, Average::Harmonic);
}

void
AmrMesh::conservativeAverage(MFAMRCellData& a_data, const std::string a_realm) const
{
  CH_TIME("AmrMesh::conservativeAverage(MFAMRCellData, string)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::conservativeAverage(MFAMRCellData, string)" << endl;
  }

  this->average(a_data, a_realm, Average::Conservative);
}

void
AmrMesh::average(MFAMRFluxData& a_data, const std::string a_realm, const Average& a_average) const
{
  CH_TIME("AmrMesh::average(MFAMRFluxData, string)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::average(MFAMRFluxData, string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::average(MFAMRFluxData, string) - could not find realm '" + a_realm + "'";

    MayDay::Abort(str.c_str());
  }

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_realms[a_realm]->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_realms[a_realm]->getEBIndexSpace(phase::solid);

  // Alias the data to regular EBFluxFABs
  EBAMRFluxData aliasGas(1 + m_finestLevel);
  EBAMRFluxData aliasSol(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    aliasGas[lvl] = RefCountedPtr<LevelData<EBFluxFAB>>(new LevelData<EBFluxFAB>());
    aliasSol[lvl] = RefCountedPtr<LevelData<EBFluxFAB>>(new LevelData<EBFluxFAB>());

    if (!ebisGas.isNull()) {
      MultifluidAlias::aliasMF(*aliasGas[lvl], phase::gas, *a_data[lvl]);
    }
    if (!ebisSol.isNull()) {
      MultifluidAlias::aliasMF(*aliasSol[lvl], phase::solid, *a_data[lvl]);
    }
  }

  if (!ebisGas.isNull()) {
    this->average(aliasGas, a_realm, phase::gas, a_average);
  }
  if (!ebisSol.isNull()) {
    this->average(aliasSol, a_realm, phase::solid, a_average);
  }
}
void
AmrMesh::arithmeticAverage(MFAMRFluxData& a_data, const std::string a_realm) const
{
  CH_TIME("AmrMesh::arithmeticAverage(MFAMRFluxData, string)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::arithmeticAverage(MFAMRFluxData, string)" << endl;
  }

  this->average(a_data, a_realm, Average::Arithmetic);
}

void
AmrMesh::harmonicAverage(MFAMRFluxData& a_data, const std::string a_realm) const
{
  CH_TIME("AmrMesh::harmonicAverage(MFAMRFluxData, string)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::harmonicAverage(MFAMRFluxData, string)" << endl;
  }

  this->average(a_data, a_realm, Average::Harmonic);
}

void
AmrMesh::conservativeAverage(MFAMRFluxData& a_data, const std::string a_realm) const
{
  CH_TIME("AmrMesh::conservativeAverage(MFAMRFluxData, string)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::conservativeAverage(MFAMRFluxData, string)" << endl;
  }

  this->average(a_data, a_realm, Average::Conservative);
}

void
AmrMesh::average(EBAMRCellData&           a_data,
                 const std::string        a_realm,
                 const phase::which_phase a_phase,
                 const Average&           a_average) const
{
  CH_TIME("AmrMesh::average(EBAMRCellData, string, phase::which_phase, Average)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::average(EBAMRCellData, string, phase::which_phase, Average)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::average(EBAMRCellData) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finestLevel; lvl > 0; lvl--) {
    const int      nComps = a_data[lvl]->nComp();
    const Interval interv(0, nComps - 1);

    EBCoarAve& aveOp = *m_realms[a_realm]->getCoarseAverage(a_phase)[lvl];

    aveOp.averageData(*a_data[lvl - 1], *a_data[lvl], interv, a_average);
  }
}

void
AmrMesh::average(EBAMRCellData&           a_data,
                 const std::string        a_realm,
                 const phase::which_phase a_phase,
                 const Average&           a_average,
                 const Interval&          a_variables) const
{
  CH_TIME("AmrMesh::average(EBAMRCellData, string, phase::which_phase, Average, Interval)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::average(EBAMRCellData, string, phase::which_phase, Average, Interval)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::average(EBAMRCellData) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finestLevel; lvl > 0; lvl--) {
    CH_assert(a_variables.end() < a_data[lvl]->nComp());

    EBCoarAve& aveOp = *m_realms[a_realm]->getCoarseAverage(a_phase)[lvl];

    aveOp.averageData(*a_data[lvl - 1], *a_data[lvl], a_variables, a_average);
  }
}

void
AmrMesh::arithmeticAverage(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::average(EBAMRCellData, string, phase::which_phase)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::average(EBAMRCellData, string, phase::which_phase)" << endl;
  }

  this->average(a_data, a_realm, a_phase, Average::Arithmetic);
}

void
AmrMesh::arithmeticAverage(EBAMRCellData&           a_data,
                           const std::string        a_realm,
                           const phase::which_phase a_phase,
                           const Interval&          a_variables) const
{
  CH_TIME("AmrMesh::average(EBAMRCellData, string, phase::which_phase, Interval)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::average(EBAMRCellData, string, phase::which_phase, Interval)" << endl;
  }

  this->average(a_data, a_realm, a_phase, Average::Arithmetic, a_variables);
}

void
AmrMesh::harmonicAverage(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::harmonicAverage(EBAMRCellData, string, phase::which_phase)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::harmonicAverage(EBAMRCellData, string, phase::which_phase)" << endl;
  }

  this->average(a_data, a_realm, a_phase, Average::Harmonic);
}

void
AmrMesh::harmonicAverage(EBAMRCellData&           a_data,
                         const std::string        a_realm,
                         const phase::which_phase a_phase,
                         const Interval&          a_variables) const
{
  CH_TIME("AmrMesh::harmonicAverage(EBAMRCellData, string, phase::which_phase, Interval)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::harmonicAverage(EBAMRCellData, string, phase::which_phase, Interval)" << endl;
  }

  this->average(a_data, a_realm, a_phase, Average::Harmonic, a_variables);
}

void
AmrMesh::conservativeAverage(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::conservativeAverage(EBAMRCellData, string, phase::which_phase)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::conservativeAverage(EBAMRCellData, string, phase::which_phase)" << endl;
  }

  this->average(a_data, a_realm, a_phase, Average::Conservative);
}

void
AmrMesh::conservativeAverage(EBAMRCellData&           a_data,
                             const std::string        a_realm,
                             const phase::which_phase a_phase,
                             const Interval&          a_variables) const
{
  CH_TIME("AmrMesh::conservativeAverage(EBAMRCellData, string, phase::which_phase, Interval)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::conservativeAverage(EBAMRCellData, string, phase::which_phase, Interval)" << endl;
  }

  this->average(a_data, a_realm, a_phase, Average::Conservative, a_variables);
}

void
AmrMesh::average(EBAMRFluxData&           a_data,
                 const std::string        a_realm,
                 const phase::which_phase a_phase,
                 const Average&           a_average) const
{
  CH_TIME("AmrMesh::arithmeticAverage(EBAMRFluxData, string, phase::which_phase");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::arithmeticAverage(EBAMRFluxData, string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::average(EBAMRFluxData) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finestLevel; lvl > 0; lvl--) {
    const int      nComps = a_data[lvl]->nComp();
    const Interval interv(0, nComps - 1);

    EBCoarAve& aveOp = *m_realms[a_realm]->getCoarseAverage(a_phase)[lvl];

    aveOp.averageData(*a_data[lvl - 1], *a_data[lvl], interv, a_average);
  }
}

void
AmrMesh::arithmeticAverage(EBAMRFluxData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::arithmeticAverage(EBAMRFluxData, string, phase::which_phase");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::arithmeticAverage(EBAMRFluxData, string, phase::which_phase)" << endl;
  }

  this->average(a_data, a_realm, a_phase, Average::Arithmetic);
}

void
AmrMesh::harmonicAverage(EBAMRFluxData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::harmonicAverage(EBAMRFluxData, string, phase::which_phase");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::harmonicAverage(EBAMRFluxData, string, phase::which_phase)" << endl;
  }

  this->average(a_data, a_realm, a_phase, Average::Harmonic);
}

void
AmrMesh::conservativeAverage(EBAMRFluxData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::conservativeAverage(EBAMRFluxData, string, phase::which_phase");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::conservativeAverage(EBAMRFluxData, string, phase::which_phase)" << endl;
  }

  this->average(a_data, a_realm, a_phase, Average::Conservative);
}

void
AmrMesh::average(EBAMRIVData&             a_data,
                 const std::string        a_realm,
                 const phase::which_phase a_phase,
                 const Average&           a_average) const
{
  CH_TIME("AmrMesh::arithmeticAverage(EBAMRIVData, string, phase::which_phase");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::arithmeticAverage(EBAMRIVData, string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::average(EBAMRIVData) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finestLevel; lvl > 0; lvl--) {
    const int      nComps = a_data[lvl]->nComp();
    const Interval interv(0, nComps - 1);

    EBCoarAve& aveOp = *m_realms[a_realm]->getCoarseAverage(a_phase)[lvl];

    aveOp.averageData(*a_data[lvl - 1], *a_data[lvl], interv, a_average);
  }
}

void
AmrMesh::arithmeticAverage(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::arithmeticAverage(EBAMRIVData, string, phase::which_phase");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::arithmeticAverage(EBAMRIVData, string, phase::which_phase)" << endl;
  }

  this->average(a_data, a_realm, a_phase, Average::Arithmetic);
}

void
AmrMesh::harmonicAverage(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::harmonicAverage(EBAMRIVData, string, phase::which_phase");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::harmonicAverage(EBAMRIVData, string, phase::which_phase)" << endl;
  }

  this->average(a_data, a_realm, a_phase, Average::Harmonic);
}

void
AmrMesh::conservativeAverage(EBAMRIVData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::conservativeAverage(EBAMRFluxData, string, phase::which_phase");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::conservativeAverage(EBAMRFluxData, string, phase::which_phase)" << endl;
  }

  this->average(a_data, a_realm, a_phase, Average::Conservative);
}

void
AmrMesh::interpGhost(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::interpGhost(EBAMRCellData, string, phase::which_phase)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::interpGhost(EBAMRCellData, string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str =
      "AmrMesh::interpGhost(EBAMRCellData, string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  this->interpGhostPwl(a_data, a_realm, a_phase);
}

void
AmrMesh::interpGhost(LevelData<EBCellFAB>&       a_fineData,
                     const LevelData<EBCellFAB>& a_coarData,
                     const int                   a_fineLevel,
                     const std::string           a_realm,
                     const phase::which_phase    a_phase) const
{
  CH_TIME("AmrMesh::interpGhost(LD<EBCellFAB>, LD<EBCellFAB>, int, string, phase::which_phase)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::interpGhost(LD<EBCellFAB>, LD<EBCellFAB>, int, string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str =
      "AmrMesh::interpGhost(LD<EBCellFAB>, LD<EBCellFAB>, int, string, phase::which_phase) - could not find realm '" +
      a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  if (a_fineLevel > 0) {

    const int      nComps = a_fineData.nComp();
    const Interval interv = Interval(0, nComps - 1);

    AggEBPWLFillPatch& fillpatch = *m_realms[a_realm]->getFillPatch(a_phase)[a_fineLevel];

    fillpatch.interpolate(a_fineData, a_coarData, a_coarData, 0.0, 0.0, 0.0, interv);
  }

  a_fineData.exchange();
}

void
AmrMesh::interpGhost(MFAMRCellData& a_data, const std::string a_realm) const
{
  CH_TIME("AmrMesh::interpGhost(MFAMRCellData, string)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::interpGhost(MFAMRCellData, string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::interpGhost(MFAMRCellData, string) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  // Do aliasing
  EBAMRCellData aliasGas(1 + m_finestLevel);
  EBAMRCellData aliasSol(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_realms[a_realm]->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_realms[a_realm]->getEBIndexSpace(phase::solid);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    aliasGas[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());
    aliasSol[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());

    if (!ebisGas.isNull()) {
      MultifluidAlias::aliasMF(*aliasGas[lvl], phase::gas, *a_data[lvl]);
    }
    if (!ebisSol.isNull()) {
      MultifluidAlias::aliasMF(*aliasSol[lvl], phase::solid, *a_data[lvl]);
    }
  }

  if (!ebisGas.isNull()) {
    this->interpGhost(aliasGas, a_realm, phase::gas);
  }
  if (!ebisSol.isNull()) {
    this->interpGhost(aliasSol, a_realm, phase::solid);
  }
}

void
AmrMesh::interpGhostPwl(MFAMRCellData& a_data, const std::string a_realm) const
{
  CH_TIME("AmrMesh::interpGhostPwl(MFAMRCellData, string)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::interpGhostPwl(MFAMRCellData, string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::interpGhostPwl(MFAMRCellData, string) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  // Do aliasing
  EBAMRCellData aliasGas(1 + m_finestLevel);
  EBAMRCellData aliasSol(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_realms[a_realm]->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_realms[a_realm]->getEBIndexSpace(phase::solid);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    aliasGas[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());
    aliasSol[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());

    if (!ebisGas.isNull())
      MultifluidAlias::aliasMF(*aliasGas[lvl], phase::gas, *a_data[lvl]);
    if (!ebisSol.isNull())
      MultifluidAlias::aliasMF(*aliasSol[lvl], phase::solid, *a_data[lvl]);
  }

  if (!ebisGas.isNull()) {
    this->interpGhostPwl(aliasGas, a_realm, phase::gas);
  }
  if (!ebisSol.isNull()) {
    this->interpGhostPwl(aliasSol, a_realm, phase::solid);
  }
}

void
AmrMesh::interpGhostPwl(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::interpGhostPwl(EBAMRCellData, string, phase::which_phase)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::interpGhostPwl(EBAMRCellData, string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::interpGhostPwl(ebamrcell, Realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = m_finestLevel; lvl > 0; lvl--) {
    const int      nComps = a_data[lvl]->nComp();
    const Interval interv(0, nComps - 1);

    AggEBPWLFillPatch& fillpatch = *m_realms[a_realm]->getFillPatch(a_phase)[lvl];

    fillpatch.interpolate(*a_data[lvl], *a_data[lvl - 1], *a_data[lvl - 1], 0.0, 0.0, 0.0, interv);
  }

  a_data.exchange();
}

void
AmrMesh::interpGhostMG(MFAMRCellData& a_data, const std::string a_realm) const
{
  CH_TIME("AmrMesh::interpGhostMG(MFAMRCellData, string)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::interpGhostMG(MFAMRCellData, string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::interpGhostMG(MFAMRCellData, string) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  // Do aliasing
  EBAMRCellData aliasGas(1 + m_finestLevel);
  EBAMRCellData aliasSol(1 + m_finestLevel);

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_realms[a_realm]->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_realms[a_realm]->getEBIndexSpace(phase::solid);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    aliasGas[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());
    aliasSol[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());

    if (!ebisGas.isNull()) {
      MultifluidAlias::aliasMF(*aliasGas[lvl], phase::gas, *a_data[lvl]);
    }
    if (!ebisSol.isNull()) {
      MultifluidAlias::aliasMF(*aliasSol[lvl], phase::solid, *a_data[lvl]);
    }
  }

  if (!ebisGas.isNull()) {
    this->interpGhostMG(aliasGas, a_realm, phase::gas);
  }
  if (!ebisSol.isNull()) {
    this->interpGhostMG(aliasSol, a_realm, phase::solid);
  }
}

void
AmrMesh::interpGhostMG(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::interpGhostPwl(EBAMRCellData, string, phase::which_phase)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::interpGhostPwl(EBAMRCellData, string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str = "AmrMesh::interpGhostPwl(ebamrcell, Realm, phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  for (int lvl = 1; lvl <= m_finestLevel; lvl++) {
    const int nComps = a_data[lvl]->nComp();

    RefCountedPtr<EBMultigridInterpolator>& interpolator = m_realms[a_realm]->getMultigridInterpolator(a_phase)[lvl];

    interpolator->coarseFineInterp(*a_data[lvl], *a_data[lvl - 1], Interval(0, nComps - 1));
  }

  a_data.exchange();
}

void
AmrMesh::interpToNewGrids(MFAMRCellData&       a_newData,
                          const MFAMRCellData& a_oldData,
                          const int            a_lmin,
                          const int            a_oldFinestLevel,
                          const int            a_newFinestLevel,
                          const bool           a_useSlopes)
{
  CH_TIME("AmrMesh::interpToNewGrids(MFAMRCellData x2, int x3, bool)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::interpToNewGrids(MFAMRCellData x2, int x3, bool)" << endl;
  }

  for (int i = 0; i < phase::numPhases; i++) {
    phase::which_phase curPhase;

    if (i == 0) {
      curPhase = phase::gas;
    }
    else {
      curPhase = phase::solid;
    }

    const RefCountedPtr<EBIndexSpace>& ebis = m_multifluidIndexSpace->getEBIndexSpace(curPhase);

    if (!(ebis.isNull())) {
      EBAMRCellData       newData = this->alias(curPhase, a_newData);
      const EBAMRCellData oldData = this->alias(curPhase, a_oldData);

      this->interpToNewGrids(newData, oldData, curPhase, a_lmin, a_oldFinestLevel, a_newFinestLevel, a_useSlopes);
    }
  }
}

void
AmrMesh::interpToNewGrids(EBAMRCellData&           a_newData,
                          const EBAMRCellData&     a_oldData,
                          const phase::which_phase a_phase,
                          const int                a_lmin,
                          const int                a_oldFinestLevel,
                          const int                a_newFinestLevel,
                          const bool               a_useSlopes)
{
  CH_TIME("AmrMesh::interpToNewGrids(EBAMRCellData x2, phase, int x3, bool)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::interpToNewGrids(EBAMRCellData x2, phase, int x3, bool)" << endl;
  }

  CH_assert(a_newData.getRealm() == a_oldData.getRealm());
  CH_assert(a_newData[0]->nComp() == a_oldData[0]->nComp());

  const int nComp = a_newData[0]->nComp();

  // These levels have not changed but ownship MIGHT have changed.
  for (int lvl = 0; lvl <= std::max(0, a_lmin - 1); lvl++) {
    a_oldData[lvl]->copyTo(*a_newData[lvl]);
  }

  // These levels have changed.
  for (int lvl = std::max(1, a_lmin); lvl <= a_newFinestLevel; lvl++) {
    RefCountedPtr<EBFineInterp>& interpolator = this->getFineInterp(a_newData.getRealm(), a_phase)[lvl];

    // Interpolate the data.
    if (a_useSlopes) {
      interpolator->regridMinMod(*a_newData[lvl], *a_newData[lvl - 1], Interval(0, nComp - 1));
    }
    else {
      interpolator->regridNoSlopes(*a_newData[lvl], *a_newData[lvl - 1], Interval(0, nComp - 1));
    }

    // There could be parts of the new grid that overlapped with the old grid (on level lvl) -- we don't want
    // to pollute the solution with interpolation there since we already have valid data.
    if (lvl <= std::min(a_oldFinestLevel, a_newFinestLevel)) {
      a_oldData[lvl]->copyTo(*a_newData[lvl]);
    }
  }
}

void
AmrMesh::interpToCentroids(EBAMRCellData& a_data, const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::interpToCentroids(EBAMRCellData, string, phase::which_phase)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::interpToCentroids(EBAMRCellData, string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    std::string str =
      "AmrMesh::interpToCentroids(EBAMRCellData, string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  const IrregAmrStencil<CentroidInterpolationStencil>& stencil =
    m_realms[a_realm]->getCentroidInterpolationStencils(a_phase);
  stencil.apply(a_data);
}

void
AmrMesh::setCoarsestGrid(const IntVect& a_nCells)
{
  CH_TIME("AmrMesh::setCoarsestGrid()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::setCoarsestGrid()" << endl;
  }

  m_numCells = a_nCells;
}

void
AmrMesh::parseVerbosity()
{
  CH_TIME("AmrMesh::parseVerbosity");

  ParmParse pp("AmrMesh");
  pp.get("verbosity", m_verbosity);

  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseVerbosity()" << endl;
  }
}

void
AmrMesh::parseCoarsestLevelNumCells()
{
  CH_TIME("AmrMesh::parseCoarsestLevelNumCells()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseCoarsestLevelNumCells()" << endl;
  }

  ParmParse   pp("AmrMesh");
  Vector<int> cells;
  cells.resize(pp.countval("coarsest_domain"));
  CH_assert(cells.size() >= SpaceDim);
  pp.getarr("coarsest_domain", cells, 0, SpaceDim);

  m_numCells = IntVect(D_DECL(cells[0], cells[1], cells[2]));
}

void
AmrMesh::parseMaxAmrDepth()
{
  CH_TIME("AmrMesh::parseMaxAmrDepth()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseMaxAmrDepth()" << endl;
  }

  ParmParse pp("AmrMesh");
  int       depth;
  pp.get("max_amr_depth", depth);
  if (depth >= 0) {
    m_maxAmrDepth = depth;
  }
  else {
    m_maxAmrDepth = 0;
  }
}

void
AmrMesh::parseMaxSimulationDepth()
{
  CH_TIME("AmrMesh::parseMaxSimulationDepth()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseMaxSimulationDepth()" << endl;
  }

  ParmParse pp("AmrMesh");
  int       depth;
  pp.get("max_sim_depth", depth);
  if (depth >= 0) {
    m_maxSimulationDepth = depth;
  }
  else {
    m_maxSimulationDepth = m_maxAmrDepth;
  }
}

void
AmrMesh::parseRefinementRatios()
{
  CH_TIME("AmrMesh::parseRefinementRatios()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseRefinementRatios()" << endl;
  }

  ParmParse   pp("AmrMesh");
  Vector<int> ratios;
  ratios.resize(pp.countval("ref_rat"));
  pp.getarr("ref_rat", ratios, 0, ratios.size());

  // Pad with whatever was last specified if user didn't supply enough refinement factors
  while (ratios.size() < m_maxAmrDepth) {
    //    ratios.push_back(2);
    ratios.push_back(ratios.back());
  }

  m_refinementRatios = ratios;
}

void
AmrMesh::parseBrFillRatio()
{
  CH_TIME("AmrMesh::parseBrFillRatio()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseBrFillRatio()" << endl;
  }

  ParmParse pp("AmrMesh");
  Real      fill_ratio = 1.0;
  pp.get("fill_ratio", fill_ratio);
  if (fill_ratio > 0.0 && fill_ratio <= 1.0) {
    m_fillRatioBR = fill_ratio;
  }
}

void
AmrMesh::setFinestLevel(const int a_finestLevel)
{
  CH_TIME("AmrMesh::setFinestLevel(int)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::setFinestLevel(int)" << endl;
  }

  m_finestLevel = a_finestLevel;
  m_finestLevel = Min(m_finestLevel, m_maxAmrDepth);        // Don't exceed m_maxAmrDepth
  m_finestLevel = Min(m_finestLevel, m_maxSimulationDepth); // Don't exceed maximum simulation depth
}

void
AmrMesh::setGrids(const Vector<Vector<Box>>&                             a_boxes,
                  const std::map<std::string, Vector<Vector<long int>>>& a_realmsAndLoads)
{
  CH_TIME("AmrMesh::setGrids(Vector<Vector<Box> >, std::map<string, Vector<Vector<long int> >)");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::setGrids(Vector<Vector<Box> >, std::map<string, Vector<Vector<long int> >)" << endl;
  }

  // TLDR: This routine is called by driver when restarting from an HDF5 file. The boxes are the same for all realms, but
  //       the restart feature can use checkpointed loads. These are read from file and sent to this routine so we can call
  //       LoadBalancing.

  const int lmin = 0;

  for (const auto& r : a_realmsAndLoads) {
    const std::string&              curRealm = r.first;
    const Vector<Vector<long int>>& curLoads = r.second;

    Vector<Vector<int>> processorIDs(1 + m_finestLevel);

    // Do load balancing.
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      LoadBalancing::makeBalance(processorIDs[lvl], curLoads[lvl], a_boxes[lvl]);
    }

    this->regridRealm(curRealm, processorIDs, a_boxes, lmin);
  }

  // Set the proxy grids, too.
  m_grids    = m_realms[Realm::Primal]->getGrids();
  m_hasGrids = true;
}

void
AmrMesh::parseMaxBoxSize()
{
  CH_TIME("AmrMesh::parseMaxBoxSize()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseMaxBoxSize()" << endl;
  }

  ParmParse pp("AmrMesh");
  int       boxSize;
  pp.get("max_box_size", boxSize);
  if (boxSize >= 4 && boxSize % 2 == 0) {
    m_maxBoxSize = boxSize;
  }
  else {
    MayDay::Error("AmrMesh::parseMaxBoxSize - must have box_size >= 4 and divisible by 2");
  }
}

void
AmrMesh::parseMaxEbisBoxSize()
{
  CH_TIME("AmrMesh::parseMaxEbisBoxSize()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseMaxEbisBoxSize()" << endl;
  }

  ParmParse pp("AmrMesh");
  int       boxSize;
  pp.get("max_ebis_box", boxSize);
  if (boxSize >= 4 && boxSize % 2 == 0) {
    m_maxEbisBoxSize = boxSize;
  }
  else {
    MayDay::Abort("AmrMesh::parseMaxEbisBoxSize - must have box_size >= 8 and divisible by 2");
  }
}

void
AmrMesh::parseBrBufferSize()
{
  CH_TIME("AmrMesh::parseBrBufferSize()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseBrBufferSize()" << endl;
  }

  ParmParse pp("AmrMesh");
  int       buffer = 2;
  pp.get("buffer_size", buffer);
  if (buffer > 0) {
    m_bufferSizeBR = buffer;
  }
  else {
    MayDay::Error("AmrMesh::parseBrBufferSize() - must have buffer > 0");
  }
}

void
AmrMesh::parseGridGeneration()
{
  CH_TIME("AmrMesh::parseGridGeneration()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseGridGeneration()" << endl;
  }

  ParmParse   pp("AmrMesh");
  std::string str;
  pp.get("grid_algorithm", str);
  if (str == "br") {
    m_gridGenerationMethod = GridGenerationMethod::BergerRigoutsous;
  }
  else if (str == "tiled") {
    m_gridGenerationMethod = GridGenerationMethod::Tiled;
  }
  else {
    MayDay::Abort("AmrMesh::parseGridGeneration - unknown grid generation method requested");
  }

  pp.get("box_sorting", str);
  if (str == "none") {
    m_boxSort = BoxSorting::None;
  }
  else if (str == "std") {
    m_boxSort = BoxSorting::Std;
  }
  else if (str == "shuffle") {
    m_boxSort = BoxSorting::Shuffle;
  }
  else if (str == "morton") {
    m_boxSort = BoxSorting::Morton;
  }
  else {
    MayDay::Abort("AmrMesh::parseGridGeneration - unknown box sorting method requested");
  }
}

void
AmrMesh::parseBlockingFactor()
{
  CH_TIME("AmrMesh::parseBlockingFactor()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseBlockingFactor()" << endl;
  }

  ParmParse pp("AmrMesh");
  int       blocking;
  pp.get("blocking_factor", blocking);
  if (blocking >= 4 && blocking % 2 == 0) {
    m_blockingFactor = blocking;
  }
}

void
AmrMesh::parseEbGhostCells()
{
  CH_TIME("AmrMesh::parseEbGhostCells()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseEbGhostCells()" << endl;
  }

  ParmParse pp("AmrMesh");

  pp.get("eb_ghost", m_numEbGhostsCells);

  if (m_numEbGhostsCells < 0)
    MayDay::Error("AmrMesh::parseEbGhostCells -- you have specified a negative number of ghost cells");
}

void
AmrMesh::parseNumGhostCells()
{
  CH_TIME("AmrMesh::parseNumGhostCells()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseNumGhostCells()" << endl;
  }

  ParmParse pp("AmrMesh");

  pp.get("num_ghost", m_numGhostCells);
  pp.get("lsf_ghost", m_numLsfGhostCells);

  if (m_numGhostCells < 0)
    MayDay::Error("AmrMesh::parseNumGhostCells -- you have specified a negative number of ghost cells for mesh data");
  if (m_numLsfGhostCells < 0)
    MayDay::Error(
      "AmrMesh::parseNumGhostCells -- you have specified a negative number of ghost cells for level-set mesh data");
}

void
AmrMesh::parseMultigridInterpolator()
{
  CH_TIME("AmrMesh::parseMultigridInterpolator()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseMultigridInterpolator()" << endl;
  }

  ParmParse pp("AmrMesh");

  pp.get("mg_interp_order", m_multigridInterpOrder);
  pp.get("mg_interp_radius", m_multigridInterpRadius);
  pp.get("mg_interp_weight", m_multigridInterpWeight);

  if (m_multigridInterpOrder < 0)
    MayDay::Error("AmrMesh::parseMultigridInterpolator -- you have specified negative order!");
  if (m_multigridInterpRadius <= 0)
    MayDay::Error("AmrMesh::parseMultigridInterpolator -- you have specified a non-positive radius!");
  if (m_multigridInterpWeight < 0)
    MayDay::Error("AmrMesh::parseMultigridInterpolator -- you have specified negative weighting");
}

void
AmrMesh::parseRedistributionRadius()
{
  CH_TIME("AmrMesh::parseRedistributionRadius()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseRedistributionRadius()" << endl;
  }

  ParmParse pp("AmrMesh");

  pp.get("redist_radius", m_redistributionRadius);

  if (m_redistributionRadius <= 0)
    MayDay::Error("AmrMesh::parseRedistributionRadius -- you have specified non-positive redistribution radius");
}

void
AmrMesh::parseCentroidStencils()
{
  CH_TIME("AmrMesh::parseCentroidStencils()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseCentroidStencils()" << endl;
  }

  ParmParse pp("AmrMesh");

  std::string str;

  pp.get("centroid_sten", str);

  // Maybe, in the future, we can change these but the user should not care about these (yet)
  m_centroidStencilRadius = 1;
  m_centroidStencilOrder  = 1;

  if (str == "linear") {
    m_centroidStencilType = IrregStencil::StencilType::Linear;
  }
  else if (str == "taylor") {
    m_centroidStencilType = IrregStencil::StencilType::TaylorExtrapolation;
  }
  else if (str == "lsq") {
    m_centroidStencilType = IrregStencil::StencilType::LeastSquares;
  }
  else if (str == "pwl") {
    m_centroidStencilType = IrregStencil::StencilType::PiecewiseLinear;
  }
  else {
    MayDay::Abort("AmrMesh::parseCentroidStencils - unknown stencil requested");
  }
}

void
AmrMesh::parseEbCentroidStencils()
{
  CH_TIME("AmrMesh::parseEbCentroidStencils()");
  if (m_verbosity > 3) {
    pout() << "AmrMesh::parseEbCentroidStencils()" << endl;
  }

  ParmParse pp("AmrMesh");

  std::string str;

  pp.get("eb_sten", str);

  // Maybe, in the future, we can change these but the user should not care about these (yet)
  m_ebCentroidStencilRadius = 1;
  m_ebCentroidStencilOrder  = 1;

  if (str == "linear") {
    m_ebCentroidStencilType = IrregStencil::StencilType::Linear;
  }
  else if (str == "taylor") {
    m_ebCentroidStencilType = IrregStencil::StencilType::TaylorExtrapolation;
  }
  else if (str == "lsq") {
    m_ebCentroidStencilType = IrregStencil::StencilType::LeastSquares;
  }
  else if (str == "pwl") {
    m_ebCentroidStencilType = IrregStencil::StencilType::PiecewiseLinear;
  }
  else {
    MayDay::Abort("AmrMesh::parseEbCentroidStencils - unknown stencil requested");
  }
}

void
AmrMesh::sanityCheck() const
{
  CH_TIME("AmrMesh::sanityCheck()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::sanityCheck()" << endl;
  }

  CH_assert(m_maxAmrDepth >= 0);
  for (int lvl = 0; lvl < m_refinementRatios.size(); lvl++) {
    CH_assert(m_refinementRatios[lvl] == 2 || m_refinementRatios[lvl] == 4);
    CH_assert(m_blockingFactor >= 4 && m_blockingFactor % m_refinementRatios[lvl] == 0);
  }

  CH_assert(m_maxBoxSize >= 8 && m_maxBoxSize % m_blockingFactor == 0);
  CH_assert(m_fillRatioBR > 0. && m_fillRatioBR <= 1.0);
  CH_assert(m_bufferSizeBR > 0);

  if (m_maxAmrDepth > 0) {
    for (int lvl = 0; lvl < m_maxAmrDepth; lvl++) {
      if (m_refinementRatios[lvl] > 2 && m_blockingFactor < 8) {
        MayDay::Abort("AmrMesh::sanityCheck -- can't use blocking factor < 8 with factor 4 refinement!");
      }
    }
  }
}

RealVect
AmrMesh::getProbLo() const
{
  CH_TIME("AmrMesh::getProbLo()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getProbLo()" << endl;
  }

  return m_probLo;
}

RealVect
AmrMesh::getProbHi() const
{
  CH_TIME("AmrMesh::getProbHi()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getProbHi()" << endl;
  }

  return m_probHi;
}

int
AmrMesh::getFinestLevel() const
{
  CH_TIME("AmrMesh::getFinestLevel()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getFinestLevel()" << endl;
  }

  return m_finestLevel;
}

int
AmrMesh::getMaxAmrDepth() const
{
  CH_TIME("AmrMesh::getMaxAmrDepth()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getMaxAmrDepth()" << endl;
  }

  return m_maxAmrDepth;
}

int
AmrMesh::getMaxSimulationDepth() const
{
  CH_TIME("AmrMesh::getMaxSimulationDepth()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getMaxSimulationDepth()" << endl;
  }

  return m_maxSimulationDepth;
}

int
AmrMesh::getNumberOfGhostCells() const
{
  CH_TIME("AmrMesh::getNumberOfGhostCells()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getNumberOfGhostCells()" << endl;
  }

  return m_numGhostCells;
}

int
AmrMesh::getNumberOfEbGhostCells() const
{
  CH_TIME("AmrMesh::getNumberOfEbGhostCells()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getNumberOfEbGhostCells()" << endl;
  }

  return m_numEbGhostsCells;
}

int
AmrMesh::getRedistributionRadius() const
{
  CH_TIME("AmrMesh::getRedistributionRadius()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getRedistributionRadius()" << endl;
  }

  return m_redistributionRadius;
}

int
AmrMesh::getBlockingFactor() const
{
  CH_TIME("AmrMesh::getBlockingFactor()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getBlockingFactor()" << endl;
  }

  return m_blockingFactor;
}

int
AmrMesh::getMaxBoxSize() const
{
  CH_TIME("AmrMesh::getMaxBoxSize()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getMaxBoxSize()" << endl;
  }

  return m_maxBoxSize;
}

int
AmrMesh::getBrBuffer() const
{
  CH_TIME("AmrMesh::getBrBuffer()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getBrBuffer()" << endl;
  }

  return m_bufferSizeBR;
}

int
AmrMesh::getMaxEbisBoxSize() const
{
  CH_TIME("AmrMesh::getMaxEbisBoxSize()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getMaxEbisBoxSize()" << endl;
  }

  return m_maxEbisBoxSize;
}

int
AmrMesh::getRefinementRatio(const int a_level1, const int a_level2) const
{
  CH_TIME("AmrMesh::getRefinementRatio(int, int)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getRefinementRatio(int, int)" << endl;
  }

  int coarLevel = Min(a_level1, a_level2);
  int fineLevel = Max(a_level1, a_level2);

  int ref = 1;
  for (int lvl = coarLevel; lvl < fineLevel; lvl++) {
    ref = ref * m_refinementRatios[lvl];
  }

  return ref;
}

ProblemDomain
AmrMesh::getFinestDomain() const
{
  CH_TIME("AmrMesh::getFinestDomain()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getFinestDomain()" << endl;
  }

  return m_domains[m_maxAmrDepth];
}

Real
AmrMesh::getFinestDx() const
{
  CH_TIME("AmrMesh::getFinestDx()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getFinestDx()" << endl;
  }

  return m_dx[m_maxAmrDepth];
}

const Vector<Real>&
AmrMesh::getDx() const
{
  CH_TIME("AmrMesh::getDx()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getDx()" << endl;
  }

  return m_dx;
}

const Vector<int>&
AmrMesh::getRefinementRatios() const
{
  CH_TIME("AmrMesh::getRefinementRatios()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getRefinementRatios()" << endl;
  }

  return m_refinementRatios;
}

const RefCountedPtr<BaseIF>&
AmrMesh::getBaseImplicitFunction(const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getBaseImplicitFunctions()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getBaseImplicitFunctions()" << endl;
  }

  return m_baseif.at(a_phase);
}

const Vector<ProblemDomain>&
AmrMesh::getDomains() const
{
  CH_TIME("AmrMesh::getDomains()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getDomains()" << endl;
  }

  return m_domains;
}

const Vector<DisjointBoxLayout>&
AmrMesh::getProxyGrids() const
{
  CH_TIME("AmrMesh::getProxyGrids()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getProxyGrids()" << endl;
  }

  return m_grids;
}

const Vector<DisjointBoxLayout>&
AmrMesh::getGrids(const std::string a_realm) const
{
  CH_TIME("AmrMesh::getGrids(string)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getGrids(string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str = "AmrMesh::getGrids - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getGrids();
}

const Vector<EBISLayout>&
AmrMesh::getEBISLayout(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getEBISLayout(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getEBISLayout(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getEBISLayout(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getEBISLayout(a_phase);
}

Vector<RefCountedPtr<LayoutData<VoFIterator>>>&
AmrMesh::getVofIterator(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getVofIterator(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getVofIterator(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getVofIterator(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getVofIterator(a_phase);
}

const Vector<RefCountedPtr<LayoutData<Vector<LayoutIndex>>>>&
AmrMesh::getNeighbors(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getNeighbors(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getNeighbors(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getNeighbors(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getNeighbors(a_phase);
}

const AMRMask&
AmrMesh::getMask(const std::string a_mask, const int a_buffer, const std::string a_realm) const
{
  CH_TIME("AmrMesh::getMask(string, int, string)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getMask(string, int, string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str = "AmrMesh::getMask(string, int, string) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getMask(a_mask, a_buffer);
}

const AMRMask&
AmrMesh::getValidCells(const std::string a_realm) const
{
  CH_TIME("AmrMesh::getValidCells(string)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getValidCells(string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str = "AmrMesh::getValidCells(string) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getValidCells();
}

const Vector<RefCountedPtr<EBLevelGrid>>&
AmrMesh::getEBLevelGrid(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getEBLevelGrid(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getEBLevelGrid(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getEBLevelGrid(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getEBLevelGrid(a_phase);
}

const Vector<RefCountedPtr<MFLevelGrid>>&
AmrMesh::getMFLevelGrid(const std::string a_realm) const
{
  CH_TIME("AmrMesh::getMFLevelGrid(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getMFLevelGrid(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getMFLevelGrid(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getMFLevelGrid();
}

const EBAMRFAB&
AmrMesh::getLevelset(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getLevelSet(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getLevelSet(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str = "AmrMesh::getLevelSet(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getLevelset(a_phase);
}

EBAMRParticleMesh&
AmrMesh::getParticleMesh(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getParticleMesh(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getParticleMesh(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str = "AmrMesh::getLevelSet(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getParticleMesh(a_phase);
}

Vector<RefCountedPtr<EBCoarAve>>&
AmrMesh::getCoarseAverage(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getCoarseAverage(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getCoarseAverage(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getCoarseAverage(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getCoarseAverage(a_phase);
}

Vector<RefCountedPtr<EBMultigridInterpolator>>&
AmrMesh::getMultigridInterpolator(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getMultigridInterpolator(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getMultigridInterpolator(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getMultigridInterpolator(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getMultigridInterpolator(a_phase);
}

Vector<RefCountedPtr<AggEBPWLFillPatch>>&
AmrMesh::getFillPatch(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getFillPatch(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getFillPatch(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getFillPatch(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getFillPatch(a_phase);
}

Vector<RefCountedPtr<EBFineInterp>>&
AmrMesh::getFineInterp(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getFineInterp(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getFineInterp(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getFineInterp(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getFineInterp(a_phase);
}

Vector<RefCountedPtr<EBFluxRegister>>&
AmrMesh::getFluxRegister(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getFluxRegister(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getFluxRegister(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getFluxRegister(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getFluxRegister(a_phase);
}

Vector<RefCountedPtr<EBLevelRedist>>&
AmrMesh::getLevelRedist(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getLevelRedist(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getLevelRedist(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getLevelRedist(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getLevelRedist(a_phase);
}

Vector<RefCountedPtr<EBCoarToFineRedist>>&
AmrMesh::getCoarToFineRedist(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getCoarToFineRedist(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getCoarToFineRedist(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getCoarToFineRedist(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getCoarToFineRedist(a_phase);
}

Vector<RefCountedPtr<EBCoarToCoarRedist>>&
AmrMesh::getCoarToCoarRedist(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getCoarToCoarRedist(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getCoarToCoarRedist(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getCoarToCoarRedist(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getCoarToCoarRedist(a_phase);
}

Vector<RefCountedPtr<EBFineToCoarRedist>>&
AmrMesh::getFineToCoarRedist(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getFineToCoarRedist(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getFineToCoarRedist(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getFineToCoarRedist(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getFineToCoarRedist(a_phase);
}

const IrregAmrStencil<CentroidInterpolationStencil>&
AmrMesh::getCentroidInterpolationStencils(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getCentroidInterpolationStencil(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getCentroidInterpolationStencil(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getCentroidInterpolationStencil(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getCentroidInterpolationStencils(a_phase);
}

const IrregAmrStencil<EbCentroidInterpolationStencil>&
AmrMesh::getEbCentroidInterpolationStencils(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getEbCentroidInterpolationStencil(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getEbCentroidInterpolationStencil(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getEbCentroidInterpolationStencil(string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getEbCentroidInterpolationStencilStencils(a_phase);
}

const IrregAmrStencil<NonConservativeDivergenceStencil>&
AmrMesh::getNonConservativeDivergenceStencils(const std::string a_realm, const phase::which_phase a_phase) const
{
  CH_TIME("AmrMesh::getNonConservativeDivergenceStencil(string, phase::which_phase)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getNonConservativeDivergenceStencil(string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::getNonConservativeDivergenceStencil(string, phase::which_phase) - could not find realm '" + a_realm +
      "'";
    MayDay::Abort(str.c_str());
  }

  return m_realms[a_realm]->getNonConservativeDivergenceStencils(a_phase);
}

bool
AmrMesh::queryRealm(const std::string a_realm) const
{
  CH_TIME("AmrMesh::queryRealm(string)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::queryRealm(string)" << endl;
  }

  bool ret = true;

  if (m_realms.find(a_realm) == m_realms.end()) {
    ret = false;
  }

  return ret;
}

bool
AmrMesh::getEbCf() const
{
  CH_TIME("AmrMesh::getEbCf()");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::getEbCf()" << endl;
  }

  return true;
}

void
AmrMesh::registerRealm(const std::string a_realm)
{
  CH_TIME("AmrMesh::registerRealm(string)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::registerRealm(string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    m_realms.emplace(a_realm, RefCountedPtr<Realm>(new Realm()));
  }
}

void
AmrMesh::registerOperator(const std::string a_operator, const std::string a_realm, const phase::which_phase a_phase)
{
  CH_TIME("AmrMesh::registerOperator(string, string, phase::which_phase)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::registerOperator(string, string, phase::which_phase)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str =
      "AmrMesh::registerOperator(string, string, phase::which_phase) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  m_realms[a_realm]->registerOperator(a_operator, a_phase);
}

void
AmrMesh::registerMask(const std::string a_mask, const int a_buffer, const std::string a_realm)
{
  CH_TIME("AmrMesh::registerMask(string, int, string)");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::registerMask(string, int, string)" << endl;
  }

  if (!this->queryRealm(a_realm)) {
    const std::string str = "AmrMesh::registerMask(string, int, string) - could not find realm '" + a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  m_realms[a_realm]->registerMask(a_mask, a_buffer);
}

void
AmrMesh::defineRealms()
{
  CH_TIME("AmrMesh::defineRealms()");
  if (m_verbosity > 5) {
    pout() << "AmrMesh::defineRealms()" << endl;
  }

  for (auto& r : m_realms) {
    r.second->define(m_grids,
                     m_domains,
                     m_refinementRatios,
                     m_dx,
                     m_probLo,
                     m_finestLevel,
                     m_numEbGhostsCells,
                     m_numGhostCells,
                     m_numLsfGhostCells,
                     m_redistributionRadius,
                     m_multigridInterpOrder,
                     m_multigridInterpRadius,
                     m_multigridInterpWeight,
                     m_centroidStencilType,
                     m_ebCentroidStencilType,
                     m_baseif,
                     m_multifluidIndexSpace);
  }
}

void
AmrMesh::regridRealm(const std::string          a_realm,
                     const Vector<Vector<int>>& a_procs,
                     const Vector<Vector<Box>>& a_boxes,
                     const int                  a_lmin)
{
  CH_TIME("AmrMesh::regridRealm(string, Vector<Vector<int> >, Vector<Vector<Box> >, int)");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::regridRealm(string, Vector<Vector<int> >, Vector<Vector<Box> >, int)" << endl;
  }

  // TLDR: This function does a base-regrid of a realm, using the specified input processor IDs and boxes.

  if (!this->queryRealm(a_realm)) {
    std::string str =
      "AmrMesh::regridRealm(string, Vector<Vector<int> >, Vector<Vector<Box> >, int) - could not find realm '" +
      a_realm + "'";
    MayDay::Abort(str.c_str());
  }

  // Make the dbl
  Vector<DisjointBoxLayout> grids(1 + m_finestLevel);

  // Levels that didn't change.
  for (int lvl = 0; lvl < a_lmin; lvl++) {
    grids[lvl] = this->getGrids(a_realm)[lvl];
  }

  // Levels that did change.
  for (int lvl = a_lmin; lvl <= m_finestLevel; lvl++) {
    grids[lvl] = DisjointBoxLayout();
    grids[lvl].define(a_boxes[lvl], a_procs[lvl], m_domains[lvl]);
    grids[lvl].close();
  }

  m_realms[a_realm]->define(grids,
                            m_domains,
                            m_refinementRatios,
                            m_dx,
                            m_probLo,
                            m_finestLevel,
                            m_numEbGhostsCells,
                            m_numGhostCells,
                            m_numLsfGhostCells,
                            m_redistributionRadius,
                            m_multigridInterpOrder,
                            m_multigridInterpRadius,
                            m_multigridInterpWeight,
                            m_centroidStencilType,
                            m_ebCentroidStencilType,
                            m_baseif,
                            m_multifluidIndexSpace);

  m_realms[a_realm]->regridBase(a_lmin);
}

std::vector<std::string>
AmrMesh::getRealms() const
{
  CH_TIME("AmrMesh::getRealms()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getRealms()" << endl;
  }

  std::vector<std::string> Realms;

  for (const auto& r : m_realms) {
    Realms.push_back(r.first);
  }

  return Realms;
}

BoxSorting
AmrMesh::getBoxSorting() const
{
  CH_TIME("AmrMesh::getBoxSorting()");
  if (m_verbosity > 1) {
    pout() << "AmrMesh::getBoxSorting()" << endl;
  }

  return m_boxSort;
}

#include <CD_NamespaceFooter.H>
