#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBEllipticLoadBalance.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBISBox.H"
#include "LoadBalance.H"
#include "EBPoissonOp.H"
#include "EBPoissonOpFactory.H"
#include "SPMD.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "EBLevelDataOps.H"
#include "NeumannPoissonDomainBC.H"
#include "NeumannPoissonEBBC.H"
#include "TimedDataIterator.H"
#include "NWOEBViscousTensorOpFactory.H"
#include "NamespaceHeader.H"




///////////////
void 
getNWOEBVTORelaxLoadsAndBoxes(Vector<unsigned long long>               &    a_loads,
                              Vector<Box>                              &    a_boxes,
                              const RefCountedPtr<BaseDomainBCFactory> &    a_domBC,
                              const RefCountedPtr<BaseEBBCFactory>     &    a_ebBC,
                              const DisjointBoxLayout                  &    a_dblOrig,
                              const ProblemDomain                      &    a_domain,
                              const EBIndexSpace                       *    a_ebisPtr )
{
  int nghost  = 4;
  EBLevelGrid eblg(a_dblOrig, a_domain, nghost, a_ebisPtr);

  //and we need to remap
  TimedDataIterator dit = a_dblOrig.timedDataIterator();
  dit.clearTime();
  dit.enableTime();

  NWOEBViscousTensorOpFactory opFact(eblg, a_domBC, a_ebBC);
  IntVect nghostPhi = opFact.getGhostPhi();
  IntVect nghostRHS = opFact.getGhostRHS();
  EBCellFactory ebcellfact(eblg.getEBISL());
  LevelData<EBCellFAB> phi(a_dblOrig, SpaceDim, nghostPhi, ebcellfact);
  LevelData<EBCellFAB> rhs(a_dblOrig, SpaceDim, nghostRHS, ebcellfact);

  EBLevelDataOps::setToZero(rhs);
  EBLevelDataOps::setToZero(phi);

  RefCountedPtr<NWOEBViscousTensorOp>  ebpo = RefCountedPtr<NWOEBViscousTensorOp>
    (opFact.newOpForTimingOnly(a_domain));

  int numRelax = 32;

  //evaluate poisson operator---homogeneous bcs so i don't have to set the value
  ebpo->relaxTimed(phi, rhs, dit, numRelax);

  Vector<Box>  boxesLocal = dit.getBoxes();
  Vector<unsigned long long> loadsLocal = dit.getTime();

  dit.disableTime();

  dit.mergeTime();
  a_loads = dit.getTime();
  a_boxes = dit.getBoxes();

}
///////////////
int
NWOEBVTORelaxLoadBalance(Vector<int>                              &  a_procs,
                         Vector<Box>                              &  a_boxes,
                         const RefCountedPtr<BaseDomainBCFactory> &  a_domBC,
                         const RefCountedPtr<BaseEBBCFactory>     &  a_ebBC,
                         const ProblemDomain                      &  a_domain,
                         const EBIndexSpace                       *  a_ebis_ptr )
{
#ifndef CH_MPI
  a_procs.resize(a_boxes.size(),0);
  int retval=0;
#else
  //first load balance the conventional way.
  Vector<Box> inBoxes = a_boxes;
  Vector<int> origProcs;
  int retval=LoadBalance(origProcs, inBoxes);

  a_procs = origProcs;

  DisjointBoxLayout dblOrig(inBoxes, origProcs, a_domain);

  Vector<unsigned long long> loads;
  getNWOEBVTORelaxLoadsAndBoxes(loads, a_boxes, a_domBC, a_ebBC, dblOrig, a_domain, a_ebis_ptr );

  retval = UnLongLongLoadBalance(a_procs,  loads, a_boxes);
#endif

  return retval;
}
#include "NamespaceFooter.H"
