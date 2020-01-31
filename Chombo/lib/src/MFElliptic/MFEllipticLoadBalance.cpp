#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MFEllipticLoadBalance.H"
#include "EBEllipticLoadBalance.H"
#include "MFIndexSpace.H"
#include "EBISLayout.H"
#include "EBISBox.H"
#include "LoadBalance.H"
#include "MFPoissonOp.H"
#include "MFPoissonOpFactory.H"
#include "SPMD.H"
#include "MFCellFAB.H"
#include "MFLevelDataOps.H"
#include "NeumannPoissonDomainBC.H"
#include "NeumannPoissonEBBC.H"
#include "TimedDataIterator.H"
#include "NamespaceHeader.H"

#define MFELB_NUM_APPLY_OPS 100
///////////////
void getPoissonLoadsAndBoxes(Vector<unsigned long long>&   a_loads,
                             Vector<Box> &                 a_boxes,
                             RefCountedPtr<MFIndexSpace>&  a_mfis,
                             const DisjointBoxLayout&      a_dblOrig,
                             const ProblemDomain&          a_domain)
{
  int nphase = a_mfis->numPhases();
  int nghost = 4;
  IntVect nghostPhi = 4*IntVect::Unit;
  Vector<EBISLayout> ebislVec;
  ebislVec.resize(nphase);
  for (int iphase=0; iphase<nphase; iphase++)
    {
      const EBIndexSpace* const ebisPtr = a_mfis->EBIS(iphase);
      ebisPtr->fillEBISLayout(ebislVec[iphase], a_dblOrig, a_domain, nghost);
    }

  //and we need to remap
  TimedDataIterator dit = a_dblOrig.timedDataIterator();
  dit.clearTime();
  dit.enableTime();
  //scoping bracket
  Vector<int> nComp(nphase, 1);
  MFCellFactory mfcellfact(ebislVec, nComp);
  LevelData<MFCellFAB> phi(a_dblOrig, 1, nghostPhi, mfcellfact);
  LevelData<MFCellFAB> lph(a_dblOrig, 1, nghostPhi, mfcellfact);

  MFLevelDataOps::setToZero(lph);
  MFLevelDataOps::setVal(phi, 1.0);

  Vector< RefCountedPtr<BaseDomainBC> > domBC(nphase);
  NeumannPoissonDomainBCFactory neumannBCFactory;
  neumannBCFactory.setValue(0.0);
  RealVect dx = RealVect::Unit;
  for (int iphase=0; iphase<nphase; iphase++)
    {
      domBC[iphase] = RefCountedPtr<BaseDomainBC>
        (neumannBCFactory.create(a_domain, ebislVec[iphase], dx));
    }

  Vector<DisjointBoxLayout> dblVec(1, a_dblOrig);
  Vector<int> refRatio(1, 2);
  Vector<Real> alpha(nphase, -1.0);
  Vector<Real> beta(nphase, 1.0);

  MFPoissonOpFactory factory(a_mfis,
                             dblVec,
                             refRatio,
                             a_domain,
                             dx,
                             RealVect::Zero,
                             domBC,
                             alpha,
                             beta,
                             1,
                             nghostPhi,
                             nghostPhi);

  factory.setJump(0.0, 0.0);

  RefCountedPtr<MFPoissonOp>  mfpo = RefCountedPtr<MFPoissonOp>
    ((MFPoissonOp*) factory.MGnewOp(a_domain, 0, false));

  //evaluate poisson operator---homogeneous bcs so i don't have to set the value
  for (int iapply = 0; iapply < MFELB_NUM_APPLY_OPS; iapply++)
    {
      mfpo->applyOp(lph, phi, dit, true);
    }
  Vector<Box>  boxesLocal = dit.getBoxes();
  Vector<unsigned long long> loadsLocal = dit.getTime();

  dit.disableTime();

  dit.mergeTime();
  a_loads = dit.getTime();
  a_boxes = dit.getBoxes();

}
///////////////
int
MFEllipticLoadBalance(Vector<int>&                  a_procs,
                      const Vector<Box>&            a_boxes,
                      RefCountedPtr<MFIndexSpace>&  a_mfis,
                      const ProblemDomain&          a_domain,
                      bool                          a_verbose)
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

  //we shall make fully covered boxes = constant load = covered load
  //we shall say that irregular points get irregular factor more load than
  //regular points.
  //compute a load for each box.   by evaluating the poisson operator
  //use one scratch space for everything
  DisjointBoxLayout dblOrig(inBoxes, origProcs, a_domain);

  Vector<unsigned long long> loads;
  Vector<Box>  boxes;
  getPoissonLoadsAndBoxes(loads, boxes, a_mfis, dblOrig, a_domain);

  resetLoadOrder(loads, boxes, inBoxes);
  if (a_verbose)
    {
      pout() << "MFEllipticLoadBalance loads:" << endl;
      for (int ibox = 0; ibox < a_boxes.size(); ibox++)
        {
          pout()
            << "   box[" << ibox << "]=" <<   a_boxes[ibox]
            << ", load[" << ibox << "]=" <<     loads[ibox] << endl;
        }
    }
  //do the load balance with our EB load estimates and the original a_boxes vector
  retval = UnLongLongLoadBalance(a_procs,  loads, a_boxes);
#endif

  return retval;
}
#include "NamespaceFooter.H"
