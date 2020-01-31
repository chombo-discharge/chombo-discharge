#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "BRMeshRefine.H"
#include "LevelOp.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "AMRIO.H"
#include "IntVectSet.H"

#include "localfunctions.H"
#include "fortranfunctions.H"
#include "DebugOut.H"

// under construction -- not ready for primetime ndk


// poor dood's search
bool findValue(const LevelData<FArrayBox>& locLDFB,
               const int x, const int y, const int z, double& value)
{
  DataIterator dit = locLDFB.dataIterator();

  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& fab = locLDFB[dit()];
    BoxIterator bit(fab.box());
    for (bit.reset(); bit.ok(); ++bit)
    {
      int xx=bit().getVect()[0];
      int yy=bit().getVect()[1];
#if CH_SPACEDIM == 3
      int zz=bit().getVect()[2];
      if (x==xx && y==yy && z==zz)
      {
#else
      if (x==xx && y==yy)
      {
#endif
        value = fab(bit(),0);
        return true;
      }
    }
  }
  return false;
}


void
prettyPrintLDFB(const LevelData<FArrayBox>*  memLDFB_Ptr)
{
  const string blank;
  prettyPrintLDFB(memLDFB_Ptr, blank);
}

void
prettyPrintLDFB(const LevelData<FArrayBox>*  memLDFB_Ptr,
               const string& label)
{

  const LevelData<FArrayBox>& memLDFB = *memLDFB_Ptr;
  Vector<Box> boxes;
  const DisjointBoxLayout& memDBL = memLDFB.getBoxes();
  LayoutIterator lit = memDBL.layoutIterator();
  for (lit.reset(); lit.ok(); ++lit)
  {
    boxes.push_back(memDBL.get(lit()));
  }

  Vector<int> assign(boxes.size(), 0);
  DisjointBoxLayout locDBL(boxes, assign);
  locDBL.close();
  LevelData<FArrayBox> locLDFB(locDBL,
                              memLDFB.nComp(),
                              memLDFB.ghostVect());
  const Interval& interv = memLDFB.interval();
  memLDFB.copyTo(interv, locLDFB, interv);


  LayoutIterator it = memDBL.layoutIterator();
  Box bound = memDBL[it()]; ++it;
  for (;it.ok() ; ++it)
  {
    bound = minBox(bound, memDBL[it()]);
  }
  const int xmin =  bound.smallEnd(0);
  const int xmax =  bound.bigEnd(0);
  const int ymin =  bound.smallEnd(1);
  const int ymax =  bound.bigEnd(1);

  double value;

  printf(" %-40s ----------------------------------------\n", label.c_str());
  int k=0;
  for (int i=xmin; i<=xmax; i++)
  {
    printf("%4d ", i);
    for (int j=ymin; j<=ymax; j++)
    {
      if (findValue(locLDFB, i, j, k, value))
      {
        printf("%11.3e ", value);
      }
      else
      {
        printf("     .      ");
      }
    }
    printf("\n");
  }

  //printf("%4d-%11.3e", x, fab(bit(),0));

}

