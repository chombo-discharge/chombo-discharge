#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef CH_SPACEDIM

#include "NodeDatasetClient.H"
#include "LayoutIterator.H"
#include "DataIterator.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

const char *defaultIntFormat = "%7.0f";

bool
ArrayViewNode(NodeFArrayBox* a_nfab)
{
  FArrayBox& fab = a_nfab->getFab();
  return ArrayView(&fab);
}


bool
ArrayViewBool(BaseFab<bool>* a_nfab)
{
  Box bx = a_nfab->box();
  BaseFab<int> intFab(bx, a_nfab->nComp());
  for (BoxIterator bit(bx); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      for (int nv = 0; nv < intFab.nComp(); nv++)
        intFab(iv, nv) = a_nfab->operator()(iv, nv) ? 1 : 0;
    }
  return ArrayViewInt(&intFab);
}


bool
MultiArrayViewLevelNodeFab(LevelData<NodeFArrayBox>* a_nfabArrayInPtr)
{
  const BoxLayout& boxesIn = a_nfabArrayInPtr->boxLayout();
  const int nvar = a_nfabArrayInPtr->nComp();
  const IntVect ghosts = a_nfabArrayInPtr->ghostVect();

  // Set boxes of boxesTemp to be nodes surrounding boxes of boxesIn.
  BoxLayout boxesTemp;
  boxesTemp.deepCopy(boxesIn);
  for (LayoutIterator lit = boxesIn.layoutIterator(); lit.ok(); ++lit)
    {
      Box nodeBox = boxesIn.get(lit());
      nodeBox.surroundingNodes();
      nodeBox.grow(ghosts);
      boxesTemp.ref(lit()) = nodeBox;
    }
  boxesTemp.close();

  // Copy a_nfabArrayInPtr to arrayTemp.
  BoxLayoutData<FArrayBox> arrayTemp(boxesTemp, nvar);
  for (DataIterator dit = a_nfabArrayInPtr->dataIterator(); dit.ok(); ++dit)
    {
      arrayTemp[dit()].copy(a_nfabArrayInPtr->operator[](dit()).getFab());
    }

  return MultiArrayViewFab(&arrayTemp);
}



bool
MultiArrayViewNodeFab(BoxLayoutData<NodeFArrayBox>* a_nfabArrayInPtr)
{
  const BoxLayout& boxesIn = a_nfabArrayInPtr->boxLayout();
  const int nvar = a_nfabArrayInPtr->nComp();

  // Set boxes of boxesTemp to be nodes surrounding boxes of boxesIn.
  BoxLayout boxesTemp;
  boxesTemp.deepCopy(boxesIn);
  for (LayoutIterator lit = boxesIn.layoutIterator(); lit.ok(); ++lit)
    {
      Box nodeBox = boxesIn.get(lit());
      nodeBox.surroundingNodes();
      boxesTemp.ref(lit()) = nodeBox;
    }
  boxesTemp.close();

  // Copy a_nfabArrayInPtr to arrayTemp.
  BoxLayoutData<FArrayBox> arrayTemp(boxesTemp, nvar);
  for (DataIterator dit = a_nfabArrayInPtr->dataIterator(); dit.ok(); ++dit)
    {
      arrayTemp[dit()].copy(a_nfabArrayInPtr->operator[](dit()).getFab());
    }

  return MultiArrayViewFab(&arrayTemp);
}



bool
MultiArrayViewInt(BoxLayoutData< BaseFab<int> >* debugLayoutData)
{
  const BoxLayout layout = debugLayoutData->boxLayout();
  LayoutData< BaseFab<Real> >* realLayoutData;
  realLayoutData = new BoxLayoutData< BaseFab<Real> >(layout, debugLayoutData->nComp());
  // realLayoutData->define(layout);
  for (DataIterator dit = debugLayoutData->dataIterator(); dit.ok(); ++dit)
    {
      BaseFab<int>& intFab = debugLayoutData->operator[](dit());
      BaseFab<Real>& realFab = realLayoutData->operator[](dit());
      for (BoxIterator bit(realFab.box()); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          for (int nv = 0; nv < intFab.nComp(); nv++)
            realFab(iv, nv) = intFab(iv, nv) * 1.0;
        }
    }
  //  return ( MultiArrayView(realLayoutData) );
  return (MultiArrayViewFormatLabel(realLayoutData,
                                    defaultIntFormat,
                                    "LayoutData<BaseFab<Int>>"));
}
#include "NamespaceFooter.H"

#endif // CH_SPACEDIM
