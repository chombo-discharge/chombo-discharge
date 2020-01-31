#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FindFace.H"

#include "NamespaceHeader.H"

FindFace::FindFace(const Box& a_baseBox,
                   const Tuple<Box, SpaceDim>& a_faces)
{
  m_faces = a_faces; 
  m_baseBox = a_baseBox;
}

FindFace::~FindFace()
{
}

Box
FindFace::operator() (const Box& a_inputBox)
{
  // yes, this is clumsy
  bool found = false;
  int idir;
  for (idir = 0; idir < SpaceDim; idir++)
    {
      Box loFace = m_baseBox;
      loFace.setBig(idir, m_baseBox.smallEnd(idir));
      if (loFace == a_inputBox)
        {
          found = true;
          break;
        }

      Box hiFace = m_baseBox;
      hiFace.setSmall(idir, m_baseBox.bigEnd(idir));
      if (hiFace == a_inputBox)
        {
          found = true;
          break;
        }
    }
  CH_assert(found);
  Box returnBox = m_faces[idir];
  return returnBox;
}

#include "NamespaceFooter.H"
