#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SimpleDivergence.H"
#include "FCDivergenceF_F.H"


#include "NamespaceHeader.H"

void simpleDivergence(FArrayBox& a_divF,
                      const FluxBox& a_F,
                      const Box& a_box,
                      RealVect a_dx)
{

  a_divF.setVal(0.0, a_box, 0, a_divF.nComp());

  for (int dir=0; dir<SpaceDim; dir++)
    {
      FORT_FCDIVERGENCE(CHF_CONST_FRA(a_F[dir]),
                        CHF_FRA(a_divF),
                        CHF_BOX(a_box),
                        CHF_CONST_REAL(a_dx[dir]),
                        CHF_INT(dir));
    }
}




#include "NamespaceFooter.H"
