#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SIMPLEDIVERGENCE_H_
#define _SIMPLEDIVERGENCE_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "RealVect.H"

#include "NamespaceHeader.H"

///
/**
   In a_divF, on a_box, stores divergence of a_F,
   where a_F is defined on a_box.
 */
void simpleDivergence(FArrayBox& a_divF,
                      const FluxBox& a_F,
                      const Box& a_box,
                      RealVect a_dx);


#include "NamespaceFooter.H"
#endif




