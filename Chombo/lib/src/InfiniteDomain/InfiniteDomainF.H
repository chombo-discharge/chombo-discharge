#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INFINITEDOMAINF_H_
#define _INFINITEDOMAINF_H_

#include "FORT_PROTO.H"

#include "NamespaceHeader.H"

extern "C"
{
#define FORT_DIFFPHI FORTRAN_NAME(DIFFPHI, diffphi)
  void FORT_DIFFPHI
  (
   Real* phiex,
   Real* phi,
   int* n,
   int* layer,
   Real* h,
   int* vrbos,
   Real* difnrm
   );


#define FORT_PRECONDRHS FORTRAN_NAME(PRECONDRHS, precondrhs)
  void FORT_PRECONDRHS
  (
   Real* rhs,
   const int* l1, const int* l2, const int* l3,
   const int* u1, const int* u2, const int* u3
   );


#define FORT_INHOMODIRICHLET3D FORTRAN_NAME(INHOMODIRICHLET3D, inhomodirichlet3d)
  void FORT_INHOMODIRICHLET3D
  (
   Real* f,
   const int* lx, const int* ly, const int* lz,
   const int* ux, const int* uy, const int* uz,
   const Real* dx, const Real* dy, const Real* dz,
   const int* op,
   int* ier
   );


#define FORT_POISSON3D FORTRAN_NAME(POISSON3D, poisson3d)
  void FORT_POISSON3D
  (
   Real* f,
   int* n1,
   int* n2,
   int* n3,
   Real* h,
   int* op
   );


#define FORT_INTERPFACE FORTRAN_NAME(INTERPFACE, interpface)
  void FORT_INTERPFACE
  (
   Real* fine,
   Real* coarse,
   Real* intrpc,
   int* n1,
   int* n2,
   int* refToCoarse,
   int* n1Coarse,
   int* n2Coarse,
   int* layer
   );


#define FORT_GETPOLYINTERPC FORTRAN_NAME(GETPOLYINTERPC, getpolyinterpc)
  void FORT_GETPOLYINTERPC
  (
   Real* intrpc,
   int* refToCoarse,
   int* layer
   );


#define FORT_GETBUFFERSIZE FORTRAN_NAME(GETBUFFERSIZE, getbuffersize)
  void FORT_GETBUFFERSIZE
  (
   int* domainLength,
   int* bufferSizeMin,
   int* bufferSize,
   int* patchSize,
   int* errCode
   );


#define FORT_GET0BUFFERSIZE FORTRAN_NAME(GET0BUFFERSIZE, get0buffersize)
  void FORT_GET0BUFFERSIZE
  (
   int* domainLength,
   int* bufferSizeMin,
   int* bufferSize,
   int* patchSize,
   int* errCode
   );


#define FORT_GET0BUFFERSIZE3D FORTRAN_NAME(GET0BUFFERSIZE3D, get0buffersize3d)
  void FORT_GET0BUFFERSIZE3D
  (
   int* domainLength,
   int* bufferSizeMin,
   int* bufferSize,
   int* patchSize,
   int* errCode
   );


#define FORT_GETBUFFERSIZE3D FORTRAN_NAME(GETBUFFERSIZE3D, getbuffersize3d)
  void FORT_GETBUFFERSIZE3D
  (
   int* domainLength,
   int* bufferSizeMin,
   int* bufferSize,
   int* patchSize,
   int* errCode
   );


#define FORT_GETCBUFFERSIZE FORTRAN_NAME(GETCBUFFERSIZE, getcbuffersize)
  void FORT_GETCBUFFERSIZE
  (
   int* domainLength,
   int* bufferSizeMin,
   int* bufferSize,
   int* patchSize,
   int* errCode
   );


#define FORT_COEFFS1SIDEDDERIV FORTRAN_NAME(COEFFS1SIDEDDERIV, coeffs1sidedderiv)
  void FORT_COEFFS1SIDEDDERIV
  (
   Real* coeff,
   int* d,
   Real* h
   );


#define FORT_COEFFSDERIV FORTRAN_NAME(COEFFSDERIV, coeffsderiv)
  void FORT_COEFFSDERIV
  (
   Real* coeff,
   int* d,
   Real* h,
   Real* start,
   Real* evalpt
   );


#define FORT_DERIV1OVERR FORTRAN_NAME(DERIV1OVERR, deriv1overr)
  void FORT_DERIV1OVERR
  (
   int* cpind,
   Real* cfac,
   int* cpow,
   int* nterms,
   int* deg
   );


#define FORT_GETRECTMATRIX FORTRAN_NAME(GETRECTMATRIX, getrectmatrix)
  void FORT_GETRECTMATRIX
  (
   Real* srcdst,
   const int* iev,
   const int* icp,
   const int* nsrc,
   const int* ndst,
   const int* refc,
   const int* ptchsz,
   const int* npchs1,
   const int* npchs2,
   const int* ndst1c,
   const int* ndst2c,
   const int* deg,
   const int* layer,
   const Real* h,
   const int* cpind,
   const Real* cfac,
   const int* cpow,
   const int* nterms,
   const int* vrbos,
   const int* ier
   );


#define FORT_GETORIRECTMATRIX FORTRAN_NAME(GETORIRECTMATRIX, getorirectmatrix)
  void FORT_GETORIRECTMATRIX
  (
   const int* ori,
   Real* srcdst,
   const int* iev,
   const int* icp,
   const int* nsrc,
   const int* ndst,
   const int* refc,
   const int* ptchsz,
   const int* npchs1,
   const int* npchs2,
   const int* ndst1c,
   const int* ndst2c,
   const int* deg,
   const int* layer,
   const Real* h,
   const int* cpind,
   const Real* cfac,
   const int* cpow,
   const int* nterms,
   const int* vrbos,
   const int* ier
   );


#define FORT_GETDIRECTRECTMATRIX FORTRAN_NAME(GETDIRECTRECTMATRIX, getdirectrectmatrix)
  void FORT_GETDIRECTRECTMATRIX
  (
   Real* srcdst,
   const Real* wtsint,
   const int* iev,
   const int* icp,
   const int* nsrc,
   const int* ndst1c,
   const int* ndst2c,
   const int* layer,
   const Real* hsrc,
   const Real* hdst,
   const Real* cpbase,
   const Real* evbase,
   const int* vrbos,
   const int* ier
   );


#define FORT_GETALLCOARSEMATRIX FORTRAN_NAME(GETALLCOARSEMATRIX, getallcoarsematrix)
  void FORT_GETALLCOARSEMATRIX
  (
   Real* sdcrse,
   int* ncpts,
   int* n1,
   int* s2,
   int* ptchsz,
   int* npchs,
   int* deg,
   Real* h,
   int* refc,
   int* cadd,
   int* cpind,
   Real* cfac,
   int* cpow,
   int* nterms,
   int* vrbos,
   int* ier
   );


#define FORT_GETRECTCOARSEMATRIX FORTRAN_NAME(GETRECTCOARSEMATRIX, getrectcoarsematrix)
  void FORT_GETRECTCOARSEMATRIX
  (
   Real* sdcrse,
   int* icp,
   int* nsrc,
   int* ndst,
   int* refc,
   int* ptchsz,
   int* npchs1,
   int* npchs2,
   int* deg,
   Real* h,
   int* cadd,
   int* cpind,
   Real* cfac,
   int* cpow,
   int* nterms,
   int* ncpts,
   int* vrbos,
   int* ier
   );


#define FORT_ADDEVAL1FACECOARSE FORTRAN_NAME(ADDEVAL1FACECOARSE, addeval1facecoarse)
  void FORT_ADDEVAL1FACECOARSE
  (
   Real* ev,
   const Real* cp,
   const int* npchs1,
   const int* npchs2,
   const int* ndst1c,
   const int* ndst2c,
   const int* nc,
   const int* layer,
   const Real* srcdst,
   const int* vrbos,
   int* ier
   );


#define FORT_ADDEVALSYM1FACECOARSE FORTRAN_NAME(ADDEVALSYM1FACECOARSE, addevalsym1facecoarse)
  void FORT_ADDEVALSYM1FACECOARSE
  (
   Real* ev,
   const Real* cp,
   const int* dflip,
   const int* sflip,
   const int* npchs1,
   const int* npchs2,
   const int* ndst1c,
   const int* ndst2c,
   const int* deg,
   const int* nc,
   const int* layer,
   const Real* srcdst,
   const int* vrbos,
   int* ier
   );


#define FORT_ADDEVAL1DIRECTFACECOARSE FORTRAN_NAME(ADDEVAL1DIRECTFACECOARSE, addeval1directfacecoarse)
  void FORT_ADDEVAL1DIRECTFACECOARSE
  (
   Real* ev,
   const Real* cp,
   const int* nsrc,
   const int* ndst1c,
   const int* ndst2c,
   const int* layer,
   const Real* srcdst,
   const int* vrbos,
   int* ier
   );


#define FORT_EVALOUTERCOARSE FORTRAN_NAME(EVALOUTERCOARSE, evaloutercoarse)
  void FORT_EVALOUTERCOARSE
  (
   Real* ev,
   const Real* cp,
   const int* ncpts,
   const int* npchs,
   const int* n2c,
   const int* cadd,
   const int* nc,
   const int* ori,
   const Real* sdcrse,
   const int* vrbos,
   int* ier
   );


#define FORT_EVALRECTOUTERCOARSE FORTRAN_NAME(EVALRECTOUTERCOARSE, evalrectoutercoarse)
  void FORT_EVALRECTOUTERCOARSE
  (
   Real* ev,
   const Real* cp,
   const int* ncpts,
   const int* npchs1,
   const int* npchs2,
   const int* ndstc,
   const int* cadd,
   const int* nc,
   const Real* sdcrse,
   const int* vrbos,
   int* ier
   );


#define FORT_EVALSYMRECTOUTERCOARSE FORTRAN_NAME(EVALSYMRECTOUTERCOARSE, evalsymrectoutercoarse)
  void FORT_EVALSYMRECTOUTERCOARSE
  (
   Real* ev,
   const Real* cp,
   const int* dflip,
   const int* ncpts,
   const int* npchs1,
   const int* npchs2,
   const int* ndstc,
   const int* cadd,
   const int* nc,
   const Real* sdcrse,
   const int* vrbos,
   int* ier
   );


#define FORT_EVALDIRECTRECTOUTERCOARSE FORTRAN_NAME(EVALDIRECTRECTOUTERCOARSE, evaldirectrectoutercoarse)
  void FORT_EVALDIRECTRECTOUTERCOARSE
  (
   Real* ev,
   const Real* cp,
   const int* ncpts,
   const int* ndstc,
   const int* cadd,
   const int* nsrc,
   const Real* sdcrse,
   const int* vrbos,
   int* ier
   );

#define FORT_GETRECTDIRECTCOARSEMATRIX FORTRAN_NAME(GETRECTDIRECTCOARSEMATRIX, getrectdirectcoarsematrix)
  void FORT_GETRECTDIRECTCOARSEMATRIX
  (
   Real* sdcrse,
   const Real* wtsint,
   const int* ncpts,
   const int* ndstc,
   const int* cadd,
   const int* nsrc,
   const Real* hsrc,
   const Real* hdst,
   const Real* cpbase,
   const Real* evbase,
   const int* vrbos,
   const int* ier
   );

}

#include "NamespaceFooter.H"

#endif
