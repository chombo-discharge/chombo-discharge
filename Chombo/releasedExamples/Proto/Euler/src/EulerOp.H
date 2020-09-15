#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EULEROP_H_
#define _EULEROP_H_

#pragma once

#include "Proto.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "FArrayBox.H"
#include "Copier.H"

constexpr unsigned int NUMCOMPS=DIM+2;
#define NGHOST 4
typedef Proto::Box Bx;
using   Proto::Point;
using   Proto::BoxData;
using     std::cout;
using     std::endl;
using     std::shared_ptr;
using   Proto::Stencil;

class EulerOp
{

public:

  static double s_dx;
  static double s_gamma;

  /// Compute the Right-Hand side, conforming to the RK4 interface.
  static double step(LevelData<FArrayBox> & a_Rhs,
                     LevelData<FArrayBox> & a_U);

  //for initial time step
  static double maxWave(LevelData<FArrayBox> & a_U);

  EulerOp(){};
  ~EulerOp(){};

private:
  // all the messy internals.
  static double proto_step(BoxData<double,NUMCOMPS>& a_Rhs,
                           const BoxData<double,NUMCOMPS>& a_U,
                           const Bx& a_dbx0);

  /// called internally from step
  static void define(const DisjointBoxLayout& a_grids,
                     const IntVect          & a_ghostVect) ;
  

  static Stencil<double> s_laplacian;
  static Stencil<double> s_deconvolve;
  static Stencil<double> s_laplacian_f[DIM];
  static Stencil<double> s_deconvolve_f[DIM];
  static Stencil<double> s_interp_H[DIM];
  static Stencil<double> s_interp_L[DIM];
  static Stencil<double> s_divergence[DIM];
  static Copier          s_exchangeCopier;

};

#endif
