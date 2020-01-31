#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKFACEREGISTER_H_
#define _MULTIBLOCKFACEREGISTER_H_

#include "LevelFluxRegister.H"
#include "MultiBlockCoordSys.H"
#include "RemoteCopier.H"
#include "NamespaceHeader.H"

//! \class MultiBlockFaceRegister
//! This subclass of LevelFluxRegister handles multi-block AMR refluxing.
class MultiBlockFaceRegister: public LevelFluxRegister
{
private:

  public:

  //! Creates a face register that handles multi-block refluxing.
  //! This capability requires a multi-block coordinate system, which
  //! handles block connectivity as well as coordinate mappings.
  MultiBlockFaceRegister(const MultiBlockCoordSys *const a_coordSys);

  //! Destructor.
  ~MultiBlockFaceRegister();

  /// Define everything (overrides base)
  void define(const DisjointBoxLayout& a_FnGrid,
              const DisjointBoxLayout& a_CrGrid,
              const ProblemDomain&     a_FnProbDom,
              int                      a_nRefine,
              int                      a_nComp,
              bool                     a_scaleFineFaces = true);

  /// Define everything (overrides base)
  void define(const DisjointBoxLayout& a_FnGrid,
              const DisjointBoxLayout& a_CrGrid,
              const ProblemDomain&     a_FnProbDom,
              IntVect                      a_ivRefine,
              int                      a_nComp,
              bool                     a_scaleFineFaces = true);


  /// Define the fine register
  void defineFine(const DisjointBoxLayout& a_FnGrid,
                  const ProblemDomain&     a_FnProbDom,
                  int                      a_nRefine,
                  int                      a_nComp,
                  bool                     a_scaleFineFaces = true);

  /// Define the fine register
  void defineFine(const DisjointBoxLayout& a_FnGrid,
                  const ProblemDomain&     a_FnProbDom,
                  IntVect                  a_ivRefine,
                  int                      a_nComp,
                  bool                     a_scaleFineFaces = true);


  /// Define the coarse register and copiers
  void defineCoarse(const DisjointBoxLayout& a_CrGrid);

  // New reface method (analagous to reflux)
  void reface(LevelData<FluxBox>& a_uCoarse,
              const Interval&     a_coarse_interval,
              const Interval&     a_face_interval,
              const Real          a_scale,
              const bool          a_directional = true);

  private:

  //! Multi-block coordinate system.
  const MultiBlockCoordSys *const m_coordSys;
                                      ///< Coordinate system from the coarse
                                      ///< level

  //! Copier for moving values from fine faces on other blocks to
  //! blocks where the corresponding coarse faces are defined.
  //! One copier per dimension.
  RemoteCopier m_remoteCopiers[SpaceDim][2];
};
#include "NamespaceFooter.H"
#endif
