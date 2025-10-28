/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LeastSquaresImplem.cpp
  @brief  Implementation of CD_LeastSquares.H
  @author Robert Marskar
*/

#ifndef CD_LeastSquaresImplem_H
#define CD_LeastSquaresImplem_H

// Our includes
#include <CD_LaPackUtils.H>
#include <CD_LeastSquares.H>
#include <CD_MultiIndex.H>
#include <CD_NamespaceHeader.H>

inline Vector<Real>
LeastSquares::makeDiagWeights(const Vector<RealVect>& a_displacements, const int a_pow)
{
  Vector<Real> ret(a_displacements.size(), 1.0);

  if (a_pow > 0) {
    for (int i = 0; i < a_displacements.size(); i++) {
      const RealVect& d   = a_displacements[i];
      const Real      len = d.vectorLength();
      const Real      w   = std::pow(len, -a_pow);

      ret[i] = w;
    }
  }

  return ret;
}

template <typename T>
std::map<IntVect, std::pair<VoFStencil, VoFStencil>>
LeastSquares::computeDualLevelStencils(const IntVectSet&       a_derivs,
                                       const IntVectSet&       a_knownTerms,
                                       const Vector<VolIndex>& a_fineVofs,
                                       const Vector<VolIndex>& a_coarVofs,
                                       const Vector<RealVect>& a_fineDisplacements,
                                       const Vector<RealVect>& a_coarDisplacements,
                                       const int               a_p,
                                       const int               a_order)
{
  const Vector<Real> fineWeights = LeastSquares::makeDiagWeights(a_fineDisplacements, a_p);
  const Vector<Real> coarWeights = LeastSquares::makeDiagWeights(a_coarDisplacements, a_p);

  return LeastSquares::computeDualLevelStencils<T>(a_derivs,
                                                   a_knownTerms,
                                                   a_fineVofs,
                                                   a_coarVofs,
                                                   a_fineDisplacements,
                                                   a_coarDisplacements,
                                                   fineWeights,
                                                   coarWeights,
                                                   a_order);
}

template <typename T>
std::map<IntVect, std::pair<VoFStencil, VoFStencil>>
LeastSquares::computeDualLevelStencils(const IntVectSet&       a_derivs,
                                       const IntVectSet&       a_knownTerms,
                                       const Vector<VolIndex>& a_fineVofs,
                                       const Vector<VolIndex>& a_coarVofs,
                                       const Vector<RealVect>& a_fineDisplacements,
                                       const Vector<RealVect>& a_coarDisplacements,
                                       const Vector<Real>&     a_fineWeights,
                                       const Vector<Real>&     a_coarWeights,
                                       const int               a_order)
{
  // TLDR: This routine does a two-level least squares solve of a system of equations, each equation describing
  //       a Taylor series expansion up to specified order (a_order). This routine "solves" this overdetermined system
  //       using weighted least squares so that we can obtain terms in the Taylor series. The user will have specified
  //       the desired unknowns (a_derivs) to be returned in stencil form as well as any terms that were already
  //       known (a_knownTerms).
  //
  //       The least squares solve is done by LaPack, and so the least squares system matrix is filled in Fortran
  //       order. The actual "solve" consists of using the Moore-Penrose pseudoinverse and once we have inverted
  //       the system we can form stencils for each of the specified terms in the Taylor series.
  //
  //       This routine does not contain much code, but can be difficult to understand.

  CH_assert(a_order > 0);
  CH_assert(a_fineVofs.size() == a_fineDisplacements.size());
  CH_assert(a_coarVofs.size() == a_coarDisplacements.size());
  CH_assert(a_fineVofs.size() == a_fineWeights.size());
  CH_assert(a_coarVofs.size() == a_coarWeights.size());

  // Initialize return stuff
  std::map<IntVect, std::pair<VoFStencil, VoFStencil>> ret;
  for (IVSIterator ivsIt(a_derivs); ivsIt.ok(); ++ivsIt) {
    ret.emplace(ivsIt(), std::make_pair(VoFStencil(), VoFStencil()));
  }

  if (a_derivs.numPts() > 0) {

    // This is because some unknowns (rows) can be eliminated.
    const int M     = LeastSquares::getTaylorExpansionSize(a_order) - a_knownTerms.numPts();
    const int Kfine = a_fineDisplacements.size();
    const int Kcoar = a_coarDisplacements.size();
    const int K     = Kfine + Kcoar;

    const IntVectSet isect = a_derivs & a_knownTerms;

    if (K < M) {
      MayDay::Abort("LeastSquares::computeDualLevelStencils -- not enough equations to achieve desired order!");
    }
    if (!isect.isEmpty()) {
      MayDay::Abort(
        "LeastSquares::computeDualLevelStencils - you have specified the same terms as both unknown and known");
    }

    // Build the A-matrix in column major order (this is what Fortran wants) so we can use LaPackUtils::computePseudoInverse.
    // ----------------------------------------------------------------------------------------------------------------------
    // If we have an (unweighted) full system then our system A*x = b is
    //
    //              A                   x             b
    //             ===                 ===           ===
    //     [1 (x-x0) (x-x0)^2 ...] [f(x)     ]     [f(x0)]
    //     [1 (x-x1) (x-x1)^2 ...] [df/dx    ]     [f(x1)]
    //     [1 (x-x2) (x-x2)^2 ...] [d^2f/dx^2]  =  [f(x2)]
    //     [:   :      :         ] [    :    ]     [  :  ]
    //     [:   :      :         ] [    :    ]     [  :  ]
    //
    // Extensions to 2D/3D simply use multi-index notation, and weights are simply multiplied into each row, i.e. we
    // solve (w*A) * x = (w*b). Inverting the system gives x = [w * A^+ * w] * b where A^+ is the Moore-Penrose
    // pseudoinverse. We put the result of [w * A^+ * w] into a stencil.
    //
    // Note that columns can be eliminated through knownTerms, in which case we remove unknowns (i.e., rows) from the system.
    // This will also correspond to a modification of the right-hand side, but the required modifications are not accesible
    // in this routine, and so the user will have to make sense of them.

    int       i = 0;                // Exists just because we fill memory linearly.
    Vector<T> linA(K * M, 0.0);     // Equal to (w*A)
    Vector<T> linAplus(M * K, 0.0); // Equal to (w*A)^+

    // Loop over column
    for (MultiIndex mi(a_order); mi.ok(); ++mi) {
      if (!a_knownTerms.contains(mi.getCurrentIndex())) {

        // Fill column, write the system using the fine vofs first, then the coarse vofs.
        for (int k = 0; k < K; k++) {
          if (k < Kfine) {
            linA[i] = a_fineWeights[k] * mi.pow(a_fineDisplacements[k]) / mi.factorial();
          }
          else {
            linA[i] = a_coarWeights[k - Kfine] * mi.pow(a_coarDisplacements[k - Kfine]) / mi.factorial();
          }
          i++;
        }
      }
    }

    // Compute the pseudo-inverse.
    const bool foundSVD = LaPackUtils::computePseudoInverse(linAplus.stdVector(), linA.stdVector(), K, M);

    if (foundSVD) {
      // When we have eliminated rows in the linear system we can't use MultiIndex to map directly, this occurs because
      // if we eliminated unknowns (rows) our system can be something like (if we eliminated term f(x)):
      //
      //     [(x-x0) (x-x0)^2 ...] [df/dx    ]     [f(x0) - f(x)]
      //     [(x-x1) (x-x1)^2 ...] [d^2f/dx^2]     [f(x1) - f(x)]
      //     [(x-x2) (x-x2)^2 ...] [         ]  =  [f(x2) - f(x)]
      //     [  :      :         ] [    :    ]     [  :         ]
      //     [  :      :         ] [    :    ]     [  :         ]
      //
      // The MultiIndex won't know about this, and will (correctly!), believe that the first row corresponds to
      // multi-index (0,0,0). Since that is truly not the case, we map the rows in A to multi-indices and use that to
      // identify the row in (w*A)^+ that corresponds to a specific unknown in the Taylor series. This is what happens below.
      std::map<IntVect, int> rowMap;
      int                    row = 0;
      for (MultiIndex mi(a_order); mi.ok(); ++mi) {

        // This is the order in which A was built.
        if (!a_knownTerms.contains(mi.getCurrentIndex())) {
          rowMap.emplace(mi.getCurrentIndex(), row);
          row++;
        }
      }

      // Recall that linAplus is M*K so the stride is always M, starting at some specified row.
      for (IVSIterator ivsIt(a_derivs); ivsIt.ok(); ++ivsIt) {
        const IntVect deriv = ivsIt();

        if (rowMap.find(ivsIt()) != rowMap.end()) {
          row = rowMap.at(ivsIt());
        }
        else {
          MayDay::Error("LeastSquares::computeDualLevelStencils -- map is out of range but this shouldn't happen!");
        }

        std::pair<VoFStencil, VoFStencil>& sten = ret.at(deriv);

        sten.first.clear();
        sten.second.clear();

        // Map the pseudoinverse into something that is usable by a stencil. Note that linAplus is (w*A)^+, but we want
        // the term [(w*A)^+ * w], so we also need to multiply in the weights here (because the right-hand side was also weighted).
        for (int k = 0; k < K; k++) {
          const int idx = row + k * M;
          if (k < Kfine) {
            sten.first.add(a_fineVofs[k], a_fineWeights[k] * linAplus[idx]);
          }
          else {
            sten.second.add(a_coarVofs[k - Kfine], a_coarWeights[k - Kfine] * linAplus[idx]);
          }
        }
      }

#warning "Should make sure that A * inv(A) == 1"
    }
    else {
      MayDay::Warning("LeastSquares::computeDualLevelStencils - could not perform singular value decomposition");
    }
  }

  return ret;
}

#include <CD_NamespaceFooter.H>

#endif
