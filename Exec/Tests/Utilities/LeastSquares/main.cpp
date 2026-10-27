/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   main.cpp
 * @brief  Unit test for LeastSquares: the stencils must reproduce polynomials up to their order.
 * @author Robert Marskar
 */

// Std includes
#include <functional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// Chombo includes
#include <BoxIterator.H>
#include <ParmParse.H>

// Our includes
#include <CD_Driver.H>
#include <CD_LeastSquares.H>
#include <CD_MultiIndex.H>

using namespace ChomboDischarge;

/**
 * @brief Polynomial sum_a c_a * x^a / a! up to a given order, so that its derivative a at the origin is c_a.
 */
struct Polynomial
{
  /**
   * @brief Terms of the polynomial: multi-index and coefficient.
   */
  std::vector<std::pair<IntVect, Real>> terms;

  /**
   * @brief Evaluate the polynomial.
   * @param[in] a_x Position
   * @return The value at a_x.
   */
  Real
  operator()(const RealVect& a_x) const
  {
    Real value = 0.0;

    for (const auto& term : terms) {
      Real monomial = term.second;

      for (int dir = 0; dir < SpaceDim; dir++) {
        monomial *= std::pow(a_x[dir], term.first[dir]);

        for (int k = 2; k <= term.first[dir]; k++) {
          monomial /= k;
        }
      }

      value += monomial;
    }

    return value;
  }

  /**
   * @brief Get the coefficient of a term, which is the derivative a_alpha at the origin.
   * @param[in] a_alpha Multi-index of the term
   * @return The coefficient, or zero if there is no such term.
   */
  Real
  coefficient(const IntVect& a_alpha) const
  {
    for (const auto& term : terms) {
      if (term.first == a_alpha) {
        return term.second;
      }
    }

    return 0.0;
  }
};

/**
 * @brief Describe a test case for the error messages.
 * @param[in] a_name  Case name
 * @param[in] a_dx    Resolution
 * @param[in] a_order Order
 * @return The description.
 */
std::string
describe(const std::string& a_name, const Real a_dx, const int a_order)
{
  std::ostringstream os;

  os << a_name << ", dx = " << a_dx << ", order " << a_order;

  return os.str();
}

/**
 * @brief Random polynomial of a given order whose terms are O(1) at distances O(a_dx).
 * @param[in]     a_order     Order
 * @param[in]     a_dx        Length scale
 * @param[in,out] a_generator Random number generator
 * @return The polynomial.
 */
Polynomial
randomPolynomial(const int a_order, const Real a_dx, std::mt19937& a_generator)
{
  std::uniform_real_distribution<Real> uniform(-1.0, 1.0);

  Polynomial polynomial;

  for (MultiIndex mi(a_order); mi.ok(); ++mi) {
    const IntVect alpha = mi.getCurrentIndex();

    polynomial.terms.emplace_back(alpha, uniform(a_generator) * std::pow(a_dx, -alpha.sum()));
  }

  return polynomial;
}

/**
 * @brief Displacements of the cells in a box of radius a_radius around the origin, without the origin, that satisfy a
 * predicate.
 * @param[in] a_radius Radius of the box
 * @param[in] a_keep   Predicate on the cell offset
 * @return The offsets.
 */
std::vector<IntVect>
neighborhood(const int a_radius, const std::function<bool(const IntVect&)>& a_keep)
{
  std::vector<IntVect> offsets;

  for (BoxIterator bit(Box(-a_radius * IntVect::Unit, a_radius * IntVect::Unit)); bit.ok(); ++bit) {
    if (bit() != IntVect::Zero && a_keep(bit())) {
      offsets.push_back(bit());
    }
  }

  return offsets;
}

/**
 * @brief Abort with a message unless a condition holds.
 * @param[in] a_condition Condition
 * @param[in] a_message   Message
 */
void
check(const bool a_condition, const std::string& a_message)
{
  if (!a_condition) {
    MayDay::Abort(("LeastSquares test -- " + a_message).c_str());
  }
}

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  ParmParse pp;

  Real tolerance;

  pp.get("tolerance", tolerance);

  std::mt19937 generator(1234);

  constexpr int numPolynomials = 20;

  const std::vector<Real> resolutions{1.0, 1.E-3, 1.E-5, 1.E-8};

  IntVectSet gradient;
  for (int dir = 0; dir < SpaceDim; dir++) {
    gradient |= BASISV(dir);
  }

  // The gradient at the origin from the cells around it, with the value at the origin known (as in EBGradient). The
  // neighborhoods are cut by planes, as next to an embedded boundary. The stencil must be exact for polynomials up to
  // the orders that the neighborhood determines. The rejected order is one that the neighborhood cannot determine,
  // because it has too few distinct coordinates in some direction, and for which no stencil may be returned.
  struct Case
  {
    std::string          name;
    std::vector<IntVect> offsets;
    int                  maxExactOrder;
    int                  rejectedOrder;
  };

  const std::vector<Case> cases{{"radius 2",
                                 neighborhood(2,
                                              [](const IntVect&) {
                                                return true;
                                              }),
                                 3,
                                 -1},
                                {"radius 2, oblique cut",
                                 neighborhood(2,
                                              [](const IntVect& iv) {
                                                return 2 * iv[0] + iv[1] >= 0;
                                              }),
                                 3,
                                 -1},
                                {"radius 2, half space",
                                 neighborhood(2,
                                              [](const IntVect& iv) {
                                                return iv[0] >= 0;
                                              }),
                                 2,
                                 3},
                                {"radius 1, half space",
                                 neighborhood(1,
                                              [](const IntVect& iv) {
                                                return iv[0] >= 0;
                                              }),
                                 1,
                                 2},
                                {"radius 2, line",
                                 neighborhood(2,
                                              [](const IntVect& iv) {
                                                for (int dir = 1; dir < SpaceDim; dir++) {
                                                  if (iv[dir] != 0) {
                                                    return false;
                                                  }
                                                }

                                                return true;
                                              }),
                                 0,
                                 1}};

  for (const auto& c : cases) {
    for (const Real dx : resolutions) {
      for (int order = 1; order <= std::max(c.maxExactOrder, c.rejectedOrder); order++) {
        const int numUnknowns = LeastSquares::getTaylorExpansionSize(order) - 1;

        if (static_cast<int>(c.offsets.size()) < numUnknowns) {
          continue;
        }

        Vector<VolIndex> vofs;
        Vector<RealVect> displacements;

        for (int i = 0; i < c.offsets.size(); i++) {
          vofs.push_back(VolIndex(static_cast<int>(i) * BASISV(0), 0));
          displacements.push_back(dx * RealVect(c.offsets[i]));
        }

        const auto stencils = LeastSquares::computeSingleLevelStencils(gradient,
                                                                       IntVectSet(IntVect::Zero),
                                                                       vofs,
                                                                       displacements,
                                                                       1,
                                                                       order);

        const std::string what = describe(c.name, dx, order);

        if (order == c.rejectedOrder) {
          for (int dir = 0; dir < SpaceDim; dir++) {
            check(stencils.at(BASISV(dir)).size() == 0, what + ": rank deficient system was not rejected");
          }

          continue;
        }

        if (order > c.maxExactOrder) {
          continue;
        }

        for (int n = 0; n < numPolynomials; n++) {
          const Polynomial polynomial = randomPolynomial(order, dx, generator);
          const Real       value0     = polynomial(RealVect::Zero);

          for (int dir = 0; dir < SpaceDim; dir++) {
            const VoFStencil& stencil = stencils.at(BASISV(dir));

            check(stencil.size() > 0, what + ": no stencil");

            Real derivative = 0.0;

            for (int i = 0; i < stencil.size(); i++) {
              derivative += stencil.weight(i) * (polynomial(displacements[stencil.vof(i).gridIndex()[0]]) - value0);
            }

            // The derivative is O(1/dx), so the error is measured in units of 1/dx.
            const Real exact = polynomial.coefficient(BASISV(dir));

            check(std::abs(derivative - exact) * dx <= tolerance * (1.0 + std::abs(exact) * dx),
                  what + ": gradient is not exact");
          }
        }
      }
    }
  }

  // Interpolation to a fine ghost cell at the origin from the fine cells next to it and the coarse cells around it, as
  // in EBLeastSquaresMultigridInterpolator. The refinement ratio is 2.
  for (const Real dx : resolutions) {
    for (int order = 1; order <= 3; order++) {
      Vector<VolIndex> fineVofs;
      Vector<VolIndex> coarVofs;
      Vector<RealVect> fineDisplacements;
      Vector<RealVect> coarDisplacements;

      for (const auto& iv : neighborhood(2, [](const IntVect& iv) {
             return iv[0] >= 1;
           })) {
        fineVofs.push_back(VolIndex(static_cast<int>(fineVofs.size()) * BASISV(0), 0));
        fineDisplacements.push_back(dx * RealVect(iv));
      }

      for (const auto& iv : neighborhood(1, [](const IntVect&) {
             return true;
           })) {
        coarVofs.push_back(VolIndex(static_cast<int>(coarVofs.size()) * BASISV(0), 0));
        coarDisplacements.push_back(dx * (2.0 * RealVect(iv) + 0.5 * RealVect::Unit));
      }

      const auto stencils = LeastSquares::computeDualLevelStencils<Real>(IntVectSet(IntVect::Zero),
                                                                         IntVectSet(),
                                                                         fineVofs,
                                                                         coarVofs,
                                                                         fineDisplacements,
                                                                         coarDisplacements,
                                                                         1,
                                                                         order);

      const VoFStencil& fineStencil = stencils.at(IntVect::Zero).first;
      const VoFStencil& coarStencil = stencils.at(IntVect::Zero).second;

      const std::string what = describe("dual level", dx, order);

      check(fineStencil.size() + coarStencil.size() > 0, what + ": no stencil");

      for (int n = 0; n < numPolynomials; n++) {
        const Polynomial polynomial = randomPolynomial(order, dx, generator);

        Real value = 0.0;

        for (int i = 0; i < fineStencil.size(); i++) {
          value += fineStencil.weight(i) * polynomial(fineDisplacements[fineStencil.vof(i).gridIndex()[0]]);
        }

        for (int i = 0; i < coarStencil.size(); i++) {
          value += coarStencil.weight(i) * polynomial(coarDisplacements[coarStencil.vof(i).gridIndex()[0]]);
        }

        const Real exact = polynomial(RealVect::Zero);

        check(std::abs(value - exact) <= tolerance * (1.0 + std::abs(exact)), what + ": interpolation is not exact");
      }
    }
  }

  pout() << "LeastSquares test passed" << endl;

  ChomboDischarge::finalize();
}
