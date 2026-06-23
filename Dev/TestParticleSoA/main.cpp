/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  Leaf-level unit test for ParticleSoA::deepCopy() and deepCopyTo().

  Payload has an FP column (phi) and an INTEGER column (tag), so the test also confirms non-FP
  payload columns copy correctly. Checks:
    - deepCopy: independent arena (different data()) with a bit-for-bit identical container;
    - mutating a copy leaves the original unchanged (no aliasing);
    - cell-sort state (isSorted + CSR offsets) is preserved by the copy;
    - deepCopyTo: writes into an existing destination, REUSES its allocation when already large
      enough (destination data() pointer unchanged), and is a no-op when copying onto self;
    - deepCopy/deepCopyTo of an empty container yields an empty container.
  Run single-rank.
*/

// Std includes
#include <cstddef>
#include <cstdint>
#include <tuple>

// Chombo includes
#include <Box.H>
#include <RealVect.H>
#include <IntVect.H>
#include <MayDay.H>
#include <parstream.H>

// Our includes
#include <CD_Initialize.H>
#include <CD_ParticleSoA.H>

using namespace ChomboDischarge;

namespace ChomboDischarge {

  /** @brief Payload with one floating-point column and one integer column. */
  struct Pay
  {
    Real         phi = 0.0;
    std::int32_t tag = 0;
  };
  template <>
  struct ParticleTraits<Pay>
  {
    static constexpr auto columns = std::make_tuple(&Pay::phi, &Pay::tag);
  };

} // namespace ChomboDischarge

namespace {

  using PSoA = ParticleSoA<Pay>;

  void
  require(const bool a_ok, const std::string& a_what)
  {
    if (!a_ok) {
      MayDay::Abort(("TestParticleSoA: check failed -- " + a_what).c_str());
    }
  }

  /** @brief Whether two containers are bit-for-bit identical (every column + sort state). */
  bool
  identical(const PSoA& a_x, const PSoA& a_y)
  {
    if (a_x.size() != a_y.size()) {
      return false;
    }
    if (a_x.isSorted() != a_y.isSorted()) {
      return false;
    }
    if (a_x.isSorted()) {
      if (a_x.numCells() != a_y.numCells()) {
        return false;
      }
      for (std::size_t c = 0; c <= a_x.numCells(); c++) {
        if (a_x.cellStart(c) != a_y.cellStart(c)) {
          return false;
        }
      }
    }
    for (std::size_t i = 0; i < a_x.size(); i++) {
      if (a_x.weight(i) != a_y.weight(i) || a_x.particleID(i) != a_y.particleID(i) || a_x.rankID(i) != a_y.rankID(i) ||
          a_x.column<&Pay::phi>()[i] != a_y.column<&Pay::phi>()[i] ||
          a_x.column<&Pay::tag>()[i] != a_y.column<&Pay::tag>()[i]) {
        return false;
      }
      for (int dir = 0; dir < SpaceDim; dir++) {
        if (a_x.position(i)[dir] != a_y.position(i)[dir]) {
          return false;
        }
      }
    }
    return true;
  }

  /** @brief Whether particle i of a_x equals particle j of a_y across every column (incl. id/rank). */
  bool
  sameParticle(const PSoA& a_x, const std::size_t a_i, const PSoA& a_y, const std::size_t a_j)
  {
    bool ok = (a_x.weight(a_i) == a_y.weight(a_j)) && (a_x.particleID(a_i) == a_y.particleID(a_j)) &&
              (a_x.rankID(a_i) == a_y.rankID(a_j)) && (a_x.column<&Pay::phi>()[a_i] == a_y.column<&Pay::phi>()[a_j]) &&
              (a_x.column<&Pay::tag>()[a_i] == a_y.column<&Pay::tag>()[a_j]);
    for (int dir = 0; dir < SpaceDim; dir++) {
      ok = ok && (a_x.position(a_i)[dir] == a_y.position(a_j)[dir]);
    }
    return ok;
  }

  /** @brief Populate a container with a_n deterministic particles (in cells [0,7]). */
  void
  fill(PSoA& a_soa, const int a_n, const int a_seed)
  {
    for (int i = 0; i < a_n; i++) {
      RealVect x;
      for (int dir = 0; dir < SpaceDim; dir++) {
        x[dir] = 0.5 + static_cast<Real>((i + a_seed) % 8);
      }
      a_soa.append(x, 1.0 + 0.25 * static_cast<double>(i + a_seed), Pay{3.0 * (i + a_seed), i + a_seed});
      a_soa.particleID(a_soa.size() - 1) = 1000 + i + a_seed;
      a_soa.rankID(a_soa.size() - 1)     = (i + a_seed) % 4;
    }
  }

} // namespace

int
main(int argc, char* argv[])
{
  initialize(argc, argv);

  {
    const RealVect probLo = RealVect::Zero;
    const Box      box(IntVect::Zero, 7 * IntVect::Unit);

    PSoA a;
    fill(a, 50, 0);
    require(a.size() == 50, "populated");

    // ---- (1) deepCopy: independent arena, identical contents ----
    PSoA b = a.deepCopy();
    require(identical(a, b), "deepCopy is bit-for-bit identical (incl. integer payload column)");
    require(a.data() != b.data(), "deepCopy owns a separate arena");

    // ---- (2) mutating the copy leaves the original untouched ----
    const RealVect     x0   = a.position(0);
    const double       w0   = a.weight(0);
    const Real         phi0 = a.column<&Pay::phi>()[0];
    const std::int32_t tag0 = a.column<&Pay::tag>()[0];
    b.setPosition(0, 99.0 * RealVect::Unit);
    b.weight(0)              = -7.0;
    b.column<&Pay::phi>()[0] = -1.0;
    b.column<&Pay::tag>()[0] = -99;
    require(!identical(a, b), "mutating the copy diverges from the original");
    bool unchanged = (a.weight(0) == w0) && (a.column<&Pay::phi>()[0] == phi0) && (a.column<&Pay::tag>()[0] == tag0);
    for (int dir = 0; dir < SpaceDim; dir++) {
      unchanged = unchanged && (a.position(0)[dir] == x0[dir]);
    }
    require(unchanged, "original is untouched by mutating the copy");

    // ---- (3) deepCopy preserves the cell-sort state ----
    a.sortByCell(box, RealVect::Unit, probLo);
    require(a.isSorted(), "source sorted");
    PSoA c = a.deepCopy();
    require(identical(a, c), "deepCopy of a sorted container matches (CSR offsets included)");
    require(a.data() != c.data(), "sorted deepCopy owns a separate arena");

    // ---- (4) deepCopyTo a fresh destination ----
    PSoA dst;
    a.deepCopyTo(dst);
    require(identical(a, dst), "deepCopyTo reproduces the source");
    require(a.data() != dst.data(), "deepCopyTo destination owns a separate arena");

    // ---- (5) deepCopyTo REUSES the destination allocation when it is already large enough ----
    void* const reused = dst.data(); // dst now has capacity >= 50
    PSoA        small;
    fill(small, 10, 100); // unsorted, smaller
    small.deepCopyTo(dst);
    require(dst.data() == reused, "deepCopyTo reuses the destination arena (no realloc when it fits)");
    require(identical(small, dst), "deepCopyTo overwrote the destination with the new source");

    // ---- (6) deepCopyTo onto self is a no-op ----
    PSoA snap = a.deepCopy();
    a.deepCopyTo(a);
    require(identical(a, snap), "deepCopyTo onto self leaves the container unchanged");

    // ---- (7) empty copies ----
    PSoA empty;
    require(empty.deepCopy().size() == 0, "deepCopy of empty is empty");
    empty.deepCopyTo(dst);
    require(dst.size() == 0, "deepCopyTo of empty empties the destination");

    // ---- (8) two-argument append (default-constructed payload) ----
    {
      PSoA p;
      p.append(RealVect::Unit, 2.0);
      require(p.size() == 1, "2-arg append adds a particle");
      require(p.column<&Pay::phi>()[0] == 0.0 && p.column<&Pay::tag>()[0] == 0, "2-arg append default payload");
    }

    // ---- (9) bulk append (catenate) preserves ALL columns incl. id/rank ----
    {
      PSoA u;
      fill(u, 7, 0);
      PSoA v;
      fill(v, 5, 100);
      PSoA uSnap = u.deepCopy();
      u.append(v);
      require(u.size() == 12, "bulk append size");
      for (std::size_t i = 0; i < 7; i++) {
        require(sameParticle(u, i, uSnap, i), "bulk append keeps the original particles (incl. id/rank)");
      }
      for (std::size_t j = 0; j < 5; j++) {
        require(sameParticle(u, 7 + j, v, j), "bulk append appends the source particles (incl. id/rank)");
      }
    }

    // ---- (9b) self bulk append ----
    {
      PSoA s;
      fill(s, 6, 3);
      PSoA sSnap = s.deepCopy();
      s.append(s);
      require(s.size() == 12, "self bulk append size");
      for (std::size_t i = 0; i < 6; i++) {
        require(sameParticle(s, i, sSnap, i) && sameParticle(s, 6 + i, sSnap, i), "self bulk append duplicates");
      }
    }

    // ---- (9c) catenate (move-append): zero-copy steal into an empty destination ----
    {
      PSoA src;
      fill(src, 8, 0);
      PSoA        srcSnap = src.deepCopy();
      void* const srcData = src.data();
      PSoA        dst; // empty
      dst.catenate(src);
      require(dst.data() == srcData, "catenate into empty steals the source arena (zero copy)");
      require(src.size() == 0, "catenate empties the source");
      require(identical(dst, srcSnap), "catenate into empty moves all particles");
    }

    // ---- (9d) catenate into a non-empty destination ----
    {
      PSoA d;
      fill(d, 6, 0);
      PSoA s;
      fill(s, 4, 100);
      PSoA dSnap = d.deepCopy();
      PSoA sSnap = s.deepCopy();
      d.catenate(s);
      require(d.size() == 10 && s.size() == 0, "catenate appends and empties the source");
      for (std::size_t i = 0; i < 6; i++) {
        require(sameParticle(d, i, dSnap, i), "catenate keeps the destination particles (incl. id/rank)");
      }
      for (std::size_t j = 0; j < 4; j++) {
        require(sameParticle(d, 6 + j, sSnap, j), "catenate appends the source particles (incl. id/rank)");
      }
      // self-catenate is a no-op
      PSoA self;
      fill(self, 5, 7);
      PSoA selfSnap = self.deepCopy();
      self.catenate(self);
      require(identical(self, selfSnap), "catenate onto self is a no-op");
    }

    // ---- (10) swap exchanges arenas + contents in O(1) ----
    {
      PSoA p;
      fill(p, 4, 0);
      PSoA q;
      fill(q, 9, 50);
      PSoA        pSnap = p.deepCopy();
      PSoA        qSnap = q.deepCopy();
      void* const pd    = p.data();
      void* const qd    = q.data();
      p.swap(q);
      require(p.size() == 9 && q.size() == 4, "swap exchanges sizes");
      require(p.data() == qd && q.data() == pd, "swap exchanges arenas (no data moved)");
      require(identical(p, qSnap) && identical(q, pSnap), "swap exchanges contents");
    }

    // ---- (11) shrinkToFit reclaims capacity, preserving contents + sort ----
    {
      PSoA p;
      p.reserve(1000);
      fill(p, 10, 0);
      require(p.capacity() >= 1000, "reserved capacity");
      PSoA snap = p.deepCopy();
      p.shrinkToFit();
      require(p.capacity() == p.size() && p.size() == 10, "shrinkToFit compacts to size");
      require(identical(p, snap), "shrinkToFit preserves contents");

      p.sortByCell(box, RealVect::Unit, probLo); // re-grows capacity internally if needed
      PSoA snap2 = p.deepCopy();
      p.shrinkToFit();
      require(p.isSorted() && identical(p, snap2), "shrinkToFit preserves sort + contents");

      PSoA e;
      e.reserve(64);
      e.shrinkToFit();
      require(e.capacity() == 0 && e.size() == 0, "shrinkToFit on empty frees the arena");
    }

    // ---- (12) cellRange matches the CSR ----
    {
      PSoA p;
      fill(p, 30, 0);
      p.sortByCell(box, RealVect::Unit, probLo);
      for (std::size_t cell = 0; cell < p.numCells(); cell++) {
        const auto r = p.cellRange(cell);
        require(r.first == p.cellStart(cell) && r.second == p.cellStart(cell + 1), "cellRange matches cellStart");
        require(r.second - r.first == p.particlesInCell(cell), "cellRange size matches particlesInCell");
      }
    }

    // ---- (13) appendParticle (single-particle copy) + transfer (appendParticle + remove) ----
    {
      PSoA src;
      fill(src, 8, 0);
      PSoA dst;
      fill(dst, 3, 50);

      // Copy src particle 5 onto the end of dst (all columns, incl. id/rank).
      dst.appendParticle(src, 5);
      require(dst.size() == 4, "appendParticle grows the destination by one");
      require(sameParticle(dst, 3, src, 5), "appendParticle copies every column (incl. id/rank)");
      require(src.size() == 8, "appendParticle does not touch the source");

      // Transfer: appendParticle + swap-pop remove == move one particle between containers.
      PSoA a;
      fill(a, 6, 0);
      PSoA              b;
      const std::size_t k     = 2;
      PSoA              snapK = a.deepCopy();
      b.appendParticle(a, k);
      a.remove(k);
      require(a.size() == 5 && b.size() == 1, "transfer moves one particle");
      require(sameParticle(b, 0, snapK, k), "transferred particle matches the original");

      // Self-append duplicates a particle.
      PSoA s;
      fill(s, 4, 7);
      s.appendParticle(s, 1);
      require(s.size() == 5 && sameParticle(s, 4, s, 1), "self appendParticle duplicates");
    }

    pout() << "All ParticleSoA deepCopy/deepCopyTo/append/appendParticle/catenate/swap/shrinkToFit/cellRange "
              "checks passed ("
           << SpaceDim << "D)." << endl;
  }

  return finalize();
}
