/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// A standalone spike for the cut-cell refinement operator, to measure what it costs.
//
// A cut cell is carried as twelve edge crossings and eight corner values -- twenty reals. The
// body is rebuilt from those on demand rather than stored, because with the stitch convention
// fixed they determine it completely, and storing the polygons instead would cost roughly a
// hundred reals per cell. Refinement clips the rebuilt body by the three midplanes.
//
// No Chombo headers, so the numbers below are the algorithm's own cost and nothing else.
// Reads the CSV that Prototypes/CutCellRefinement/main.cpp writes.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

constexpr int numEdges   = 12;
constexpr int numCorners = 8;
// Measured maxima over every geometry tested are ten vertices and thirteen polygons; these
// leave headroom and every write to them is guarded, because overflowing them silently drops
// faces and leaves a body that is not closed. Keeping them tight matters: the body is copied
// and scanned constantly, so its size is what the cost is made of.
constexpr int    maxVerts = 14;
constexpr int    maxPolys = 20;
constexpr double planeTol = 1.0e-12;
constexpr double edgeEps  = 1.0e-12;
constexpr double nullArea = 1.0e-18;

// Transverse directions of each coordinate direction.
constexpr int trans[3][2] = {{1, 2}, {0, 2}, {0, 1}};

long   gReason[4]    = {0, 0, 0, 0};
double gWorstClosure = 0.0;
long   gOverflow     = 0;

struct Vec
{
  double v[3];

  double&
  operator[](const int a_i)
  {
    return v[a_i];
  }

  double
  operator[](const int a_i) const
  {
    return v[a_i];
  }
};

inline Vec
operator-(const Vec& a, const Vec& b)
{
  return Vec{{a[0] - b[0], a[1] - b[1], a[2] - b[2]}};
}

inline Vec
operator+(const Vec& a, const Vec& b)
{
  return Vec{{a[0] + b[0], a[1] + b[1], a[2] + b[2]}};
}

inline Vec
operator*(const double s, const Vec& a)
{
  return Vec{{s * a[0], s * a[1], s * a[2]}};
}

inline Vec
cross(const Vec& a, const Vec& b)
{
  return Vec{{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}};
}

inline double
dot(const Vec& a, const Vec& b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline double
norm(const Vec& a)
{
  return std::sqrt(dot(a, a));
}

// A boundary polygon. face is 2*dir+side when the polygon lies in a cell face, and -1 when it
// is part of the interface.
struct Poly
{
  Vec         v[maxVerts];
  signed char vkey[maxVerts]; // edge this vertex came from, or -1 for a cell corner
  int         n    = 0;
  int         face = -1;
};

struct Body
{
  Poly p[maxPolys];
  int  n = 0;
};

// The seven moment families an EBISBox exposes, plus the true interface area, which is the
// weight the boundary centroid has to coarsen under.
struct Moments
{
  double volFrac = 0.0;
  Vec    volCentroid{{0.0, 0.0, 0.0}};
  double areaFrac[6]     = {0.0};
  Vec    faceCentroid[6] = {};
  double bndryArea       = 0.0;
  double bndryAreaTrue   = 0.0;
  Vec    normal{{0.0, 0.0, 0.0}};
  Vec    bndryCentroid{{0.0, 0.0, 0.0}};
};

inline bool
isFluid(const double a_v)
{
  // One predicate for corners and crossings alike, matching the sign test in
  // GeometryShop::insideOutsideFromNodes: the side is read off the sign bit, so a node at
  // negative zero is fluid and one at positive zero is solid.
  return std::copysign(1.0, a_v) < 0.0;
}

inline Vec
cornerPos(const int a_c)
{
  return Vec{{-0.5 + ((a_c >> 0) & 1), -0.5 + ((a_c >> 1) & 1), -0.5 + ((a_c >> 2) & 1)}};
}

// Edge e runs along direction e/4; the low two bits of e%4 give its offsets on the transverse
// directions, in the order trans[dir] lists them.
inline void
edgeGeom(const int a_e, int& a_dir, int a_off[3])
{
  a_dir    = a_e / 4;
  a_off[0] = a_off[1] = a_off[2] = 0;
  a_off[trans[a_dir][0]]         = (a_e >> 0) & 1;
  a_off[trans[a_dir][1]]         = (a_e >> 1) & 1;
}

inline int
edgeIndex(const int a_dir, const int a_off[3])
{
  return 4 * a_dir + (a_off[trans[a_dir][0]] | (a_off[trans[a_dir][1]] << 1));
}

inline void
edgeCorners(const int a_e, int& a_lo, int& a_hi)
{
  int dir, off[3];
  edgeGeom(a_e, dir, off);

  a_lo = 0;
  for (int k = 0; k < 3; k++) {
    if (k != dir) {
      a_lo |= (off[k] & 1) << k;
    }
  }

  a_hi = a_lo | (1 << dir);
}

// The four corners of face (dir,side), walked as a circuit of the square.
inline void
faceCorners(const int a_dir, const int a_side, int a_c[4])
{
  const int t0       = trans[a_dir][0];
  const int t1       = trans[a_dir][1];
  const int ab[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};

  for (int i = 0; i < 4; i++) {
    a_c[i] = (a_side << a_dir) | (ab[i][0] << t0) | (ab[i][1] << t1);
  }
}

// The four edges bounding face (dir,side), in the same circuit order: edge i joins corner i to
// corner i+1.
inline void
faceEdges(const int a_dir, const int a_side, int a_e[4])
{
  const int t0         = trans[a_dir][0];
  const int t1         = trans[a_dir][1];
  const int spec[4][3] = {{0, 0, t0}, {1, 0, t1}, {1, 1, t0}, {0, 1, t1}};

  for (int i = 0; i < 4; i++) {
    int off[3] = {0, 0, 0};
    off[a_dir] = a_side;
    off[t0]    = spec[i][0];
    off[t1]    = spec[i][1];

    const int run = spec[i][2];
    if (run == t0) {
      off[t0] = 0;
    }
    else {
      off[t1] = 0;
    }

    a_e[i] = edgeIndex(run, off);
  }
}

inline Vec
cutPoint(const int a_e, const double a_t)
{
  int dir, off[3];
  edgeGeom(a_e, dir, off);

  const double t = std::min(std::max(a_t, edgeEps), 1.0 - edgeEps);

  Vec p{{-0.5 + off[0], -0.5 + off[1], -0.5 + off[2]}};
  p[dir] = -0.5 + t;

  return p;
}

// The area vector alone. Most callers want only this or its magnitude, and computing the
// centroid as well costs a second pass over the vertices and a square root neither needs.
inline Vec
polyVec(const Poly& a_p)
{
  Vec v{{0.0, 0.0, 0.0}};

  for (int i = 1; i < a_p.n - 1; i++) {
    v = v + (0.5 * cross(a_p.v[i] - a_p.v[0], a_p.v[i + 1] - a_p.v[0]));
  }

  return v;
}

// Area, area vector and centroid of a simple planar polygon. Areas are signed about the
// polygon's own normal: a face contour that is a polyline rather than one chord leaves a
// non-convex polygon, and summing unsigned fan areas over-counts those.
inline void
polyMoments(const Poly& a_p, double& a_area, Vec& a_vec, Vec& a_centroid)
{
  a_area     = 0.0;
  a_vec      = Vec{{0.0, 0.0, 0.0}};
  a_centroid = Vec{{0.0, 0.0, 0.0}};

  if (a_p.n < 3) {
    return;
  }

  for (int i = 1; i < a_p.n - 1; i++) {
    a_vec = a_vec + (0.5 * cross(a_p.v[i] - a_p.v[0], a_p.v[i + 1] - a_p.v[0]));
  }

  a_area = norm(a_vec);
  if (a_area <= 0.0) {
    return;
  }

  const Vec u = (1.0 / a_area) * a_vec;

  Vec c{{0.0, 0.0, 0.0}};
  for (int i = 1; i < a_p.n - 1; i++) {
    const double s = 0.5 * dot(cross(a_p.v[i] - a_p.v[0], a_p.v[i + 1] - a_p.v[0]), u);
    c              = c + (s / 3.0) * (a_p.v[0] + a_p.v[i] + a_p.v[i + 1]);
  }

  a_centroid = (1.0 / a_area) * c;
}

// The chords on a face, each an ordered pair of edge indices. A face whose corners alternate
// fluid and solid has four crossings and two ways to join them; the bilinear interpolant's
// saddle value shares its sign with the pair that is connected, and it needs nothing but data
// both cells adjoining the face already hold, so they cannot disagree.
inline int
facePairs(const int    a_dir,
          const int    a_side,
          const bool   a_present[numEdges],
          const double a_phi[numCorners],
          int          a_pair[2][2])
{
  int ce[4], cc[4];
  faceEdges(a_dir, a_side, ce);
  faceCorners(a_dir, a_side, cc);

  int hit[4];
  int nhit = 0;
  for (int i = 0; i < 4; i++) {
    if (a_present[ce[i]]) {
      hit[nhit++] = i;
    }
  }

  if (nhit == 0) {
    return 0;
  }

  if (nhit == 2) {
    a_pair[0][0] = ce[hit[0]];
    a_pair[0][1] = ce[hit[1]];
    return 1;
  }

  if (nhit != 4) {
    return -1;
  }

  const double f00 = a_phi[cc[0]];
  const double f10 = a_phi[cc[1]];
  const double f11 = a_phi[cc[2]];
  const double f01 = a_phi[cc[3]];

  const double den    = f00 + f11 - f10 - f01;
  const double saddle = (std::abs(den) > 0.0) ? (f00 * f11 - f10 * f01) / den : (f00 + f11);

  // The corners alternate, so the saddle shares its side with exactly one of the two diagonals.
  // That diagonal meets through the middle; the other is the pair the two chords cut off,
  // whether it is the fluid pair or the solid one.
  if (isFluid(saddle) == isFluid(f00)) {
    a_pair[0][0] = ce[0]; // chords cut off face corners 1 and 3
    a_pair[0][1] = ce[1];
    a_pair[1][0] = ce[2];
    a_pair[1][1] = ce[3];
  }
  else {
    a_pair[0][0] = ce[3]; // chords cut off face corners 0 and 2
    a_pair[0][1] = ce[0];
    a_pair[1][0] = ce[1];
    a_pair[1][1] = ce[2];
  }

  return 2;
}

inline void
orientOutward(Poly& a_p, const int a_dir, const int a_side)
{
  Vec nout{{0.0, 0.0, 0.0}};
  nout[a_dir] = (a_side == 0) ? -1.0 : 1.0;

  Vec cr{{0.0, 0.0, 0.0}};
  for (int i = 1; i < a_p.n - 1; i++) {
    cr = cr + cross(a_p.v[i] - a_p.v[0], a_p.v[i + 1] - a_p.v[0]);
  }

  if (dot(cr, nout) < 0.0) {
    for (int i = 0; i < a_p.n / 2; i++) {
      std::swap(a_p.v[i], a_p.v[a_p.n - 1 - i]);
      std::swap(a_p.vkey[i], a_p.vkey[a_p.n - 1 - i]);
    }
  }
}

// The fluid parts of one face, as outward-oriented polygons. Returns the number written.
inline int
faceWalk(const int    a_dir,
         const int    a_side,
         const double a_cross[numEdges],
         const bool   a_present[numEdges],
         const double a_phi[numCorners],
         Poly         a_out[2])
{
  int ce[4], cc[4];
  faceEdges(a_dir, a_side, ce);
  faceCorners(a_dir, a_side, cc);

  int nhit = 0;
  for (int i = 0; i < 4; i++) {
    nhit += a_present[ce[i]] ? 1 : 0;
  }

  int       pair[2][2] = {{-1, -1}, {-1, -1}};
  const int npair      = facePairs(a_dir, a_side, a_present, a_phi, pair);
  if (npair < 0) {
    return -1;
  }

  if (nhit != 4) {
    Poly p;
    for (int i = 0; i < 4; i++) {
      if (isFluid(a_phi[cc[i]])) {
        p.vkey[p.n] = -1;
        p.v[p.n++]  = cornerPos(cc[i]);
      }
      if (a_present[ce[i]]) {
        p.vkey[p.n] = static_cast<signed char>(ce[i]);
        p.v[p.n++]  = cutPoint(ce[i], a_cross[ce[i]]);
      }
    }

    if (p.n < 3) {
      return 0;
    }

    p.face = 2 * a_dir + a_side;
    orientOutward(p, a_dir, a_side);

    const Vec pv = polyVec(p);
    if (dot(pv, pv) <= 0.0) {
      return 0;
    }

    a_out[0] = p;

    return 1;
  }

  // Four crossings: the chords cut off one diagonal pair of corners, and which pair that is is
  // the decider's whole output. Reading the fluid region off the corner signs instead lets the
  // face disagree with the loop assembly, which trusts the pairing.
  bool isolated[4] = {false, false, false, false};
  for (int k = 0; k < 2; k++) {
    for (int i = 0; i < 4; i++) {
      const int prev = ce[(i + 3) % 4];
      if ((pair[k][0] == prev && pair[k][1] == ce[i]) || (pair[k][0] == ce[i] && pair[k][1] == prev)) {
        isolated[i] = true;
      }
    }
  }

  int first = -1;
  int niso  = 0;
  for (int i = 0; i < 4; i++) {
    if (isolated[i]) {
      niso++;
      if (first < 0) {
        first = i;
      }
    }
  }

  if (niso != 2) {
    return -1;
  }

  // The cut-off corners are solid, so the fluid is one polygon wrapping around both of them.
  if (!isFluid(a_phi[cc[first]])) {
    Poly p;
    for (int i = 0; i < 4; i++) {
      if (isFluid(a_phi[cc[i]])) {
        p.vkey[p.n] = -1;
        p.v[p.n++]  = cornerPos(cc[i]);
      }
      p.vkey[p.n] = static_cast<signed char>(ce[i]);
      p.v[p.n++]  = cutPoint(ce[i], a_cross[ce[i]]);
    }

    p.face = 2 * a_dir + a_side;
    orientOutward(p, a_dir, a_side);

    const Vec pv = polyVec(p);
    if (dot(pv, pv) <= 0.0) {
      return 0;
    }

    a_out[0] = p;

    return 1;
  }

  // The cut-off corners are the fluid ones, so the fluid is two corner triangles.
  int nout = 0;
  for (int i = 0; i < 4; i++) {
    if (!isolated[i]) {
      continue;
    }

    const int prev = ce[(i + 3) % 4];

    Poly p;
    p.vkey[p.n] = static_cast<signed char>(prev);
    p.v[p.n++]  = cutPoint(prev, a_cross[prev]);
    p.vkey[p.n] = -1;
    p.v[p.n++]  = cornerPos(cc[i]);
    p.vkey[p.n] = static_cast<signed char>(ce[i]);
    p.v[p.n++]  = cutPoint(ce[i], a_cross[ce[i]]);
    p.face      = 2 * a_dir + a_side;
    orientOutward(p, a_dir, a_side);

    const Vec pv = polyVec(p);
    if (dot(pv, pv) > 0.0) {
      a_out[nout++] = p;
    }
  }

  return nout;
}

// Connected components of the fluid corners, along cell edges.
inline int
fluidComponents(const double a_phi[numCorners])
{
  bool seen[numCorners] = {false};
  int  n                = 0;

  for (int c = 0; c < numCorners; c++) {
    if (seen[c] || !isFluid(a_phi[c])) {
      continue;
    }

    n++;

    int stack[numCorners];
    int top      = 0;
    stack[top++] = c;

    while (top > 0) {
      const int x = stack[--top];
      if (seen[x]) {
        continue;
      }

      seen[x] = true;
      for (int d = 0; d < 3; d++) {
        const int y = x ^ (1 << d);
        if (isFluid(a_phi[y]) && !seen[y]) {
          stack[top++] = y;
        }
      }
    }
  }

  return n;
}

// The crossings ordered into closed loops. Each crossing lies on exactly two faces and each
// face pairs its crossings up, so every node has degree two and the graph is a union of cycles.
// Several loops just mean the interface enters as several sheets, which is legal; what must be
// refused is disconnected fluid, since that is a cell with more than one VoF.
inline int
crossingLoops(const bool   a_present[numEdges],
              const double a_phi[numCorners],
              int          a_loop[numEdges],
              int          a_start[numEdges + 1])
{
  int adj[numEdges][2];
  int deg[numEdges] = {0};

  for (int d = 0; d < 3; d++) {
    for (int s = 0; s < 2; s++) {
      int       pair[2][2];
      const int npair = facePairs(d, s, a_present, a_phi, pair);
      if (npair < 0) {
        return -1;
      }

      for (int k = 0; k < npair; k++) {
        const int x = pair[k][0];
        const int y = pair[k][1];
        if (deg[x] > 1 || deg[y] > 1) {
          return -1;
        }

        adj[x][deg[x]++] = y;
        adj[y][deg[y]++] = x;
      }
    }
  }

  int nx = 0;
  for (int e = 0; e < numEdges; e++) {
    if (a_present[e]) {
      nx++;
      if (deg[e] != 2) {
        return -1;
      }
    }
  }

  if (nx < 3) {
    return -1;
  }

  bool used[numEdges] = {false};
  int  nloop          = 0;
  int  put            = 0;
  a_start[0]          = 0;

  for (int e = 0; e < numEdges; e++) {
    if (!a_present[e] || used[e]) {
      continue;
    }

    const int begin = put;
    int       prev  = -1;
    int       cur   = e;

    while (true) {
      used[cur]     = true;
      a_loop[put++] = cur;

      const int nxt = (adj[cur][0] != prev) ? adj[cur][0] : adj[cur][1];
      if (nxt == e) {
        break;
      }
      if (used[nxt] || put > numEdges) {
        return -1;
      }

      prev = cur;
      cur  = nxt;
    }

    if (put - begin < 3) {
      return -1;
    }

    a_start[++nloop] = put;
  }

  return nloop;
}

enum class Kind
{
  Cut,
  Regular,
  Covered
};

// Whether the cell is cut at all. A cut cell needs corners strictly on both sides. A facet
// lying exactly in a cell-face plane leaves every corner at or below zero with one whole face
// at exactly zero: the cell is entirely fluid and that face is covered, because the cell beyond
// it is solid. Reading those zeros as solid instead invents crossings that bound no fluid.
inline Kind
classify(const double a_phi[numCorners])
{
  double lo = a_phi[0];
  double hi = a_phi[0];

  for (int c = 1; c < numCorners; c++) {
    lo = std::min(lo, a_phi[c]);
    hi = std::max(hi, a_phi[c]);
  }

  if (hi <= 0.0) {
    return Kind::Regular;
  }
  if (lo >= 0.0) {
    return Kind::Covered;
  }

  return Kind::Cut;
}

// Moments of a cell the interface only touches, never enters. a_B and the normal come from the
// same divergence identity as everywhere else, and the boundary is the covered face itself.
inline Moments
degenerateMoments(const double a_phi[numCorners], const Kind a_kind)
{
  Moments m;
  m.volFrac = (a_kind == Kind::Regular) ? 1.0 : 0.0;

  for (int d = 0; d < 3; d++) {
    for (int s = 0; s < 2; s++) {
      if (a_kind == Kind::Covered) {
        m.areaFrac[2 * d + s] = 0.0;
        continue;
      }

      int cc[4];
      faceCorners(d, s, cc);

      bool allZero = true;
      for (int k = 0; k < 4; k++) {
        allZero = allZero && (a_phi[cc[k]] == 0.0);
      }

      // the face lies in the interface exactly when all four of its corners do
      m.areaFrac[2 * d + s] = allZero ? 0.0 : 1.0;
    }
  }

  Vec v{{0.0, 0.0, 0.0}};
  for (int d = 0; d < 3; d++) {
    v[d] = m.areaFrac[2 * d + 1] - m.areaFrac[2 * d];
  }

  const double nrm = norm(v);
  if (nrm > 0.0) {
    m.bndryArea     = nrm;
    m.bndryAreaTrue = nrm;
    m.normal        = (1.0 / nrm) * v;

    if (a_kind == Kind::Regular) {
      for (int d = 0; d < 3; d++) {
        for (int s = 0; s < 2; s++) {
          if (m.areaFrac[2 * d + s] == 0.0) {
            m.bndryCentroid    = Vec{{0.0, 0.0, 0.0}};
            m.bndryCentroid[d] = -0.5 + s;
          }
        }
      }
    }
  }

  return m;
}

// Build the closed body bounding the fluid. Returns false when the cell is refused.
bool
buildBody(const double a_cross[numEdges], const bool a_present[numEdges], const double a_phi[numCorners], Body& a_body)
{
  a_body.n = 0;

  // Fluid arriving in several pieces is not refused: the finest level is single-valued by
  // construction, so the cell is one VoF holding their sum, exactly as GeometryShop emits one
  // IrregNode per irregular cell without ever testing connectivity. The sum is still a
  // partition of the parent, so every coarsening relation is untouched.
  gReason[0] += (fluidComponents(a_phi) > 1) ? 1 : 0;

  for (int d = 0; d < 3; d++) {
    for (int s = 0; s < 2; s++) {
      if (a_body.n + 2 > maxPolys) {
        gOverflow++;
        return false;
      }

      const int nf = faceWalk(d, s, a_cross, a_present, a_phi, &a_body.p[a_body.n]);
      if (nf < 0) {
        gReason[1]++;
        return false;
      }

      a_body.n += nf;
    }
  }

  // the face polygons occupy the front of the body; the loops are oriented against them, so
  // there is no need to walk the faces a second time for every loop
  const int nFace = a_body.n;

  int loop[numEdges];
  int start[numEdges + 1];

  const int nloop = crossingLoops(a_present, a_phi, loop, start);
  if (nloop <= 0) {
    gReason[2]++;
    return false;
  }

  // Orient each loop against the face polygons that share its chords: a shared edge is
  // traversed once each way on a closed surface.
  for (int L = 0; L < nloop; L++) {
    const int b = start[L];
    const int e = start[L + 1];

    bool done = false;
    for (int fp = 0; fp < nFace && !done; fp++) {
      const Poly& face = a_body.p[fp];

      for (int i = 0; i < face.n && !done; i++) {
        const int x = face.vkey[i];
        const int y = face.vkey[(i + 1) % face.n];
        if (x < 0 || y < 0) {
          continue;
        }

        int ix = -1, iy = -1;
        for (int j = b; j < e; j++) {
          if (loop[j] == x) {
            ix = j;
          }
          if (loop[j] == y) {
            iy = j;
          }
        }
        if (ix < 0 || iy < 0) {
          continue;
        }

        // the face polygon runs x -> y, so the patch must run y -> x
        const int nxt = b + ((iy - b + 1) % (e - b));
        if (loop[nxt] != x) {
          for (int q = 0; q < (e - b) / 2; q++) {
            std::swap(loop[b + q], loop[e - 1 - q]);
          }
        }

        done = true;
      }
    }

    // fan the loop from its centroid
    Vec ctr{{0.0, 0.0, 0.0}};
    for (int i = b; i < e; i++) {
      ctr = ctr + cutPoint(loop[i], a_cross[loop[i]]);
    }
    ctr = (1.0 / (e - b)) * ctr;

    for (int i = b; i < e; i++) {
      Poly t;
      t.v[0] = ctr;
      t.v[1] = cutPoint(loop[i], a_cross[loop[i]]);
      t.v[2] = cutPoint(loop[b + ((i - b + 1) % (e - b))], a_cross[loop[b + ((i - b + 1) % (e - b))]]);
      t.n    = 3;
      t.face = -1;

      const Vec tv = polyVec(t);
      if (dot(tv, tv) >= nullArea * nullArea) {
        if (a_body.n >= maxPolys) {
          gOverflow++;
          return false;
        }
        a_body.p[a_body.n++] = t;
      }
    }
  }

  // verify rather than assume: a folded patch leaves the body open or the volume negative, and
  // this is what keeps such a cell from being emitted as a good one
  Vec    closure{{0.0, 0.0, 0.0}};
  double vol = 0.0;
  for (int i = 0; i < a_body.n; i++) {
    closure = closure + polyVec(a_body.p[i]);

    for (int k = 1; k < a_body.p[i].n - 1; k++) {
      vol += dot(a_body.p[i].v[0], cross(a_body.p[i].v[k], a_body.p[i].v[k + 1])) / 6.0;
    }
  }

  const bool ok = norm(closure) <= 1.0e-9 && vol >= -1.0e-12 && vol <= 1.0 + 1.0e-12;
  if (!ok) {
    gReason[3]++;
    gWorstClosure = std::max(gWorstClosure, norm(closure));
  }
  return ok;
}

// Quantised vertex identity, for matching the cap segments up. Compared componentwise rather
// than hashed: a hash collision would splice two cap loops together.
struct Key
{
  long q[3];

  bool
  operator==(const Key& a_o) const
  {
    return q[0] == a_o.q[0] && q[1] == a_o.q[1] && q[2] == a_o.q[2];
  }

  bool
  operator!=(const Key& a_o) const
  {
    return !(*this == a_o);
  }
};

// Order two points so that an edge is always interpolated from the same end, whichever
// polygon is walking it. Two polygons sharing an edge traverse it in opposite directions, and
// a + t(b-a) is only algebraically equal to b + t'(a-b), not bitwise equal: the two differ in
// the last bits, and that is enough to put them in different quantisation buckets and break
// the cap loop that has to join them.
inline bool
lexLess(const Vec& a, const Vec& b)
{
  if (a[0] != b[0]) {
    return a[0] < b[0];
  }
  if (a[1] != b[1]) {
    return a[1] < b[1];
  }

  return a[2] < b[2];
}

inline Key
key(const Vec& a_p)
{
  return Key{{std::lround(a_p[0] * 1.0e11), std::lround(a_p[1] * 1.0e11), std::lround(a_p[2] * 1.0e11)}};
}

// Clip a closed body by an axis-aligned half-space, capping the opening it leaves.
void
clip(const Body& a_in, const int a_dim, const double a_coord, const bool a_keepLo, Body& a_out)
{
  a_out.n = 0;

  Vec          segA[maxPolys];
  Vec          segB[maxPolys];
  Key          keyA[maxPolys];
  Key          keyB[maxPolys];
  int          nseg = 0;
  const double sgn  = a_keepLo ? 1.0 : -1.0;

  for (int ip = 0; ip < a_in.n; ip++) {
    const Poly& p = a_in.p[ip];

    if (a_out.n >= maxPolys) {
      gOverflow++;
      return;
    }

    Poly& q = a_out.p[a_out.n];
    q.n     = 0;
    q.face  = p.face;

    for (int i = 0; i < p.n; i++) {
      const Vec&   a  = p.v[i];
      const Vec&   b  = p.v[(i + 1) % p.n];
      const double fa = sgn * (a[a_dim] - a_coord);
      const double fb = sgn * (b[a_dim] - a_coord);

      if (q.n + 2 > maxVerts) {
        gOverflow++;
        return;
      }

      if (fa <= planeTol) {
        q.vkey[q.n] = p.vkey[i];
        Vec keep    = a;
        if (fa >= -planeTol) {
          // place it exactly on the plane, or the cap it belongs to stops being recognised as
          // lying in a cell face and its area is charged to the interface
          keep[a_dim] = a_coord;
        }
        q.v[q.n++] = keep;
      }

      if ((fa < -planeTol && fb > planeTol) || (fb < -planeTol && fa > planeTol)) {
        const bool   ab = lexLess(a, b);
        const Vec&   p0 = ab ? a : b;
        const Vec&   p1 = ab ? b : a;
        const double f0 = ab ? fa : fb;
        const double f1 = ab ? fb : fa;

        const double t = f0 / (f0 - f1);
        Vec          x = p0 + t * (p1 - p0);
        x[a_dim]       = a_coord;
        q.vkey[q.n]    = -1;
        q.v[q.n++]     = x;
      }
    }

    if (q.n < 3) {
      continue;
    }

    const Vec qv = polyVec(q);
    if (dot(qv, qv) < nullArea * nullArea) {
      continue;
    }

    a_out.n++;

    bool allOn = true;
    for (int i = 0; i < q.n; i++) {
      allOn = allOn && (std::abs(q.v[i][a_dim] - a_coord) < planeTol);
    }
    if (allOn) {
      continue;
    }

    for (int i = 0; i < q.n; i++) {
      const Vec& a = q.v[i];
      const Vec& b = q.v[(i + 1) % q.n];

      if (nseg >= maxPolys) {
        gOverflow++;
        return;
      }
      if (std::abs(a[a_dim] - a_coord) < planeTol && std::abs(b[a_dim] - a_coord) < planeTol && key(a) != key(b)) {
        segA[nseg] = b;
        segB[nseg] = a;
        keyA[nseg] = key(b);
        keyB[nseg] = key(a);
        nseg++;
      }
    }
  }

  // Assemble the cap segments into closed loops, consuming SEGMENTS rather than vertices: a
  // cross-section can have several components and two of them may share a vertex, and marking
  // that vertex used drops one loop and splices the other.
  bool used[maxPolys] = {false};

  for (int s0 = 0; s0 < nseg; s0++) {
    if (used[s0]) {
      continue;
    }

    used[s0] = true;

    if (a_out.n >= maxPolys) {
      gOverflow++;
      return;
    }

    Poly& loop        = a_out.p[a_out.n];
    loop.n            = 0;
    loop.face         = 2 * a_dim + (a_keepLo ? 1 : 0);
    loop.vkey[loop.n] = -1;
    loop.v[loop.n++]  = segA[s0];

    Vec       cur    = segB[s0];
    Key       curKey = keyB[s0];
    const Key endKey = keyA[s0];
    bool      ok     = false;

    for (int guard = 0; guard <= nseg + 1; guard++) {
      if (curKey == endKey) {
        ok = true;
        break;
      }

      int nxt = -1;
      for (int j = 0; j < nseg; j++) {
        if (!used[j] && keyA[j] == curKey) {
          nxt = j;
          break;
        }
      }

      if (nxt < 0 || loop.n >= maxVerts) {
        break;
      }

      used[nxt]         = true;
      loop.vkey[loop.n] = -1;
      loop.v[loop.n++]  = cur;
      cur               = segB[nxt];
      curKey            = keyB[nxt];
    }

    if (!ok || loop.n < 3) {
      continue;
    }

    const Vec lv = polyVec(loop);
    if (dot(lv, lv) >= nullArea * nullArea) {
      a_out.n++;
    }
  }
}

// The eight children, each in its own frame.
void
refine(const Body& a_in, Body a_kid[8])
{
  // Clip as a tree. Splitting on x gives two bodies, not eight; splitting those on y gives
  // four. Doing it per octant would repeat the first two planes six extra times each.
  static thread_local Body half[2];
  static thread_local Body quarter[4];

  for (int i = 0; i < 2; i++) {
    clip(a_in, 0, 0.0, i == 0, half[i]);
  }

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      clip(half[i], 1, 0.0, j == 0, quarter[2 * i + j]);
    }
  }

  for (int o = 0; o < 8; o++) {
    clip(quarter[2 * (o & 1) + ((o >> 1) & 1)], 2, 0.0, ((o >> 2) & 1) == 0, a_kid[o]);

    // into the child's own [-0.5,0.5]^3
    Body& cur = a_kid[o];
    for (int i = 0; i < cur.n; i++) {
      for (int k = 0; k < cur.p[i].n; k++) {
        for (int d = 0; d < 3; d++) {
          cur.p[i].v[k][d] = 2.0 * (cur.p[i].v[k][d] - (-0.25 + 0.5 * ((o >> d) & 1)));
        }
      }
    }
  }
}

Moments
moments(const Body& a_b)
{
  Moments m;

  double accArea[6] = {0.0};
  Vec    accMom[6]  = {};
  Vec    ebVec{{0.0, 0.0, 0.0}};
  Vec    ebCen{{0.0, 0.0, 0.0}};
  double ebArea = 0.0;
  Vec    volMom{{0.0, 0.0, 0.0}};
  double vol = 0.0;

  for (int i = 0; i < a_b.n; i++) {
    const Poly& p = a_b.p[i];

    double area;
    Vec    vec, cen;
    polyMoments(p, area, vec, cen);

    if (area > 0.0) {
      int onFace = p.face;
      if (onFace < 0) {
        for (int d = 0; d < 3 && onFace < 0; d++) {
          for (int s = 0; s < 2 && onFace < 0; s++) {
            bool all = true;
            for (int k = 0; k < p.n; k++) {
              all = all && (std::abs(p.v[k][d] - (-0.5 + s)) < planeTol);
            }
            if (all) {
              onFace = 2 * d + s;
            }
          }
        }
      }

      if (onFace >= 0) {
        // signed: where the interface lies in a cell face it bounds a hole in that face, and
        // the aperture is the net area open to flux
        const int    d  = onFace / 2;
        const double sa = vec[d] * ((onFace % 2 == 0) ? -1.0 : 1.0);

        accArea[onFace] += sa;
        accMom[onFace] = accMom[onFace] + sa * cen;
      }
      else {
        ebArea += area;
        ebVec = ebVec + vec;
        ebCen = ebCen + area * cen;
      }
    }

    for (int k = 1; k < p.n - 1; k++) {
      const double vt = dot(p.v[0], cross(p.v[k], p.v[k + 1])) / 6.0;
      vol += vt;
      volMom = volMom + (vt * 0.25) * (p.v[0] + p.v[k] + p.v[k + 1]);
    }
  }

  m.volFrac = vol;
  if (std::abs(vol) > 0.0) {
    m.volCentroid = (1.0 / vol) * volMom;
  }

  for (int f = 0; f < 6; f++) {
    if (std::abs(accArea[f]) <= 1.0e-15) {
      m.areaFrac[f] = 0.0;
    }
    else {
      m.areaFrac[f]            = accArea[f];
      m.faceCentroid[f]        = (1.0 / accArea[f]) * accMom[f];
      m.faceCentroid[f][f / 2] = 0.0;
    }
  }

  const double nrm = norm(ebVec);
  if (nrm > 1.0e-15) {
    m.bndryArea     = nrm;
    m.bndryAreaTrue = ebArea;
    m.normal        = (-1.0 / nrm) * ebVec;
    m.bndryCentroid = (1.0 / ebArea) * ebCen;
  }

  return m;
}

// ---------------------------------------------------------------------------- driver

struct Cell
{
  double cross[numEdges];
  bool   present[numEdges];
  double phi[numCorners];
  double kappa;
  double alpha[6];
};

std::vector<Cell>
readCsv(const std::string& a_file)
{
  std::vector<Cell> out;

  FILE* fp = std::fopen(a_file.c_str(), "r");
  if (fp == nullptr) {
    return out;
  }

  char line[8192];
  if (std::fgets(line, sizeof(line), fp) == nullptr) {
    std::fclose(fp);
    return out;
  }

  std::vector<std::string> head;
  for (char* tok = std::strtok(line, ",\n"); tok != nullptr; tok = std::strtok(nullptr, ",\n")) {
    head.push_back(tok);
  }

  auto col = [&head](const std::string& a_name) {
    for (size_t i = 0; i < head.size(); i++) {
      if (head[i] == a_name) {
        return static_cast<int>(i);
      }
    }
    return -1;
  };

  const int cLevel  = col("level");
  const int cNumVoF = col("numVoFs");
  const int cKappa  = col("volFrac");

  int cEdge[numEdges];
  int cCorner[numCorners];
  for (int e = 0; e < numEdges; e++) {
    cEdge[e] = col("edgeCut" + std::to_string(e));
  }
  for (int c = 0; c < numCorners; c++) {
    cCorner[c] = col("corner" + std::to_string(c));
  }

  int cFaces[6];
  int cAlpha[6];
  for (int f = 0; f < 6; f++) {
    cFaces[f] = col("numFaces" + std::to_string(f));
    cAlpha[f] = col("areaFrac" + std::to_string(f));
  }

  int                              topLevel = 0;
  std::vector<std::vector<double>> rows;

  while (std::fgets(line, sizeof(line), fp) != nullptr) {
    std::vector<double> r;
    for (char* tok = std::strtok(line, ",\n"); tok != nullptr; tok = std::strtok(nullptr, ",\n")) {
      r.push_back(std::atof(tok));
    }
    if (r.size() < head.size()) {
      continue;
    }
    topLevel = std::max(topLevel, static_cast<int>(r[cLevel]));
    rows.push_back(std::move(r));
  }
  std::fclose(fp);

  for (const auto& r : rows) {
    if (static_cast<int>(r[cLevel]) != topLevel || static_cast<int>(r[cNumVoF]) != 1) {
      continue;
    }

    bool ok = true;
    for (int f = 0; f < 6; f++) {
      ok = ok && (static_cast<int>(r[cFaces[f]]) <= 1);
    }
    if (!ok) {
      continue;
    }

    Cell c;
    for (int e = 0; e < numEdges; e++) {
      const double t = r[cEdge[e]];
      c.present[e]   = (t >= 0.0);
      c.cross[e]     = t;
    }
    for (int k = 0; k < numCorners; k++) {
      c.phi[k] = r[cCorner[k]];
    }
    c.kappa = r[cKappa];
    for (int f = 0; f < 6; f++) {
      c.alpha[f] = (static_cast<int>(r[cFaces[f]]) == 1) ? r[cAlpha[f]] : 0.0;
    }

    // keep the same predicate everywhere: a crossing whose end corners agree bounds no fluid
    for (int e = 0; e < numEdges; e++) {
      if (c.present[e]) {
        int lo, hi;
        edgeCorners(e, lo, hi);
        c.present[e] = (isFluid(c.phi[lo]) != isFluid(c.phi[hi]));
      }
    }

    out.push_back(c);
  }

  return out;
}

int
descend(const Body& a_b, const int a_depth, const int a_maxDepth, double& a_worst)
{
  const Moments parent = moments(a_b);

  if (a_depth == a_maxDepth) {
    return 1;
  }

  Body kid[8];
  refine(a_b, kid);

  double sumK    = 0.0;
  double sumA[6] = {0.0};
  int    count   = 1;

  for (int o = 0; o < 8; o++) {
    const Moments m = moments(kid[o]);
    sumK += m.volFrac;

    for (int d = 0; d < 3; d++) {
      for (int s = 0; s < 2; s++) {
        if (((o >> d) & 1) == s) {
          sumA[2 * d + s] += m.areaFrac[2 * d + s];
        }
      }
    }

    count += descend(kid[o], a_depth + 1, a_maxDepth, a_worst);
  }

  const double kerr = std::abs(sumK / 8.0 - parent.volFrac);

  a_worst = std::max(a_worst, kerr);
  for (int f = 0; f < 6; f++) {
    a_worst = std::max(a_worst, std::abs(sumA[f] / 4.0 - parent.areaFrac[f]));
  }

  return count;
}

} // namespace

int
main(int argc, char* argv[])
{
  if (argc < 2) {
    std::printf("usage: %s <cutcells.csv> [depth]\n", argv[0]);
    return 1;
  }

  const int  depth = (argc > 2) ? std::atoi(argv[2]) : 3;
  const bool dump  = (argc > 3) && (std::string(argv[3]) == "--dump");

  const std::vector<Cell> cells = readCsv(argv[1]);
  if (cells.empty()) {
    std::printf("no cells read from %s\n", argv[1]);
    return 1;
  }

  if (dump) {
    int idx = 0;
    for (const auto& c : cells) {
      Body b;
      if (!buildBody(c.cross, c.present, c.phi, b)) {
        std::printf("%d REFUSED\n", idx);
      }
      else {
        Body kid[8];
        refine(b, kid);
        std::printf("%d CUT %.12f", idx, moments(b).volFrac);
        for (int o = 0; o < 8; o++) {
          std::printf(" %.9f", moments(kid[o]).volFrac);
        }
        std::printf("\n");
      }
      idx++;
    }
    return 0;
  }

  // build
  auto t0 = std::chrono::steady_clock::now();

  std::vector<Body> bodies;
  bodies.reserve(cells.size());
  int    refused    = 0;
  int    nRegular   = 0;
  int    nCovered   = 0;
  double worstDegen = 0.0;

  for (const auto& c : cells) {
    const Kind k = classify(c.phi);
    if (k != Kind::Cut) {
      // a cell the interface only touches refines into eight of the same, so its coarsening
      // relations hold trivially; what is worth checking is that it matches what Chombo stored
      const Moments m = degenerateMoments(c.phi, k);
      worstDegen      = std::max(worstDegen, std::abs(m.volFrac - c.kappa));
      for (int f = 0; f < 6; f++) {
        worstDegen = std::max(worstDegen, std::abs(m.areaFrac[f] - c.alpha[f]));
      }

      (k == Kind::Regular) ? nRegular++ : nCovered++;
      continue;
    }

    Body b;
    if (buildBody(c.cross, c.present, c.phi, b)) {
      bodies.push_back(b);
    }
    else {
      refused++;
    }
  }

  auto t1 = std::chrono::steady_clock::now();

  // agreement with what Chombo stored, on the cells that built
  double worstKappa = 0.0;
  for (const auto& c : cells) {
    if (classify(c.phi) != Kind::Cut) {
      continue;
    }

    Body b;
    if (!buildBody(c.cross, c.present, c.phi, b)) {
      continue;
    }

    worstKappa = std::max(worstKappa, std::abs(moments(b).volFrac - c.kappa));
  }

  // refine
  auto   t2       = std::chrono::steady_clock::now();
  double worstCon = 0.0;
  long   cellsOut = 0;

  for (const auto& b : bodies) {
    cellsOut += descend(b, 0, depth, worstCon);
  }

  auto t3 = std::chrono::steady_clock::now();

  const double buildMs  = std::chrono::duration<double, std::milli>(t1 - t0).count();
  const double refineMs = std::chrono::duration<double, std::milli>(t3 - t2).count();

  std::printf("  %-34s %7zu cut cells, %d refused\n", argv[1], cells.size(), refused);
  std::printf("    build        %9.2f ms   %8.0f ns/cell\n",
              buildMs,
              1.0e6 * buildMs / std::max<size_t>(cells.size(), 1));
  std::printf("    refine d=%d   %9.2f ms   %8.0f ns/output cell   (%ld cells)\n",
              depth,
              refineMs,
              1.0e6 * refineMs / std::max<long>(cellsOut, 1),
              cellsOut);
  std::printf("    worst |kappa - chombo| %.2e     worst conservation %.2e\n", worstKappa, worstCon);
  std::printf("    touch-only cells: %d regular, %d covered, worst gap to chombo %.2e\n",
              nRegular,
              nCovered,
              worstDegen);
  std::printf("    capacity overflows %ld\n", gOverflow);
  std::printf("    merged multi-piece cells %ld;  refusals: face-count %ld  loops %ld  invalid-body %ld"
              "  (worst closure %.2e)\n",
              gReason[0] / 2,
              gReason[1] / 2,
              gReason[2] / 2,
              gReason[3] / 2,
              gWorstClosure);

  return 0;
}
