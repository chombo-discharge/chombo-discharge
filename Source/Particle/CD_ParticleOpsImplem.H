/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ParticleOpsImplem.H
  @brief  Implementation of CD_ParticleOps.H
  @author Robert Marskar
*/

#ifndef CD_ParticleOpsImplem_H
#define CD_ParticleOpsImplem_H

// Chombo includes
#include <CH_Timer.H>
#include <PolyGeom.H>

// Our includes
#include <CD_ParticleOps.H>
#include <CD_PolyUtils.H>
#include <CD_Random.H>
#include <CD_NamespaceHeader.H>

inline IntVect
ParticleOps::getParticleCellIndex(const RealVect& a_particlePosition,
                                  const RealVect& a_probLo,
                                  const Real&     a_dx) noexcept
{
  return IntVect(D_DECL(std::floor((a_particlePosition[0] - a_probLo[0]) / a_dx),
                        std::floor((a_particlePosition[1] - a_probLo[1]) / a_dx),
                        std::floor((a_particlePosition[2] - a_probLo[2]) / a_dx)));
}

inline IntVect
ParticleOps::getParticleCellIndex(const RealVect& a_particlePosition,
                                  const RealVect& a_probLo,
                                  const RealVect& a_dx) noexcept
{
  return IntVect(D_DECL(std::floor((a_particlePosition[0] - a_probLo[0]) / a_dx[0]),
                        std::floor((a_particlePosition[1] - a_probLo[1]) / a_dx[1]),
                        std::floor((a_particlePosition[2] - a_probLo[2]) / a_dx[2])));
}

template <typename P, const Real& (P::*weight)() const>
inline void
ParticleOps::getPhysicalParticlesPerCell(EBAMRCellData& a_ppc, const ParticleContainer<P>& a_src) noexcept
{
  CH_TIME("ParticleOps::getPhysicalParticlesPerCell");

  const RealVect probLo = a_src.getProbLo();

  for (int lvl = 0; lvl <= a_src.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_src.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();
    const RealVect           dx  = a_src.getDx()[lvl];

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const List<P>& particles = a_src[lvl][din].listItems();

      FArrayBox& ppc = (*a_ppc[lvl])[din].getFArrayBox();
      ppc.setVal(0.0);

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        const P&       p  = lit();
        const RealVect x  = p.position();
        const Real     w  = (p.*weight)();
        const IntVect  iv = ParticleOps::getParticleGridCell(x, probLo, dx);

        ppc(iv, 0) += w;
      }
    }
  }
}

template <typename P>
inline void
ParticleOps::getComputationalParticlesPerCell(EBAMRCellData& a_ppc, const ParticleContainer<P>& a_src) noexcept
{
  CH_TIME("ParticleOps::getComputationalParticlesPerCell");

  const RealVect probLo = a_src.getProbLo();

  for (int lvl = 0; lvl <= a_src.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_src.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();
    const RealVect           dx  = a_src.getDx()[lvl];

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const List<P>& particles = a_src[lvl][din].listItems();

      FArrayBox& ppc = (*a_ppc[lvl])[din].getFArrayBox();
      ppc.setVal(0.0);

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        const P&       p  = lit();
        const RealVect x  = p.position();
        const IntVect  iv = ParticleOps::getParticleGridCell(x, probLo, dx);

        ppc(iv, 0) += 1.0;
      }
    }
  }
}

inline IntVect
ParticleOps::getParticleGridCell(const RealVect& a_particlePosition,
                                 const RealVect& a_probLo,
                                 const RealVect& a_dx) noexcept
{
  return IntVect(D_DECL(std::floor((a_particlePosition[0] - a_probLo[0]) / a_dx[0]),
                        std::floor((a_particlePosition[1] - a_probLo[1]) / a_dx[1]),
                        std::floor((a_particlePosition[2] - a_probLo[2]) / a_dx[2])));
}

inline bool
ParticleOps::domainIntersection(const RealVect& a_oldPos,
                                const RealVect& a_newPos,
                                const RealVect& a_probLo,
                                const RealVect& a_probHi,
                                Real&           a_s)
{

  // TLDR: This code does a boundary intersection test and returns where on the interval [oldPos, newPos] the intersection
  //       happened. We do this by checking if the particle moves towards a particular domain side and ends up outside of it.

  a_s = std::numeric_limits<Real>::max();

  bool crossedDomainBoundary = false;

  const RealVect path = a_newPos - a_oldPos;

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const Side::LoHiSide side      = sit();
      const RealVect       wallPoint = (side == Side::Lo) ? a_probLo : a_probHi; // A point on the domain side
      const RealVect       n0        = sign(side) *
                          RealVect(
                            BASISV(dir)); // Normal vector pointing OUT of the domain on side sit and direction dir.
      const Real normPath = PolyGeom::dot(n0, path); // Component of the path that is normal to the domain edge/face.

      // If normPath > 0 then the particle trajectory points towards the domain edge/face and we can have an intersection.
      if (normPath > 0.0) {

        // s determines the intersection point between the particle path and the plane corresponding to the domain edge/face. Note that
        // we consider the edge/face to be an infinite plane and we just compute the intersection point between each edge/face and select the
        // closest intersection point.
        const Real s = PolyGeom::dot(wallPoint - a_oldPos, n0) / normPath;
        if (s >= 0.0 && s <= 1.0) {
          crossedDomainBoundary = true;
          if (s < a_s) {
            a_s = s;
          }
        }
      }
    }
  }

  return crossedDomainBoundary;
}

inline bool
ParticleOps::ebIntersectionBisect(const RefCountedPtr<BaseIF>& a_impFunc,
                                  const RealVect&              a_oldPos,
                                  const RealVect&              a_newPos,
                                  const Real&                  a_bisectStep,
                                  Real&                        a_s)
{

  // TLDR: We compute the intersection point using a bisection algorithm. We divide the full path into intervals and check if an interval
  //       has a root. If it does, we compute it using Brent's algorithm.

  a_s = std::numeric_limits<Real>::max();

  bool crossedEB = false;

  const Real     pathLen = (a_newPos - a_oldPos).vectorLength(); // Total path len
  const int      nsteps  = ceil(pathLen / a_bisectStep);         // Number of bisection intervals
  const RealVect dxStep  = (a_newPos - a_oldPos) / nsteps;       // Physical length of each bisection interval

  // Check each interval
  RealVect curPos = a_oldPos;
  for (int istep = 0; istep < nsteps; istep++) {
    const Real fa = a_impFunc->value(curPos); // Value of the implicit function at the start of the bisection interval
    const Real fb = a_impFunc->value(curPos +
                                     dxStep); // Value of the implicit function at the end of the bisection interval

    if (fa * fb <= 0.0) {

      // If this triggered we happen to know that f(pos+dxStep) > 0.0 and f(pos) < 0.0 and so we must have a root on the interval. We now compute the precise location
      // where the particle crossed the EB. For that we use a Brent root finder on the interval [pos, pos+dxStep]. This is a 1D problem.
      const RealVect intersectionPos = PolyUtils::brentRootFinder(a_impFunc, curPos, curPos + dxStep);
      a_s                            = (intersectionPos - a_oldPos).vectorLength() / pathLen;
      crossedEB                      = true;

      break;
    }
    else { // Move to next interval
      curPos += dxStep;
    }
  }

  return crossedEB;
}

inline bool
ParticleOps::ebIntersectionRaycast(const RefCountedPtr<BaseIF>& a_impFunc,
                                   const RealVect&              a_oldPos,
                                   const RealVect&              a_newPos,
                                   const Real&                  a_tolerance,
                                   Real&                        a_s)
{

  a_s = std::numeric_limits<Real>::max();

  bool ret = false;

  // Absolute distance to EB.
  auto dist = [&](const RealVect& x) -> Real {
    return std::abs(a_impFunc->value(x));
  };

  const Real D  = (a_newPos - a_oldPos).vectorLength(); // Total particle path length
  const Real D0 = dist(a_oldPos);                       // Distance to EB from starting position

  // If the distance to the EB from the starting position is smaller than the total path length, we need to check for intersections.
  if (D > D0) {

    const RealVect t = (a_newPos - a_oldPos) / D; // Particle trajectory.

    // Move a_oldPos along +t. If we end up too close to the boundary the particle has intersected the BC. Note that this does NOT check for whether or not
    // the particle moves tangential to the EB surface. The length of each step is the distance to the EB, so if the particle is close to the EB but moves
    // tangentially to it, this routine will be EXTREMELY slow.
    RealVect xa = a_oldPos;
    Real     r  = D;
    Real     d  = dist(xa);

    while (d < r) {

      if (d < a_tolerance) { // We collided.
        a_s = (xa - a_oldPos).vectorLength() / D;
        ret = true;
        break;
      }
      else { // We did not collide.
        xa += t * d;
        r -= d;
        d = dist(xa);
      }
    }
  }

  return ret;
}

template <typename P>
inline void
ParticleOps::copy(ParticleContainer<P>& a_dst, const ParticleContainer<P>& a_src) noexcept
{
  CH_TIME("ParticleOps::copy(ParticleContainer<P> x2)");

  CH_assert(a_dst.getRealm() == a_src.getRealm());

  for (int lvl = 0; lvl <= a_dst.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_dst.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      a_dst[lvl][din].listItems() = a_src[lvl][din].listItems();
    }
  }
}

template <typename P>
inline void
ParticleOps::copyDestructive(ParticleContainer<P>& a_dst, ParticleContainer<P>& a_src) noexcept
{
  CH_TIME("ParticleOps::copyDestructive(ParticleContainer<P> x2)");

  CH_assert(a_dst.getRealm() == a_src.getRealm());

  for (int lvl = 0; lvl <= a_dst.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_dst.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      a_dst[lvl][din].listItems() = a_src[lvl][din].listItems();

      a_src[lvl][din].listItems().clear();
    }
  }
}

template <typename P, const Real& (P::*scalarQuantity)() const>
inline Real
ParticleOps::sum(const ParticleContainer<P>& a_particles) noexcept
{
  CH_TIME("ParticleOps::sum(ParticleContainer<P>)");

  Real particleSum = 0.0;

  for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : particleSum)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const List<P>& particles = a_particles[lvl][din].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        particleSum += (lit().*scalarQuantity)();
      }
    }
  }

  return ParallelOps::sum(particleSum);
}

template <typename P, Real (P::*scalarQuantity)()>
inline Real
ParticleOps::sum(const ParticleContainer<P>& a_particles) noexcept
{
  CH_TIME("ParticleOps::sum(ParticleContainer<P>)");

  Real particleSum = 0.0;

  for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : particleSum)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const List<P>& particles = a_particles[lvl][din].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        particleSum += (lit().*scalarQuantity)();
      }
    }
  }

  return ParallelOps::sum(particleSum);
}

template <typename P>
inline void
ParticleOps::removeParticles(ParticleContainer<P>&                a_particles,
                             const std::function<bool(const P&)>& a_removeCriterion) noexcept
{
  CH_TIME("ParticleOps::removeParticles(ParticleContainer<P>, std::function<bool(const P&)>)");

  for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& particles = a_particles[lvl][din].listItems();

      for (ListIterator<P> lit(particles); lit.ok();) {
        if (a_removeCriterion(lit())) {
          particles.remove(lit);
        }
        else {
          ++lit;
        }
      }
    }
  }
}

template <typename P>
inline void
ParticleOps::transferParticles(ParticleContainer<P>&                a_dstParticles,
                               ParticleContainer<P>&                a_srcParticles,
                               const std::function<bool(const P&)>& a_transferCrit) noexcept
{
  CH_TIME("ParticleOps::transferParticles(ParticleContainer<P>, ParticleContainer<P>& std::function<bool(const P&)>)");

  CH_assert(a_dstParticles.getRealm() == a_srcParticles.getRealm());
  CH_assert(a_dstParticles.getFinestLevel() == a_srcParticles.getFinestLevel());

  for (int lvl = 0; lvl <= a_srcParticles.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_srcParticles.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& dstParticles = a_dstParticles[lvl][din].listItems();
      List<P>& srcParticles = a_srcParticles[lvl][din].listItems();

      for (ListIterator<P> lit(srcParticles); lit.ok();) {
        if (a_transferCrit(lit())) {
          dstParticles.transfer(lit);
        }
        else {
          ++lit;
        }
      }
    }
  }
}

template <typename P>
inline void
ParticleOps::setData(ParticleContainer<P>& a_particles, const std::function<void(P&)>& a_functor) noexcept
{
  CH_TIME("ParticleOps::setData(ParticleContainer<P>, std::function<void(P&)>)");

  for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& particles = a_particles[lvl][din].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        a_functor(lit());
      }
    }
  }
}

template <typename P, Real& (P::*particleScalarField)()>
inline void
ParticleOps::setValue(ParticleContainer<P>& a_particles, const Real a_value) noexcept
{
  CH_TIME("ParticleOps::setValue(ParticleContainer<P>, Real");

  for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& particles = a_particles[lvl][din].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        (p.*particleScalarField)() = a_value;
      }
    }
  }
}

template <typename P, RealVect& (P::*particleVectorField)()>
inline void
ParticleOps::setValue(ParticleContainer<P>& a_particles, const RealVect a_value) noexcept
{
  CH_TIME("ParticleOps::setValue(ParticleContainer<P>, RealVect");

  for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& particles = a_particles[lvl][din].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        (p.*particleVectorField)() = a_value;
      }
    }
  }
}

#ifdef CH_MPI
template <typename P>
inline void
ParticleOps::scatterParticles(
  std::map<std::pair<unsigned int, unsigned int>, List<P>>&              a_receivedParticles,
  std::vector<std::map<std::pair<unsigned int, unsigned int>, List<P>>>& a_sentParticles) noexcept
{
  CH_TIME("ParticleOps::scatterParticles");

  const int numRanks = numProc();

  int mpiErr;

  CH_assert(a_sentParticles.size() == numProc());

  const size_t linearSize = P().size();

  std::vector<unsigned int> sendSizes(numProc(), 0);
  std::vector<unsigned int> recvSizes(numProc(), 0);

  // Figure out the size of a_sentParticles[irank]
  for (int irank = 0; irank < numRanks; irank++) {
    sendSizes[irank] = 0;

    // Figure out the message size sent from this rank.
    const std::map<std::pair<unsigned int, unsigned int>, List<P>>& particlesToRank = a_sentParticles[irank];

    // Size of particles along.
    for (const auto& m : particlesToRank) {
      sendSizes[irank] += m.second.length();
    }
    sendSizes[irank] *= linearSize;

    // Add in the map keys and the list length
    sendSizes[irank] += particlesToRank.size() * (2 * sizeof(unsigned int) + sizeof(size_t));
  }

  mpiErr = MPI_Alltoall(&sendSizes[0], 1, MPI_UNSIGNED, &recvSizes[0], 1, MPI_UNSIGNED, Chombo_MPI::comm);
  if (mpiErr != MPI_SUCCESS) {
    MayDay::Error("ParticleOps::scatterParticles - MPI communication error in MPI_AlltoAll");
  }

  // Sanity check
  if (recvSizes[procID()] != 0) {
    MayDay::Error("ParticleOps::scatterParticles - 'recvSizes[procID()] != 0' failed");
  }

  std::vector<MPI_Request> recvReq(numRanks - 1);
  std::vector<char*>       recvBuf(numRanks, nullptr);

  std::vector<MPI_Request> sendReq(numRanks - 1);
  std::vector<char*>       sendBuf(numRanks, nullptr);

  int recvReqCount = 0;
  int sendReqCount = 0;

  for (int irank = 0; irank < numRanks; irank++) {
    if (recvSizes[irank] > 0) {
      recvBuf[irank] = new char[recvSizes[irank]];

      MPI_Irecv(recvBuf[irank], recvSizes[irank], MPI_CHAR, irank, 1, Chombo_MPI::comm, &recvReq[recvReqCount]);

      recvReqCount++;
    }
  }

  // Pack data into buffers
  for (int irank = 0; irank < numRanks; irank++) {
    if (sendSizes[irank] > 0) {
      sendBuf[irank] = new char[sendSizes[irank]];

      char* data = sendBuf[irank];

      const std::map<std::pair<unsigned int, unsigned int>, List<P>>& particles = a_sentParticles[irank];

      // Linearize each entry in the map onto the send buffer. We encode this by (level, index, list length, List<P>)
      for (const auto& cur : particles) {

        // and grid index
        *((unsigned int*)data) = cur.first.first;
        data += sizeof(unsigned int);
        *((unsigned int*)data) = cur.first.second;
        data += sizeof(unsigned int);

        // List length aka number of particles
        *((size_t*)data) = cur.second.length();
        data += sizeof(size_t);

        // Linearize the particle list
        for (ListIterator<P> lit(cur.second); lit.ok(); ++lit) {
          const P& p = lit();

          p.linearOut((void*)data);
          data += linearSize;
        }
      }

      // Send from me to the other rank.
      MPI_Isend(sendBuf[irank], sendSizes[irank], MPI_CHAR, irank, 1, Chombo_MPI::comm, &sendReq[sendReqCount]);

      sendReqCount++;
    }
  }

  std::vector<MPI_Status> recvStatus(numRanks - 1);
  std::vector<MPI_Status> sendStatus(numRanks - 1);

  if (sendReqCount > 0) {
    mpiErr = MPI_Waitall(sendReqCount, &sendReq[0], &sendStatus[0]);

    if (mpiErr != MPI_SUCCESS) {
      MayDay::Error("ParticleOps::scatterParticles: send communication failed");
    }
  }

  if (recvReqCount > 0) {
    mpiErr = MPI_Waitall(recvReqCount, &recvReq[0], &recvStatus[0]);

    if (mpiErr != MPI_SUCCESS) {
      MayDay::Error("ParticleOps::scatterParticles: receive communication failed");
    }
  }

  // Unpack the send buffers
  for (int irank = 0; irank < numRanks; irank++) {
    P p;

    if (recvSizes[irank] > 0) {

      char* data = recvBuf[irank];

      unsigned int in = 0;

      while (in < recvSizes[irank]) {

        // Level and grid index
        unsigned int lvl = *((unsigned int*)data);
        data += sizeof(unsigned int);

        unsigned int idx = *((unsigned int*)data);
        data += sizeof(unsigned int);

        // List length aka number of particles
        size_t numParticles = *((size_t*)data);
        data += sizeof(size_t);

        for (size_t ipart = 0; ipart < numParticles; ipart++) {
          p.linearIn((void*)data);
          data += linearSize;

          a_receivedParticles[std::pair(lvl, idx)].add(p);
        }

        in += 2 * sizeof(unsigned int) + sizeof(size_t) + numParticles * linearSize;
      }
    }
  }

  // Delete buffers
  for (int irank = 0; irank < numRanks; irank++) {
    a_sentParticles[irank].clear();

    if (sendSizes[irank] > 0) {
      delete[] sendBuf[irank];
    }
    if (recvSizes[irank] > 0) {
      delete[] recvBuf[irank];
    }
  }
}
#endif

#include <CD_NamespaceFooter.H>

#endif
