/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ParticleMergePartitioners.H
  @brief  Implementation of various partitioners for superparticle handling of ItoSolver
  @author Robert Marskar
*/

#ifndef CD_ParticleMergePartitioners_H
#define CD_ParticleMergePartitioners_H

// Std includes
#include <limits>

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_ParticleMerge.H>
#include <CD_NamespaceHeader.H>

namespace ParticleMerge {

  /*!
    @brief Partitioning function which partitions masses equally in the kd-tree. 
  */
  template <class T>
  NodePartitioner<T> NodePartitionEqualMass = [](std::shared_ptr<Node<T>>& a_node, const int a_dir) {
    CH_TIME("NodePartitionEqualMass<T>");

    CH_assert(a_node != nullptr);
    CH_assert(a_dir >= 0);
    CH_assert(a_dir < SpaceDim);

    // Get a handle to the node data.
    std::vector<T>& data = a_node->getData();
    const Real      mass = a_node->mass();

    // Some basic assertions -- the loop should evaluate to an empty loop when DEBUG=TRUE.
    CH_assert(data.size() > 0);
    CH_assert(mass > 1.0);

    for (const auto& n : data) {
      CH_assert(n.mass() >= 1.0);
    }

    // TLDR: This function splits the primitives in the node into two subvolumes. A splitting plane is made along the input coordinate and it is positioned
    //       such that the total mass on the left/right side of the splitting plane differ by at most one computational particle. The guts of the code mainly
    //       consists of tests for how to divide the "median" particle -- there are various cases that we need to cover here (such as the median particle being
    //       a physical particle or a superparticle).

    // 1. Sort data along the input coordinate. We need this because we need to have something sensible to iterate through
    //    when we place the splitting plane.
    std::sort(std::begin(data), std::end(data), [a_dir](const T& p1, const T& p2) { return p1[a_dir] < p2[a_dir]; });

    // 2. Figure out where to place the splitting plane -- this is not the physical coordinate but index in data (the list of primitives). We are
    //    actually looking for the particle on some median where the masses between the two half-planes can be made approximately the same. We begin
    //    by setting the median on i = 0 in the primitives vector (data). We then go through the primitives and check if adding the median particle
    //    ends up placing MORE mass in the left node than the right node. Once that happens we have the inflection point where we need to split. The median
    //    particle will be assigned later because it may be a superparticle and we can equilibriate the masses even further.
    //
    //
    // In all of the below. massLeft and massRight are the masses on each side of the median particle. As remarked above, the median particle could be split
    // later.
    int  medianParticleIndex = 0;
    Real massLeft            = 0.0;                           // Zero mass, as it should be.
    Real massRight = mass - data[medianParticleIndex].mass(); // So, entire mass minus the median particle mass.

    for (int i = 1; i < data.size(); i++) {
      const Real& medianParticleMass = data[medianParticleIndex].mass(); // This is the current median particle mass.

      // If adding the median particle to the mass in the left node still yields a lower mass than the right now, we just shift the median particle index. So, the previous
      // "median particle" goes to the left node but we then have a new median particle.
      if (massLeft + medianParticleMass < massRight) {
        massLeft += medianParticleMass;

        medianParticleIndex = i;

        massRight = mass - massLeft -
                    data[medianParticleIndex].mass(); // Right mass = Total mass - left mass - median particle mass
      }
      else {
        break;
      }
    }

    // 2. Make the two leaves and push data into them. Currently, we know that all particles to the left of the "median" particle go into the left half-plane
    //    and all particles to the right of the "median" particle go into the right half-plane. We do not yet know how to place the median particle -- it could
    //    go into either the left node or the right node, or it could be split into and placed into either node.
    std::shared_ptr<Node<T>>& left  = a_node->getLeft();
    std::shared_ptr<Node<T>>& right = a_node->getRight();

    if (!left)
      left = std::make_shared<Node<T>>();
    if (!right)
      right = std::make_shared<Node<T>>();

    std::vector<T>& leftData  = left->getData();
    std::vector<T>& rightData = right->getData();

    leftData.assign(data.begin(), data.begin() + medianParticleIndex);
    rightData.assign(data.begin() + medianParticleIndex + 1, data.end());

    // 3. We next need to assign the median particle. If we can split up the median point into several ones, we do it. This would be the case if the median point
    //    is a superparticle; in this case we divide the weights as best we can and then assign them to the left/right halves. If we can't then we just assign the
    //    particle to whichever left/right child node has the lowest mass.
    const T& medianParticle = data[medianParticleIndex];

    if (medianParticle.canSplit()) {
      const RealVect& medianPos  = medianParticle.pos();  // Physical position
      const Real&     medianMass = medianParticle.mass(); // Always positive
      const Real&     medianEnergy =
        medianParticle.energy(); // Note: This is average energy so it remain the same for subparticles.
      const Real massDiff = massRight - massLeft; // Not necessarily positive.

      // This hook is for when all mass goes into one part.
      if (medianMass < std::abs(massDiff)) {
        if (massDiff > 0) { // Right mass is bigger than left so put entire particle in left mass.
          massLeft += medianMass;
          leftData.emplace_back(T(medianPos, medianMass, medianEnergy));
        }
        else {                     // Left mass is bigger than right so put everything into the right node.
          massRight += medianMass; // Goes into right
          rightData.emplace_back(T(medianPos, medianMass, medianEnergy));
        }
      }
      else {
        // Distribute mass as best we can. This might seem weird because we can have massLeft > massRight in which
        // case massDiff < 0. This still works since we are simply subtracting the difference from the left leaf first and
        // then divide a larger total mass between the two leaves. Essentially, we are subtracting the part of the median
        // particle that gives massLeft > massRight from the left particle first, and then we reassign a larger total mass.
        Real tmpMassLeft  = massDiff;
        Real tmpMassRight = 0.0;
        Real tmpMassRem   = medianMass - massDiff;

        const long long N = llround(tmpMassRem);
        if (N > 0LL) {
          const long long NR = N / 2;
          const long long NL = N - NR;

          tmpMassLeft += (tmpMassRem / N) * NL;
          tmpMassRight += (tmpMassRem / N) * NR;
        }

        if (tmpMassLeft > 0.0) {
          leftData.emplace_back(T(medianPos, tmpMassLeft, medianEnergy));
          massLeft += tmpMassLeft;
        }

        if (tmpMassRight > 0.0) {
          rightData.emplace_back(T(medianPos, tmpMassRight, medianEnergy));
          massRight += tmpMassRight;
        }
      }
    }
    else {
      // Not a computational particle -- put the particle in the child node with the lowest mass.
      if (massLeft <= massRight) {
        leftData.emplace_back(medianParticle);
        massLeft += medianParticle.mass();
      }
      else {
        rightData.emplace_back(medianParticle);
        massRight += medianParticle.mass();
      }
    }

    // Update the tree structure.
    left->setMass(massLeft);
    right->setMass(massRight);

    left->setParent(a_node);
    right->setParent(a_node);

    data.resize(0);

    // Some assertions in case this code every breaks.
    CH_assert((std::abs(mass - (massLeft + massRight)) <= 1.E-12));
    CH_assert(massLeft >= std::numeric_limits<Real>::min());
    CH_assert(massRight >= std::numeric_limits<Real>::min());

    // These should be evaluated to empty loopps by the compiler if DEBUG=FALSE.
    for (const auto& l : leftData) {
      CH_assert(l.mass() >= std::numeric_limits<Real>::min());
    }

    for (const auto& r : rightData) {
      CH_assert(r.mass() >= std::numeric_limits<Real>::min());
    }
  };
} // namespace ParticleMerge

#include <CD_NamespaceFooter.H>

#endif
