#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ZeldovichIBC.H"
#include "ParmParse.H"

ZeldovichIBC::ZeldovichIBC()
{}

ZeldovichIBC::ZeldovichIBC(const ProblemDomain& a_domain,
                           const Real& a_dx,
                           const RealVect& a_offset,
                           const Real& a_initialA,
                           const Real& a_collapseA,
                           const Real& a_partPerCell,
                           const Real& a_wavenum,
                           const Cosmology& a_cosmology)
{
  define(a_domain,
         a_dx,
         a_offset,
         a_initialA,
         a_collapseA,
         a_partPerCell,
         a_wavenum,
         a_cosmology);
}

void ZeldovichIBC::define(const ProblemDomain& a_domain,
                          const Real& a_dx,
                          const RealVect& a_offset,
                          const Real& a_initialA,
                          const Real& a_collapseA,
                          const Real& a_partPerCell,
                          const Real& a_wavenum,
                          const Cosmology& a_cosmology)
{
  m_domain = a_domain;
  m_dx = a_dx;
  m_offset = a_offset;
  m_initialA = a_initialA;
  m_collapseA = a_collapseA;
  m_partPerCell = a_partPerCell;
  m_wavenum = a_wavenum;
  m_cosmology = a_cosmology;
}

ZeldovichIBC* ZeldovichIBC::new_zeldovichIBC() const
{
  ZeldovichIBC* retval = new ZeldovichIBC(m_domain,
                                          m_dx,
                                          m_offset,
                                          m_initialA,
                                          m_collapseA,
                                          m_partPerCell,
                                          m_wavenum,
                                          m_cosmology);
  
  return retval;
}

void ZeldovichIBC::initialize(ParticleData<Particle>& a_data)
{

  // compute domain extent
  IntVect domainDimensions = m_domain.size();
  RealVect domainLeftEdge = m_offset;
  RealVect domainRightEdge = ((RealVect)domainDimensions)*m_dx;
  RealVect Lbox = domainRightEdge - domainLeftEdge;
  
  // figure out number of particles and total simulation volume
  Real numParticles = 1.0;
  Real volume = 1.0;
  for (int dim = 0; dim < CH_SPACEDIM; dim++)
    {
      numParticles *= m_partPerCell*domainDimensions[dim];
      volume *= Lbox[dim];
    }
  
  // setup the orientation, wave number, and amplitude of the 
  // perturbation here.
  ParmParse pp;
  Vector<int> n;
  pp.getarr("wavenums", n, 0, CH_SPACEDIM);

  RealVect orientation = RealVect(D_DECL((Real) n[0], (Real) n[1], (Real) n[2]));
  Real orientationMag  = sqrt((orientation*orientation).sum());
  RealVect khat = orientation / orientationMag;
  RealVect k = m_wavenum * orientation;
  Real kMag = sqrt((k*k).sum());
  Real amplitude = 1.0 / kMag / m_collapseA;
  
  // the initial proper time and adot (time derivate of expansion factor)
  Real tInitial = m_cosmology.t(m_initialA);
  Real adot     = m_cosmology.adot(tInitial);
  
  // mass of each particle. Chosen so that the mean density is 1.0
  Real mp = volume / numParticles;
  
  // Each proc loops over its own boxes, setting the initial
  // particle positions in each cell.
  const BoxLayout& BL = a_data.getBoxes();
  DataIterator dit(BL);
  
  for (dit.begin(); dit.ok(); ++dit)
    {
      CH_XD::List<Particle> thisList;
      Real dx = m_dx;
      
      const Box thisBox = BL.get(dit);
      Box partBox = thisBox;
      
      // coarsen or refine box so we can iterate over cell centers
      // and get the right number of particles per cell
      if (m_partPerCell < 1.0)
        {
          const int partRef = (int)(1.0 / m_partPerCell);
          partBox.coarsen(partRef);
          dx /= m_partPerCell;
        }
      else if (m_partPerCell > 1.0)
        {
          const int partRef = (int)(m_partPerCell);
          partBox.refine(partRef);
          dx /= m_partPerCell;
        }
      
      BoxIterator bit(partBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          // index of this cell
          IntVect iv = bit();
          
          // Lagrangian coordinate
          RealVect q = m_offset + ((RealVect)iv + 0.5)*dx;
          
          // compute dot product of k and q
          Real kq = 0.0;
          for (int dim=0; dim < CH_SPACEDIM; dim++)
            {
              kq += k[dim] * q[dim];
            }
          
          // solve for the initial particle position and velocity
          RealVect x = q + amplitude * m_initialA * sin(kq) * khat;
          RealVect v = amplitude * m_initialA * adot * sin(kq) * khat;
          
          // enforce periodic boundary conditions
          for (int dim=0; dim < CH_SPACEDIM; ++dim)
            {
              if (x[dim] < domainLeftEdge[dim])
                {
                  x[dim] += Lbox[dim];
                }
              if (x[dim] > domainRightEdge[dim])
                {
                  x[dim] -= Lbox[dim];
                }
            }
          
          // create this particle and append to list
          Particle particle(mp, x, v);
          thisList.append(particle);
        }

      // finally, add particles destructively to this ListBox. Those that are
      // left behind are outcasts.
      a_data[dit].addItemsDestructive(thisList);
      
      // If all the particles were contained on the Box that 
      // created them, thisList would now be empty. But, 
      // because we perturb the initial particle distribution,
      // particles can move off the box that created them. We 
      // add these to the outcast list to be re-distributed later. 
      a_data.outcast().catenate(thisList);
    }
}
