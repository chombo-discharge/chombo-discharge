#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "localfunctions.H"
#include "parstream.H"

#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "BRMeshRefine.H"
#include "FineInterp.H"


// set domain boundary conditions from input file
int setDomainBC(DomainGhostBC& domghostbc,
                const int* ibclo,
                const int* ibchi,
                const bool verbose)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    //lo bc
    {
      NeumannBC neumbc(idir, Side::Lo);
      DirichletBC dircbc(idir, Side::Lo);
      if (ibclo[idir] == 0)
      {
        if (verbose&& procID()==0)
          pout() << "dirchelet bcs in low direction for side " << idir << endl;
        domghostbc.setBoxGhostBC(dircbc);
      }
      else if (ibclo[idir] == 1)
      {
        if (verbose&& procID()==0)
          pout() << "neumann bcs in low direction for side " << idir << endl;
        domghostbc.setBoxGhostBC(neumbc);
      }
      else
      {
        if (verbose)
        {
          cerr << "setDomainBC:: bogus input bc_lo flag" << endl;
        }
        return(1);
      }
    }
    //hi bc
    {
      NeumannBC neumbc(idir, Side::Hi);
      DirichletBC dircbc(idir, Side::Hi);
      if (ibchi[idir] == 0)
      {
        if (verbose&& procID()==0)
          pout() << "dirchelet bcs in high direction for side " << idir << endl;
        domghostbc.setBoxGhostBC(dircbc);
      }
      else if (ibchi[idir] == 1)
      {
        if (verbose&& procID()==0)
          pout() << "neumann bcs in high direction for side " << idir << endl;
        domghostbc.setBoxGhostBC(neumbc);
      }
      else
      {
        if (verbose)
          cerr << "setDomainBC:: bogus input bc_hi flag" << endl;
        return(2);
      }
    }
  }
  return(0);
}


