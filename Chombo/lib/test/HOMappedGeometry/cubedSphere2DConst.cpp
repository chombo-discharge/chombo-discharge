#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>
#include "mappedGridIO.H"
#include "newMappedGridIO.H"
#include "UGIO.H"
using std::cout;
using std::endl;
using std::ofstream;

#include "parstream.H"
#include "ParmParse.H"
#include "CubedSphere2DCS.H"
#include "CubedSphere2DPanelCS.H"
#include "BRMeshRefine.H"
#include "BoxIterator.H"
#include "computeNorm.H"
#include "FABView.H"
#include "LoadBalance.H"
#include <cfloat>

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testCubedSphere" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = false ;
static bool writePlotFiles = true;
static bool writeGridsFile = true;

// #ifdef CH_USE_DOUBLE
// static Real precision = 1.0e-15;
// #else
// static Real precision = 1.0e-7;
// #endif

///
// Parse the standard test options (-v -q -h) and
// app-specific options (-S <domain_size>) out of the command line
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
    {
      if ( argv[i][0] == '-' ) //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
            {
              verbose = true ;
              //argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              //argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-h" ,3 ) == 0 )
            {
              pout() << "usage: " << pgmname << " [-hqv]" << std::endl ;
              exit( 99 ) ;
            }
        }
    }
  return;
}


void
writeMeshToFile(char* a_filestr,
                MultiBlockCoordSys& a_MBCS,
                const DisjointBoxLayout& a_grids)
{
  ofstream os(a_filestr, ios::out);

  if (os.fail())
    {
      pout() << "cannot open grid output file " << a_filestr << endl;
      MayDay::Error();
    }

  // now create node-centered data for geometric info
  const NewCoordSys* baseCoordSysPtr = a_MBCS.getCoordSys(0);
  int dimension = baseCoordSysPtr->dimension();
  // no ghost cells here
  IntVect ghostVect = IntVect::Zero;
  LevelData<NodeFArrayBox> nodeLoc(a_grids,
                                   dimension,
                                   ghostVect);

  DataIterator dit = a_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& gridBox = a_grids[dit];
      NodeFArrayBox& thisNodeFAB = nodeLoc[dit];
      // node-centered FAB
      FArrayBox& thisFAB = thisNodeFAB.getFab();

      const NewCoordSys* thisCoordSys = a_MBCS.getCoordSys(gridBox);
      CH_assert(thisCoordSys);
      thisCoordSys->getNodeRealCoordinates(thisFAB, thisFAB.box());

      IntVect loVect = gridBox.smallEnd();
      IntVect hiVect = gridBox.bigEnd();
      // account for cell centering by adding the Unit vector
      hiVect += IntVect::Unit;

      int ilo = loVect[0];
      int ihi = hiVect[0];
      int jlo = loVect[1];
      int jhi = hiVect[1];

      IntVect iv(loVect);
      for (int i=ilo; i<= ihi; i++)
        {
          iv[0] = i;

          for (int j=jlo; j<=jhi; j++)
            {
              iv[1] = j;

              os << thisFAB(iv,0) << "  "
                 << thisFAB(iv,1) << "  "
                 << thisFAB(iv,2) << "  "
                 << "   "  << endl;
            } // end loop over j

          // add a blank line between coordinate lines
          os << endl;
        } // end loop over i
      // end x-direction

      os << endl;
      // theta-direction
      for (int j=jlo; j<=jhi; j++)
        {
          iv[1] = j;

          for (int i=ilo; i<= ihi; i++)
            {
              iv[0] = i;

              os << thisFAB(iv,0) << "  "
                 << thisFAB(iv,1) << "  "
                 << thisFAB(iv,2) << "  "
                 << "   "  << endl;
            } // end loop over i

          // add a blank line between coordinate lines
          os << endl;
        } // end loop over j
      // end y-direction


      os << endl;
    } // end loop over boxes

  os.close();

} // end if we're writing out the mesh

int
testCubedSphere2D()
{
  int returnCode = 0;
#if CH_SPACEDIM >= 2 && CH_SPACEDIM <= 3
  // Resolution of each panel on the cubed sphere
  //static const int nResolution = 40;
  static const int nResolution = 8;

  // Construct a new cubed sphere
  CubedSphere2DCS cubedSphereCS(nResolution);
  Real rMin = 0.3;
  Real rMax = 1.0;
  cubedSphereCS.setRadii(rMin, rMax);

  Vector<Box> mappingBlocks = cubedSphereCS.mappingBlocks();

  Vector<int> procAssign(mappingBlocks.size(), 0);

#ifdef CH_MPI
  LoadBalance(procAssign, mappingBlocks);
#endif

  DisjointBoxLayout grids(mappingBlocks, procAssign);

  // Three components of data
  LevelData<FArrayBox> data(grids, 3);

  // now create node-centered data for geometric info
  const NewCoordSys *pBaseCoordSys = cubedSphereCS.getCoordSys(0);

  CH_assert(pBaseCoordSys);

  int dimension = pBaseCoordSys->dimension();

  IntVect ghostVect = data.ghostVect();

  LevelData<NodeFArrayBox> vectNodeLoc(grids, dimension, ghostVect);

  Vector<Real> blockVal(6);
  // Pick these to make colors look good.
  blockVal[0] = 0.;
  blockVal[1] = 1.25;
  blockVal[2] = 2.5;
  blockVal[3] = 3.75;
  blockVal[4] = 4.1;
  blockVal[5] = 5.;

  // for now, initialize data to block number
  DataIterator dit = data.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& gridBox = grids[dit];
      int blockNumber = cubedSphereCS.whichBlock(gridBox);

      const CubedSphere2DPanelCS *pPanel =
        dynamic_cast<const CubedSphere2DPanelCS*>(
          cubedSphereCS.getCoordSys(blockNumber));

      CH_assert(pPanel);

      FArrayBox& nodeFAB = vectNodeLoc[dit].getFab();
      pPanel->getNodeMappedCoordinates(nodeFAB, nodeFAB.box());

      data[dit].setVal(blockVal[blockNumber]);
    }

  // write grids file
  if (writeGridsFile)
    {
      char filestr[80];
      int level = 0;
      sprintf(filestr, "grid%d.out", level);
      writeMeshToFile(filestr,cubedSphereCS, grids);
    }

  // call IO function
  if (writePlotFiles)
    {
#ifdef CH_USE_HDF5
      string fileRoot("cubedSphere2D.");

      // bit of a hack to get a reasonable domainBox
      // this is silly, but quick and easy
      IntVectSet gridIVS;
      LayoutIterator lit = grids.layoutIterator();
      for (lit.begin(); lit.ok(); ++lit)
        {
          gridIVS |= grids[lit];
        }
      Box domainBox = gridIVS.minBox();

      WriteMappedUGHDF5(fileRoot,
                        grids,
                        data,
                        cubedSphereCS,
                        domainBox,1);
#endif
    }
#endif
  return returnCode;
}

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif

  parseTestOptions(argc, argv) ;

  if ( verbose )
  {
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl;
  }

  int status = testCubedSphere2D();

  if (status == 0)
  {
    cout << indent << "CubedSphere2DCS test" << " passed." << endl ;
  }
  else
  {
    cout << indent << "CubedSphere2DCS test" << " failed with return code " << status << endl ;
  }

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize ();
#endif

  return status ;
}

