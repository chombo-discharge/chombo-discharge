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

class WilliamsonTestCase2
{
  public:
    // Constructor
    WilliamsonTestCase2() :
      m_dAlpha(0.0),
      m_dH0(4.7057e-4),
      m_dU0(M_PI / 6.0),
      m_dOmega(7.292e-5 * 86400.0),
      m_dG(9.80616 * 86400.0 * 86400.0 / 6.37122e6)
    { }

  public:
    // Compute the height and momentum field at the given point in
    // mapped coordinates.
    void calculateState(
      Real dAlpha,
      Real dBeta,
      int nPanel,
      Real& dH,
      Real& dHU,
      Real& dHV
    )
    {
      CH_assert((nPanel >= 0) && (nPanel <= 5));

      // Convert to Gnomonic coordinates
      Real dX = tan(dAlpha);
      Real dY = tan(dBeta);

      Real dDelta = sqrt(1.0 + dX * dX + dY * dY);

      Real dA;

      switch(nPanel)
      {
        case 0:
          dA = dY / dDelta * cos(m_dAlpha)
             - 1.0 / dDelta * sin(m_dAlpha);
          break;

        case 1:
          dA = dY / dDelta * cos(m_dAlpha)
             + dX / dDelta * sin(m_dAlpha);
          break;

        case 2:
          dA = dY / dDelta * cos(m_dAlpha)
             + 1.0 / dDelta * sin(m_dAlpha);
          break;

        case 3:
          dA = dY / dDelta * cos(m_dAlpha)
             - dX / dDelta * sin(m_dAlpha);
          break;

        case 4:
          dA = 1.0 / dDelta * cos(m_dAlpha)
             + dY / dDelta * sin(m_dAlpha);
          break;

        case 5:
          dA = - 1.0 / dDelta * cos(m_dAlpha)
               - dY / dDelta * sin(m_dAlpha);
          break;
      }

      dH = (m_dH0 -
        (m_dOmega * m_dU0 + 0.5 * m_dU0 * m_dU0) / m_dG * (dA * dA)) /
        (dDelta * dDelta * dDelta);

      // These quantities should be computed from the position in RLL
      // coordinates
      dHU = 0.0;
      dHV = 0.0;

    }

    // Compute the height and momentum field over the given FArrayBox.
    void calculateState(FArrayBox& data,
                        const NodeFArrayBox& Xi,
                        const Box& a_box,
                        int nPanel)
    {
      BoxIterator bit(a_box);

      for (bit.begin(); bit.ok(); ++bit)
      {
        IntVect iv[2];
        iv[0] = bit();
        iv[1] = iv[0] + BASISV(0) + BASISV(1);

        // For now use the center-point of a given element
        calculateState(
          0.5 * (Xi(iv[0],0) + Xi(iv[1],0)),
          0.5 * (Xi(iv[0],1) + Xi(iv[1],1)),
          nPanel,
          data(iv[0],0),
          data(iv[0],1),
          data(iv[0],2));
      }
    }

  public:
    // Rotation angle of grid with respect to the equator
    Real m_dAlpha;

    // Background height field (in Earth radii)
    Real m_dH0;

    // Background reference flow velocity (in Earth radii / day)
    Real m_dU0;

    // Angular rotation rate of the Earth (in 1 / day)
    Real m_dOmega;

    // Gravitational constant (in Earth radii / day^2)
    Real m_dG;
};

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


// I/O function --
int writeCubedSphereIO(LevelData<FArrayBox>& a_data,
                       MultiBlockCoordSys& a_MBCS,
                       const string& a_fileRoot)
{
  const NewCoordSys* baseCoordSysPtr = a_MBCS.getCoordSys(0);
  int dimension = baseCoordSysPtr->dimension();

  // first, write out data to "regular" hdf5 file:
  char iter_str[80];

  sprintf(iter_str, "%s%dd.hdf5", a_fileRoot.c_str(), SpaceDim);
  string dataFileName(iter_str);

  RealVect dx = baseCoordSysPtr->dx();
  // construct a domain box which encloses all of the boxes
  const DisjointBoxLayout& grids = a_data.getBoxes();

  // this is silly, but quick and easy
  IntVectSet gridIVS;
  LayoutIterator lit = grids.layoutIterator();
  for (lit.begin(); lit.ok(); ++lit)
    {
      gridIVS |= grids[lit];
    }
  Box domainBox = gridIVS.minBox();

#ifdef CH_USE_HDF5
  WriteUGHDF5(dataFileName, grids,
              a_data, domainBox);
#endif

  // now create node-centered data for geometric info
  Vector<LevelData<NodeFArrayBox>* > vectNodeLoc(1, NULL);

  // use same ghosting as cell-centered data used
  IntVect ghostVect = a_data.ghostVect();
  vectNodeLoc[0] = new LevelData<NodeFArrayBox>(grids,
                                                dimension,
                                                ghostVect);

  LevelData<NodeFArrayBox>& levelNodeData = *vectNodeLoc[0];

  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& gridBox = grids[dit];
      // this is gonna be slow, since we're doing it node-by-node
      NodeFArrayBox& thisNodeFAB = levelNodeData[dit];
      // node-centered FAB
      FArrayBox& thisFAB = thisNodeFAB.getFab();

      const NewCoordSys* thisCoordSys = a_MBCS.getCoordSys(gridBox);
      CH_assert(thisCoordSys);
      thisCoordSys->getNodeRealCoordinates(thisFAB, thisFAB.box());


    } // end loop over grids



  // create names
  Vector<string> locationNames(dimension);
  for (int dir=0; dir<dimension; dir++)
    {
      if (dir==0)
        {
          locationNames[0] = "x";
        }
      else if (dir == 1)
        {
          locationNames[1] = "y";
        }
      else if (dir == 2)
        {
          locationNames[2] = "z";
        }
      else if (dir == 3)
        {
          locationNames[3] = "u";
        }
      else if (dir == 4)
        {
          locationNames[4] = "v";
        }
      else if (dir == 5)
        {
          locationNames[5] = "w";
        }
    }

  sprintf(iter_str, "%s%dd.map.hdf5", a_fileRoot.c_str(), SpaceDim);
  string gridInfoFileName(iter_str);

  int numLevels = 1;

#ifdef CH_USE_HDF5
  // set up for NodeAMRIO call
  Vector<DisjointBoxLayout> vectGrids(1,grids);
  Real dt = 0;
  Real time = 0;
  Vector<int> vectRatio(1,1);

  // now call nodal WriteAMRHierarchy function...
  WriteAMRHierarchyHDF5(gridInfoFileName,
                        vectGrids,
                        vectNodeLoc,
                        locationNames,
                        domainBox,
                        dx[0],
                        dt,
                        time,
                        vectRatio,
                        numLevels);
#endif

  // clean up after ourselves here
  for (int level=0; level<numLevels; level++)
    {
      if (vectNodeLoc[level] != NULL)
        {
          delete vectNodeLoc[level];
          vectNodeLoc[level] = NULL;
        }
    }

  return 0;
}

int
testCubedSphere2D()
{
  int returnCode = 0;
#if CH_SPACEDIM >= 2 && CH_SPACEDIM <= 3
  // Resolution of each panel on the cubed sphere
  //static const int nResolution = 40;
  static const int nResolution = 20;

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

  // Construct the williamson test case 2
  WilliamsonTestCase2 test;

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

      test.calculateState(data[dit], vectNodeLoc[dit], gridBox, blockNumber);
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
                        domainBox,0);

      //writeCubedSphereIO(data, cubedSphereCS, fileRoot );
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

