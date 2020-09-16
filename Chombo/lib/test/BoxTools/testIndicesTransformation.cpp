#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <string>
#include <iostream>
using std::string;
using std::endl;

#include "parstream.H"
#include "ParmParse.H"
#include "IndicesFunctions.H"
#include "IndicesTransformation.H"
#include "BoxIterator.H"
#include "FABView.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testIndicesTransformation" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = false ;
static bool dataTransformations = false ;

///
// Parse the standard test options (-v -q -h) out of the command line.
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
              pout() << "usage: " << pgmname << " [-hqvd]" << std::endl ;
              pout() << "-h: help" << endl;
              pout() << "-q: quiet" << endl;
              pout() << "-v: verbose" << endl;
              pout() << "-d: include data transforms (slow)" << endl;
              exit( 99 ) ;
            }
          if ( strncmp( argv[i] ,"-d" ,3 ) == 0 )
            {
              dataTransformations = true ;
              //argv[i] = "" ;
            }
        }
    }
  return;
}


int
testAllIndicesTransformation()
{
  int returnCode = 0;

  /*
    First get a list of transformations with
    all possible permutations,
    all possible signs,
    and several different translations.
   */

  // allPermsVect will hold all possible permutation vectors.
  Vector<IntVect> allPermsVect;
  Box allPermsBox(IntVect::Zero, (SpaceDim-1)*IntVect::Unit);
  for (BoxIterator bit(allPermsBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (isPermutationVect(iv))
        {
          allPermsVect.push_back(iv);
        }
    }

  // allSignsVect will hold all possible sign vectors.
  Vector<IntVect> allSignsVect;
  Box allSignsBox(-IntVect::Unit, IntVect::Unit);
  for (BoxIterator bit(allSignsBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      if (isSignVect(iv))
        {
          allSignsVect.push_back(iv);
        }
    }

  IntVect translators = IntVect(D_DECL6(7, 11, 13, 17, 19, 23));
  // allTranslationsVect will hold all translations we will use.
  Vector<IntVect> allTranslationsVect;
  Box translationIndexBox(-SpaceDim*IntVect::Unit, SpaceDim*IntVect::Unit);
  for (BoxIterator bit(translationIndexBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      // Check that no component of iv is 0.
      IntVect ivAbs = absolute(iv);
      if (ivAbs > IntVect::Zero)
        {
          IntVect ivSign = iv / ivAbs;
          IntVect translation;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              // index into translators
              int ind = ivAbs[idir] - 1;
              translation[idir] = ivSign[idir] * translators[ind];
            }
          allTranslationsVect.push_back(translation);
        }
    }

  int numPerms = allPermsVect.size();
  int numSigns = allSignsVect.size();
  int numTranslations = allTranslationsVect.size();

  Vector<IndicesTransformation> allTransformsVect;
  for (int iPerm = 0; iPerm < numPerms; iPerm++)
    {
      for (int iSign = 0; iSign < numSigns; iSign++)
        {
          for (int iTrans = 0; iTrans < numTranslations; iTrans++)
            {
              IntVect perm = allPermsVect[iPerm];
              IntVect sign = allSignsVect[iSign];
              IntVect translation = allTranslationsVect[iTrans];
              IndicesTransformation it(perm, sign, translation);
              allTransformsVect.push_back(it);
            }
        }
    }
  int numTransformations = allTransformsVect.size();

  if (verbose)
    {
      pout() << "Count " << numTransformations << " transformations: "
             << numPerms << " permutations, "
             << numSigns << " signs, "
             << numTranslations << " translations."
             << endl;
    }

  /*
    Check that composition with inverse is identity.
  */

  for (int iTransform = 0; iTransform < numTransformations; iTransform++)
    {
      IndicesTransformation it = allTransformsVect[iTransform];
      IndicesTransformation itInverse = it.inverse();

      IndicesTransformation FwdBack = it.compose(itInverse);
      if (FwdBack != IndicesTransformation::Identity)
        {
          return -1;
        }
      IndicesTransformation BackFwd = itInverse.compose(it);
      if (BackFwd != IndicesTransformation::Identity)
        {
          return -2;
        }
    }
  if (verbose)
    {
      pout() << "All transformations passed inverse test."
             << endl;
    }

  /*
    Test that when each transformation is refined and coarsened,
    the end result is the original transformation.
  */
  int refRatio = 3;
  for (int iTransform = 0; iTransform < numTransformations; iTransform++)
    {
      IndicesTransformation it = allTransformsVect[iTransform];
      IndicesTransformation itRefined = it.refine(refRatio);
      IndicesTransformation itRefinedCoarse = itRefined.coarsen(refRatio);
      if (itRefinedCoarse != it)
        {
          return -3;
        }
    }
  if (verbose)
    {
      pout() << "All transformations passed refinement and coarsening test."
             << endl;
    }

  /*
    Test transformations of boxes with different centerings.
  */
  IntVect srcLo = -2*IntVect::Unit;
  IntVect srcHi = 5*IntVect::Unit;
  Box centeringBox(IntVect::Zero, IntVect::Unit);
  for (BoxIterator bitCen(centeringBox); bitCen.ok(); ++bitCen)
    {
      IntVect centering = bitCen();
      Box srcBox(srcLo, srcHi, centering);

      /*
        Test that Box transformation composed with its inverse
        is identity transformation.
       */
      for (int iTransform = 0; iTransform < numTransformations; iTransform++)
        {
          IndicesTransformation it = allTransformsVect[iTransform];
          Box FwdBackBox = it.transformBack(it.transformFwd(srcBox));
          if (FwdBackBox != srcBox)
            {
              pout() << "Failed composite transformation with centering "
                     << centering
                     << ", forward and then back:" << endl;
              pout() << it << endl;
              pout() << srcBox << " to " << FwdBackBox << endl;
              pout() << "intermediate " << it.transformFwd(srcBox) << endl;
              return -4;
            }
          Box BackFwdBox = it.transformFwd(it.transformBack(srcBox));
          if (BackFwdBox != srcBox)
            {
              pout() << "Failed composite transformation with centering "
                     << centering
                     << ", back and then forward:" << endl;
              pout() << it << endl;
              pout() << srcBox << " to " << BackFwdBox << endl;
              pout() << "intermediate " << it.transformBack(srcBox) << endl;
              return -5;
            }
        }
      if (verbose)
        {
          pout() << "Transformations with centering " << centering
                 << " passed inverse test."
                 << endl;
        }

      BaseFab<int> srcFab;
      Interval srcIntvl(0, 0);
      Interval dstIntvl(0, 0);
      if (dataTransformations)
        {
          /*
            Define the BaseFab<int> to be transformed.
            FIXME: Try testing with more than one component.
          */
          srcFab.define(srcBox, 1);
          for (BoxIterator bit(srcBox); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              IntVect ivOff = iv - srcLo;
              // In 3D, set srcFab(iv, 0) = 10*(10*ivOff[0] + ivOff[1]) + ivOff[2].
              int val = ivOff[0];
              for (int idir = 1; idir < SpaceDim; idir++)
                {
                  val *= 10;
                  val += ivOff[idir];
                }
              srcFab(iv, 0) = val;
            }

          /*
            Test forward transformation of BaseFab.
          */
          for (int iTransform = 0; iTransform < numTransformations; iTransform++)
            {
              IndicesTransformation it = allTransformsVect[iTransform];
              Box dstBox = it.transformFwd(srcBox);
              BaseFab<int> dstFab(dstBox, 1);
              it.transformFwd(dstFab, srcFab, srcBox, dstIntvl, srcIntvl);
              for (BoxIterator bit(srcBox); bit.ok(); ++bit)
                {
                  IntVect ivSrc = bit();
                  IntVect ivDst = it.transformWithType(ivSrc, centering, true);
                  CH_assert(dstBox.contains(ivDst));
                  if (dstFab(ivDst, 0) != srcFab(ivSrc, 0))
                    {
                      return -11;
                    }
                }
            }
          if (verbose)
            {
              pout() << "Forward transformations of BaseFab with centering "
                     << centering
                     << " passed."
                     << endl;
            }

          /*
            Test backward transformation of BaseFab.
          */
          for (int iTransform = 0; iTransform < numTransformations; iTransform++)
            {
              IndicesTransformation it = allTransformsVect[iTransform];
              Box dstBox = it.transformBack(srcBox);
              BaseFab<int> dstFab(dstBox, 1);
              it.transformBack(dstFab, srcFab, srcBox, dstIntvl, srcIntvl);
              for (BoxIterator bit(srcBox); bit.ok(); ++bit)
                {
                  IntVect ivSrc = bit();
                  IntVect ivDst = it.transformWithType(ivSrc, centering, false);
                  CH_assert(dstBox.contains(ivDst));
                  if (dstFab(ivDst, 0) != srcFab(ivSrc, 0))
                    {
                      return -12;
                    }
                }
            }
          if (verbose)
            {
              pout() << "Backward transformations of BaseFab with centering "
                     << centering
                     << " passed."
                     << endl;
            }
        }

      /*
        Now try collapsing srcBox in each node-centered dimension.
       */
      for (int faceDir = 0; faceDir < SpaceDim; faceDir++)
        {
          if (centering == BASISV(faceDir))
            {
              IntVect srcFaceHi = srcHi;
              srcFaceHi[faceDir] = srcLo[faceDir];
              Box srcFaceBox(srcLo, srcFaceHi, centering);

              /*
                Test that Box transformation composed with its inverse
                is identity transformation.
              */
              for (int iTransform = 0; iTransform < numTransformations; iTransform++)
                {
                  IndicesTransformation it = allTransformsVect[iTransform];
                  Box FwdBackBox = it.transformBack(it.transformFwd(srcFaceBox));
                  if (FwdBackBox != srcFaceBox)
                    {
                      return -13;
                    }
                  Box BackFwdBox = it.transformFwd(it.transformBack(srcFaceBox));
                  if (BackFwdBox != srcFaceBox)
                    {
                      return -14;
                    }
                }
              if (verbose)
                {
                  pout() << "Transformations with centering " << centering
                         << " passed inverse test on face-centered box."
                         << endl;
                }

              if (dataTransformations)
                {
                  BaseFab<int> srcFaceFab(srcFaceBox, 1);
                  srcFaceFab.copy(srcFab);

                  /*
                    Test forward transformation of face-centered BaseFab.
                  */
                  for (int iTransform = 0; iTransform < numTransformations; iTransform++)
                    {
                      IndicesTransformation it = allTransformsVect[iTransform];
                      Box dstFaceBox = it.transformFwd(srcFaceBox);
                      BaseFab<int> dstFaceFab(dstFaceBox, 1);
                      it.transformFwd(dstFaceFab, srcFaceFab, srcFaceBox, dstIntvl, srcIntvl);
                      for (BoxIterator bit(srcFaceBox); bit.ok(); ++bit)
                        {
                          IntVect ivSrc = bit();
                          IntVect ivDst = it.transformWithType(ivSrc, centering, true);
                          CH_assert(dstFaceBox.contains(ivDst));
                          if (dstFaceFab(ivDst, 0) != srcFaceFab(ivSrc, 0))
                            {
                              return -15;
                            }
                        }
                    }
                  if (verbose)
                    {
                      pout() << "Forward transformations of face-centered BaseFab with centering "
                             << centering
                             << " passed."
                             << endl;
                    }

                  /*
                    Test backward transformation of face-centered BaseFab.
                  */
                  for (int iTransform = 0; iTransform < numTransformations; iTransform++)
                    {
                      IndicesTransformation it = allTransformsVect[iTransform];
                      Box dstFaceBox = it.transformBack(srcFaceBox);
                      BaseFab<int> dstFaceFab(dstFaceBox, 1);
                      it.transformBack(dstFaceFab, srcFaceFab, srcFaceBox, dstIntvl, srcIntvl);
                      for (BoxIterator bit(srcFaceBox); bit.ok(); ++bit)
                        {
                          IntVect ivSrc = bit();
                          IntVect ivDst = it.transformWithType(ivSrc, centering, false);
                          CH_assert(dstFaceBox.contains(ivDst));
                          if (dstFaceFab(ivDst, 0) != srcFaceFab(ivSrc, 0))
                            {
                              return -16;
                            }
                        }
                    }
                  if (verbose)
                    {
                      pout() << "Backward transformations of face-centered BaseFab with centering "
                             << centering
                             << " passed."
                             << endl;
                    }

                }
            }
        }
    }

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

  if ( dataTransformations )
    {
      pout() << "Including transformations of BaseFab." << endl;
    }
  else
    {
      pout() << "NO transformations of BaseFab; use flag -d to include them." << endl;
    }

  int status = testAllIndicesTransformation();

  if (status == 0)
  {
    pout() << indent << "IndicesTransformation test"
           << " passed." << endl ;
  }
  else
  {
    pout() << indent << "IndicesTransformation test"
           << " failed with return code " << status << endl ;
  }

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize ();
#endif

  return status ;
}
