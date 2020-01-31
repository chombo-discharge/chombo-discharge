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
using std::cout;
using std::endl;
using std::ofstream;


#include "parstream.H"
#include "ParmParse.H"
#include "FourthOrderUtil.H"
#include "BoxIterator.H"
#include "computeNorm.H"
#include "FABView.H"

/// Global variables for handling output:
static const char* pgmname = "testUtil" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = false ;

#ifdef CH_USE_DOUBLE
static Real precision = 1.0e-15;
#else
static Real precision = 1.0e-7;
#endif

enum probtypeEnum
{
  constant = 0,
  linear,
  quadratic,
  sinusoidal,
  num_probtype
};

//int probtype = constant;
//int probtype = linear;
//int probtype = quadratic;
int probtype = sinusoidal;

int phi2Type = sinusoidal;

RealVect oldfaceFerrL1[SpaceDim];
RealVect oldfaceFerrL2[SpaceDim];
RealVect oldfaceFerrMax[SpaceDim];


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

void initFaceData(LevelData<FluxBox>& a_F,
                  const ProblemDomain& a_levelDomain,
                  const RealVect& a_dxLevel)
{
  Real Pi = 4.0*atan(1.0);

  DataIterator dit = a_F.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisF = a_F[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          // set offest (location of face center relative
          // to dx*intVect
          RealVect offset = a_dxLevel;
          offset *= 0.5;
          offset[dir] = 0.0;

          FArrayBox& thisFdir = thisF[dir];

          // this is the really slow way to do this
          BoxIterator bit(thisFdir.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc = iv*a_dxLevel + offset;

              //first shot at things --  set up a constant field
              if (probtype == constant)
                {
                  thisFdir(iv, 0) = 1;
                }
              else if (probtype == linear)
                {
                  thisFdir(iv,0) = loc[0];
                }
              else if (probtype == linear)
                {
                  thisFdir(iv,0) = D_TERM(loc[0],+loc[1],+loc[2]);
                }
              else if (probtype == quadratic)
                {
                  thisFdir(iv,0) = D_TERM(loc[0]*loc[0],
                                          +loc[1]*loc[1],
                                          +loc[2]*loc[2]);
                } // end if quadratic

              else if (probtype == sinusoidal)
                {
                  thisFdir(iv,0) = D_TERM(sin(2*Pi*loc[0]),
                                          +sin(2.*Pi*loc[1]),
                                          +sin(2.*Pi*loc[2]));

                } // end if sinusoidal
              else
                {
                  // bad probtype
                  pout() << "Bad probtype == " << probtype << endl;
                  MayDay::Error("Don't know what to do");
                }

            } // end loop over face cells
        } // end loop over dir
    } // end loop over boxes

}

void initFaceData2(LevelData<FluxBox>& a_F,
                  const ProblemDomain& a_levelDomain,
                  const RealVect& a_dxLevel)
{
  Real Pi = 4.0*atan(1.0);

  DataIterator dit = a_F.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisF = a_F[dit];
      for (int dir=0; dir<SpaceDim; dir++)
        {
          // set offest (location of face center relative
          // to dx*intVect
          RealVect offset = a_dxLevel;
          offset *= 0.5;
          offset[dir] = 0.0;

          FArrayBox& thisFdir = thisF[dir];

          // this is the really slow way to do this
          BoxIterator bit(thisFdir.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc = iv*a_dxLevel + offset;

              //first shot at things --  set up a constant field
              if (probtype == constant)
                {
                  thisFdir(iv, 0) = 1;
                }
              else if (probtype == linear)
                {
                  thisFdir(iv,0) = loc[0];
                }
              else if (probtype == linear)
                {
                  thisFdir(iv,0) = D_TERM(loc[0],+loc[1],+loc[2]);
                }
              else if (probtype == quadratic)
                {
                  thisFdir(iv,0) = D_TERM(loc[0]*loc[0],
                                          +loc[1]*loc[1],
                                          +loc[2]*loc[2]);
                } // end if quadratic

              else if (probtype == sinusoidal)
                {
                  thisFdir(iv,0) = D_TERM(sin(2*Pi*loc[0]),
                                          *sin(2.*Pi*loc[1]),
                                          *sin(2.*Pi*loc[2]));

                } // end if sinusoidal
              else
                {
                  // bad probtype
                  pout() << "Bad probtype == " << probtype << endl;
                  MayDay::Error("Don't know what to do");
                }

            } // end loop over face cells
        } // end loop over dir
    } // end loop over boxes

}


void initCellData(LevelData<FArrayBox>& a_F,
                  const ProblemDomain& a_levelDomain,
                  const RealVect& a_dxLevel)
{
  Real Pi = 4.0*atan(1.0);

  // set offest (location of cell center relative
  // to dx*intVect
  RealVect offset = a_dxLevel;
  offset *= 0.5;

  DataIterator dit = a_F.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisF = a_F[dit];

      // this is the really slow way to do this
      BoxIterator bit(thisF.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc = iv*a_dxLevel + offset;

          //first shot at things --  set up a constant field
          if (probtype == constant)
            {
              thisF(iv, 0) = 1;
            }
          else if (probtype == linear)
            {
              thisF(iv,0) = loc[0];
            }
          else if (probtype == linear)
            {
              thisF(iv,0) = D_TERM(loc[0],+loc[1],+loc[2]);
            }
          else if (probtype == quadratic)
            {
              thisF(iv,0) = D_TERM(loc[0]*loc[0],
                                   +loc[1]*loc[1],
                                   +loc[2]*loc[2]);
            } // end if quadratic

          else if (probtype == sinusoidal)
            {
              // make it positive everywhere
              thisF(iv,0) = D_TERM(sin(2*Pi*loc[0]),
                                   +sin(2.*Pi*loc[1]),
                                   +sin(2.*Pi*loc[2]));

            } // end if sinusoidal
          else
            {
              // bad probtype
              pout() << "Bad probtype == " << probtype << endl;
              MayDay::Error("Don't know what to do");
            }

        } // end loop over  cells
    } // end loop over boxes
}


void initCellData2(LevelData<FArrayBox>& a_F,
                   const ProblemDomain& a_levelDomain,
                   const RealVect& a_dxLevel)
{
  Real Pi = 4.0*atan(1.0);

  // set offest (location of cell center relative
  // to dx*intVect
  RealVect offset = a_dxLevel;
  offset *= 0.5;

  DataIterator dit = a_F.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisF = a_F[dit];

      // this is the really slow way to do this
      BoxIterator bit(thisF.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc = iv*a_dxLevel + offset;

          //first shot at things --  set up a constant field
          if (phi2Type == constant)
            {
              thisF(iv, 0) = 1;
            }
          else if (phi2Type == linear)
            {
              thisF(iv,0) = loc[0];
            }
          else if (phi2Type == linear)
            {
              thisF(iv,0) = D_TERM(loc[0],+loc[1],+loc[2]);
            }
          else if (phi2Type == quadratic)
            {
              thisF(iv,0) = D_TERM(loc[0]*(1.0-loc[0]),
                                   +loc[1]*(1.0-loc[1]),
                                   +loc[2]*(1.0-loc[2]));
            } // end if quadratic

          else if (phi2Type == sinusoidal)
            {
              thisF(iv,0) = 0.5 + SpaceDim + D_TERM(cos(2*Pi*(0.25+loc[0])),
                                                    +cos(2.*Pi*(0.25+loc[1])),
                                                    +cos(2.*Pi*(0.25+loc[2])));

            } // end if sinusoidal
          else
            {
              // bad probtype
              pout() << "Bad phi2Type == " << phi2Type << endl;
              MayDay::Error("Don't know what to do");
            }

        } // end loop over  cells
    } // end loop over boxes
}



void exactCellF(LevelData<FArrayBox>& a_exactF,
                const ProblemDomain& a_levelDomain,
                const RealVect& a_dxLevel)
{
  Real Pi = 4.0*atan(1.0);

  // set offest (location of cell center relative
  // to dx*intVect
  RealVect offset = a_dxLevel;
  offset *= 0.5;

  DataIterator dit = a_exactF.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisF = a_exactF[dit];

      // this is the really slow way to do this
      BoxIterator bit(thisF.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc = iv*a_dxLevel + offset;

          if (probtype == constant)
            {
              thisF(iv, 0) = 0.0;
            }
          else if (probtype == linear)
            {
              thisF(iv,0) = 3.0;
            }
          else if (probtype == quadratic)
            {

            }
          else if (probtype == sinusoidal)
            {
              // first remove offset (reference from lower left corner)
              loc -= offset;
              RealVect hiLoc = loc;
              hiLoc += a_dxLevel;

              loc *= 2*Pi;
              hiLoc *= 2.*Pi;

              Real mult = 1.0/(2.0*Pi*a_dxLevel[0]);

              thisF(iv,0) = mult*(D_TERM((cos(loc[0])-cos(hiLoc[0])),
                                         +(cos(loc[1])-cos(hiLoc[1])),
                                         +(cos(loc[2])-cos(hiLoc[2]))) );

            }
          else
            {
              // bad probtype
              pout() << "Bad probtype == " << probtype << endl;
              MayDay::Error("Don't know what to do");
            }
        }
    } // end loop over boxes
}



void exactFaceF(LevelData<FluxBox>& a_exactF,
                const ProblemDomain& a_levelDomain,
                const RealVect& a_dxLevel)
{
  Real Pi = 4.0*atan(1.0);

  // set offest (location of cell center relative
  // to dx*intVect
  RealVect offset = a_dxLevel;
  offset *= 0.5;

  DataIterator dit = a_exactF.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisF = a_exactF[dit];

      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisFdir = thisF[dir];

          // this is the really slow way to do this
          BoxIterator bit(thisFdir.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc = iv*a_dxLevel + offset;

              if (probtype == constant)
                {
                  thisFdir(iv, 0) = 0.0;
                }
              else if (probtype == linear)
                {
                  thisFdir(iv,0) = 3.0;
                }
              else if (probtype == quadratic)
                {

                }
              else if (probtype == sinusoidal)
                {
                  // first remove offset (reference from lower left corner)
                  loc -= offset;
                  RealVect hiLoc = loc;
                  hiLoc += a_dxLevel;
                  hiLoc[dir] -= a_dxLevel[dir];

                  loc *= 2*Pi;
                  hiLoc *= 2.*Pi;

                  Real mult = 1.0/(2.0*Pi*a_dxLevel[0]);

                  Real val = 0.0;
                  for (int dir1=0; dir1<SpaceDim; dir1++)
                    {
                      if (dir1 != dir)
                        {
                          val += mult*(cos(loc[dir1])-cos(hiLoc[dir1]));
                        }
                      else
                        {
                          val += sin(loc[dir]);
                        }
                    }

                  thisFdir(iv,0) = val;

                }
              else
                {
                  // bad probtype
                  pout() << "Bad probtype == " << probtype << endl;
                  MayDay::Error("Don't know what to do");
                }
            } // loop over cells
        } // loop over face directions
    } // end loop over boxes
}


void exactFaceProduct(LevelData<FluxBox>& a_exactF,
                      const ProblemDomain& a_levelDomain,
                      const RealVect& a_dxLevel)
{
  Real Pi = 4.0*atan(1.0);

  // set offest (location of cell center relative
  // to dx*intVect
  RealVect offset = a_dxLevel;
  offset *= 0.5;

  DataIterator dit = a_exactF.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisF = a_exactF[dit];

      for (int dir=0; dir<SpaceDim; dir++)
        {
          FArrayBox& thisFdir = thisF[dir];

          Real faceArea = 1.0;
          for (int sideDir = 0; sideDir<SpaceDim; sideDir++)
            {
              if (sideDir != dir)
                {
                  faceArea *= a_dxLevel[sideDir];
                }
            }

          // this is the really slow way to do this
          BoxIterator bit(thisFdir.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc = iv*a_dxLevel + offset;

              if (probtype == constant)
                {
                  MayDay::Error("constant probtype not impmlemented for product");
                }
              else if (probtype == linear)
                {
                  MayDay::Error("linear probtype not implemented for product");
                }
              else if (probtype == quadratic)
                {
                  MayDay::Error("linear probtype not implemented for product");
                }
              else if (probtype == sinusoidal)
                {
                  // first remove offset (reference from lower left corner)
                  loc -= offset;
                  RealVect hiLoc = loc;
                  hiLoc += a_dxLevel;
                  hiLoc[dir] -= a_dxLevel[dir];

                  loc *= 2*Pi;
                  hiLoc *= 2.*Pi;

                  Real mult = 1.0/(8.0*Pi);

                  Real val = 1.0;
                  for (int dir1=0; dir1<SpaceDim; dir1++)
                    {
                      if (dir1 != dir)
                        {
                          val *= 0.5*a_dxLevel[dir] - mult*(sin(2.0*hiLoc[dir1])-sin(2.0*loc[dir1]));
                        }
                      else
                        {
                          val *= sin(loc[dir])*sin(loc[dir]);
                        }
                    }

                  thisFdir(iv,0) = val/faceArea;

                }
              else
                {
                  // bad probtype
                  pout() << "Bad probtype == " << probtype << endl;
                  MayDay::Error("Don't know what to do");
                }
            } // loop over cells
        } // loop over face directions
    } // end loop over boxes
}


int
testUtil();

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testUtil();

  if ( status == 0 )
    cout << indent << "Util test" << " passed." << endl ;
  else
    cout << indent << "Util test" << " failed with return code "
         << status << endl ;



  if ( status == 0 )
    cout << indent << pgmname << " passed." << endl ;
  else
    cout << indent << pgmname << " failed with return code " << status << endl ;


#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return status ;
}

int
testUtil()
{
  int returnCode = 0;
  // error tolerance
  Real eps = 1.0e-7;

  int numCells = 8;

  if (probtype == constant)
    {
      pout() << "constant F test problem" << endl;
    }
  else if (probtype == linear)
    {
      pout() << "linear test problem" << endl;
    }
  else if (probtype == quadratic)
    {
      pout() << "quadratic test problem" << endl;
    }
  else if (probtype == sinusoidal)
    {
      pout() << "sinusoidal test problem" << endl;
    }
  else
    {
      pout() << "unknown test problem type" << endl;
    }

  // initialize oldfaceFerr to zero
  for (int dir=0; dir<SpaceDim; dir++)
    {
      oldfaceFerrL1[dir] = RealVect::Zero;
      oldfaceFerrL2[dir] = RealVect::Zero;
      oldfaceFerrMax[dir] = RealVect::Zero;
    }

  IntVect loVect(IntVect::Zero);
  IntVect hiVect( (numCells-1)*IntVect::Unit);

  Box domainBox(loVect, hiVect);
  ProblemDomain baseDomain(domainBox);

  // set to periodic for cases where we don't want to handle
  // boundary conditions
  for (int dir=0; dir<SpaceDim; dir++)
    {
      baseDomain.setPeriodic(dir,true);
    }

  int numLevels = 5;

  Vector<int> nRefVect(numLevels-1, 2);

  Real dx = 1.0/numCells;
  RealVect dxVect(D_DECL(dx, dx, dx));

  // simplest possible test
  ProblemDomain levelDomain = baseDomain;
  Vector<Box> boxes(1, domainBox);
  Vector<int> procAssign(1,0);

  IntVect ghostVect = IntVect::Unit;

  Real L1Crse, L2Crse, maxCrse;
  Real maxFaceErr[SpaceDim];
  Real maxFaceErrCrse[SpaceDim];
  Real maxCellToFaceErrCrse[SpaceDim];
  Real maxFaceProdErrCrse[SpaceDim];
  Real L1DivideErrCrse, L2DivideErrCrse, maxDivideErrCrse;


  RealVect dxLevel = dxVect;
  for (int level = 0; level<numLevels; level++)
    {
      DisjointBoxLayout levelGrids(boxes, procAssign, levelDomain);

      // test cell-centered fourth-order averages
      {
        LevelData<FArrayBox> F(levelGrids, 1, ghostVect);

        // initialize cell-centered point values
        initCellData(F, levelDomain, dxLevel);

        // now compute fourth-order cell-average values
        fourthOrderAverage(F);

        // now evaluate error
        LevelData<FArrayBox> error(levelGrids, 1, IntVect::Zero);


        exactCellF(error, levelDomain, dxLevel);

        DataIterator dit = error.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            FArrayBox& thisF = F[dit];
            FArrayBox& thisErr = error[dit];

            thisErr.minus(thisF,
                          levelGrids[dit],
                          0,0,1);
          }

        DisjointBoxLayout* fineGridPtr = NULL;
        int nRefFine = -1;
        Interval errInterval(0,0);

        Real L1Err = computeNorm(error, fineGridPtr,
                                 nRefFine, dxLevel[0],
                                 errInterval, 1);

        Real L2Err = computeNorm(error, fineGridPtr,
                                 nRefFine, dxLevel[0],
                                 errInterval, 2);

        Real maxErr = computeNorm(error, fineGridPtr,
                                  nRefFine, dxLevel[0],
                                  errInterval, 0);

        pout() << "   L1(err) = " << L1Err
               << "   L2(err) = " << L2Err
               << "   Max(err) = " << maxErr << endl;

        if (level > 0)
          {
            Real L1rate, L2rate, maxRate;
            L1rate = log(Abs(L1Crse/L1Err))/log(2.0);
            L2rate = log(Abs(L2Crse/L2Err))/log(2.0);
            maxRate = log(Abs(maxCrse/maxErr))/log(2.0);
            pout() << "         rate:  L1 = " << L1rate
                   << "   L2 = " << L2rate
                   << "   Max = " << maxRate << endl;

            // check convergence rates for pass/fail
            Real baselineRate = 3.8;
            // it can take a little while to get into
            // the asymptotic regime
            if (level == 1) baselineRate = 3.5;
;
            bool fail = false;
            if ((L1rate < baselineRate) && (L1Err > precision))
              {
                fail = true;
              }
            if ((L2rate < baselineRate) && (L2Err > precision))
              {
                fail = true;
              }
            if ((maxRate < baselineRate) && (maxErr > precision))
              {
                fail = true;
              }

            if (fail)
              {
                returnCode += 1;
                pout() << "CC fourth-order averaging convergence test FAILS"
                       << endl;
              }
          }
        L1Crse = L1Err;
        L2Crse = L2Err;
        maxCrse = maxErr;

      }


      // test face-centered fourth-order averages
      {
        LevelData<FluxBox> F(levelGrids, 1, ghostVect);

        // initialize cell-centered point values
        initFaceData(F, levelDomain, dxLevel);

        // now compute fourth-order cell-average values
        fourthOrderAverage(F);

        // now evaluate error
        LevelData<FluxBox> error(levelGrids, 1, IntVect::Zero);


        exactFaceF(error, levelDomain, dxLevel);

        D_TERM(maxFaceErr[0] = 0.0;,
               maxFaceErr[1] = 0.0;,
               maxFaceErr[2] = 0.0;)

        DataIterator dit = error.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            for (int dir=0; dir<SpaceDim; dir++)
              {
                FArrayBox& thisF = F[dit][dir];
                FArrayBox& thisErr = error[dit][dir];

                thisErr.minus(thisF,
                              thisErr.box(),
                              0,0,1);
                Real thisMaxErr = thisErr.max();
                Real thisMinErr = thisErr.min();

                if (abs(thisMaxErr) > maxFaceErr[dir]) maxFaceErr[dir] = abs(thisMaxErr);
                if (abs(thisMinErr) > maxFaceErr[dir]) maxFaceErr[dir] = abs(thisMinErr);

              }
          }


        pout() << "  Max(faceErr) = "
               << D_TERM(maxFaceErr[0],
                         << ",   " << maxFaceErr[1],
                         << ",   " << maxFaceErr[2])
                         << endl;

        if (level > 0)
          {
            Real rate;
            // also check convergence rates for pass/fail
            Real baselineRate = 3.8;
            // it can take a little while to get into
            // the asymptotic regime
            if (level == 1) baselineRate = 3.5;

            bool fail = false;

            pout() << "         rate: ";

              for (int dir=0; dir<SpaceDim; dir++)
                {
                  rate = log(Abs(maxFaceErrCrse[dir]/maxFaceErr[dir]))/log(2.0);
                  pout() << "   " << rate;

                  if ((rate <baselineRate) && (maxFaceErr[dir] > precision))
                    {
                      fail = true;
                    }
                }

              pout() << endl;



            if (fail)
              {
                returnCode += 10;
                pout() << "FC fourth-order averaging convergence test FAILS"
                       << endl;
              }
          }
        D_TERM(maxFaceErrCrse[0] = maxFaceErr[0];,
               maxFaceErrCrse[1] = maxFaceErr[1];,
               maxFaceErrCrse[2] = maxFaceErr[2];)

      }



      // test cell-to-face fourth-order averaging
      {
        IntVect cellGhost = ghostVect + 2*IntVect::Unit;
        LevelData<FArrayBox> F(levelGrids, 1, cellGhost);
        LevelData<FluxBox> faceF(levelGrids, 1, ghostVect);

        // initialize cell-centered point values
        initCellData(F, levelDomain, dxLevel);

        // now compute fourth-order cell-average values
        fourthOrderAverage(F);

        // fourth-order cell-to-face
        fourthOrderCellToFace(faceF, F);


        // now compute fourth-order face-average values
        //fourthOrderAverage(faceF);

        // now evaluate error
        LevelData<FluxBox> error(levelGrids, 1, IntVect::Zero);


        exactFaceF(error, levelDomain, dxLevel);

        D_TERM(maxFaceErr[0] = 0.0;,
               maxFaceErr[1] = 0.0;,
               maxFaceErr[2] = 0.0;)

        DataIterator dit = error.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            for (int dir=0; dir<SpaceDim; dir++)
              {
                FArrayBox& thisF = faceF[dit][dir];
                FArrayBox& thisErr = error[dit][dir];

                thisErr.minus(thisF,
                              thisErr.box(),
                              0,0,1);
                Real thisMaxErr = thisErr.max();
                Real thisMinErr = thisErr.min();

                if (abs(thisMaxErr) > maxFaceErr[dir]) maxFaceErr[dir] = abs(thisMaxErr);
                if (abs(thisMinErr) > maxFaceErr[dir]) maxFaceErr[dir] = abs(thisMinErr);

              }
          }


        pout() << "  Max(cellToFaceErr) = "
               << D_TERM(maxFaceErr[0],
                         << ",   " << maxFaceErr[1],
                         << ",   " << maxFaceErr[2])
                         << endl;

        if (level > 0)
          {
            Real rate;
            // also check convergence rates for pass/fail
            Real baselineRate = 3.8;
            // it can take a little while to get into
            // the asymptotic regime
            if (level == 1) baselineRate = 3.5;
            bool fail = false;

            pout() << "         rate: ";

              for (int dir=0; dir<SpaceDim; dir++)
                {
                  rate = log(Abs(maxCellToFaceErrCrse[dir]/maxFaceErr[dir]))/log(2.0);
                  pout() << "   " << rate;

                  if ((rate <baselineRate) && (maxFaceErr[dir] > precision))
                    {
                      fail = true;
                    }
                }

              pout() << endl;



            if (fail)
              {
                returnCode += 100;
                pout() << "Cell-To-Face fourth-order averaging convergence test FAILS"
                       << endl;
              }
          }
        D_TERM(maxCellToFaceErrCrse[0] = maxFaceErr[0];,
               maxCellToFaceErrCrse[1] = maxFaceErr[1];,
               maxCellToFaceErrCrse[2] = maxFaceErr[2];)

      }

      // test fourth-order product averages
      {
        LevelData<FluxBox> F(levelGrids, 1, ghostVect);
        LevelData<FluxBox> G(levelGrids, 1, ghostVect);
        LevelData<FluxBox> FtimesG(levelGrids, 1, ghostVect);

        // initialize face-centered point values
        initFaceData2(F, levelDomain, dxLevel);
        // initially try setting both F and G the same
        initFaceData2(G, levelDomain, dxLevel);

        // compute 4th-order averages
        fourthOrderAverage(F);
        fourthOrderAverage(G);

        fourthOrderMultFace(FtimesG, F, G);

        // now compute error
        LevelData<FluxBox> error(levelGrids, 1, IntVect::Zero);

        exactFaceProduct(error, levelDomain, dxLevel);

        D_TERM(maxFaceErr[0] = 0.0;,
               maxFaceErr[1] = 0.0;,
               maxFaceErr[2] = 0.0;)

        DataIterator dit = error.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            for (int dir=0; dir<SpaceDim; dir++)
              {
                FArrayBox& thisF = FtimesG[dit][dir];
                FArrayBox& thisErr = error[dit][dir];

                thisErr.minus(thisF,
                              thisErr.box(),
                              0,0,1);
                Real thisMaxErr = thisErr.max();
                Real thisMinErr = thisErr.min();

                if (abs(thisMaxErr) > maxFaceErr[dir]) maxFaceErr[dir] = abs(thisMaxErr);
                if (abs(thisMinErr) > maxFaceErr[dir]) maxFaceErr[dir] = abs(thisMinErr);

              }
          }


        pout() << "  Max(FaceProductErr) = "
               << D_TERM(maxFaceErr[0],
                         << ",   " << maxFaceErr[1],
                         << ",   " << maxFaceErr[2])
                         << endl;

        if (level > 0)
          {
            Real rate;
            // also check convergence rates for pass/fail
            Real baselineRate = 3.8;
            // it can take a little while to get into
            // the asymptotic regime
            if (level == 1) baselineRate = 3.5;
            bool fail = false;

            pout() << "         rate: ";

              for (int dir=0; dir<SpaceDim; dir++)
                {
                  rate = log(Abs(maxFaceProdErrCrse[dir]/maxFaceErr[dir]))/log(2.0);
                  pout() << "   " << rate;

                  if ((rate <baselineRate) && (maxFaceErr[dir] > precision))
                    {
                      fail = true;
                    }
                }

              pout() << endl;



            if (fail)
              {
                returnCode += 100;
                pout() << "Face-centered product fourth-order averaging convergence test FAILS"
                       << endl;
              }
          }
        D_TERM(maxFaceProdErrCrse[0] = maxFaceErr[0];,
               maxFaceProdErrCrse[1] = maxFaceErr[1];,
               maxFaceProdErrCrse[2] = maxFaceErr[2];)

      }

      // test cell-centered fourth-order division
      {
        LevelData<FArrayBox> F(levelGrids, 1, ghostVect);
        LevelData<FArrayBox> G(levelGrids, 1, ghostVect);
        LevelData<FArrayBox> FG(levelGrids, 1, ghostVect);
        LevelData<FArrayBox> computedF(levelGrids, 1, IntVect::Zero);

        // initialize cell-centered point values
        initCellData(F, levelDomain, dxLevel);
        initCellData2(G, levelDomain, dxLevel);

        // multiply the two together to get F*G point values
        DataIterator dit = levelGrids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            FG[dit].copy(F[dit]);
            FG[dit].mult(G[dit]);
          }

        // now compute fourth-order cell-average values
        fourthOrderAverage(FG);
        fourthOrderAverage(G);

        FG.exchange();
        G.exchange();

        // now divide G, which should give us F back
        for (dit.begin(); dit.ok(); ++dit)
          {
            cellFGToCellF(computedF[dit], FG[dit], G[dit], levelGrids[dit]);
          }

        // now evaluate error
        LevelData<FArrayBox> error(levelGrids, 1, IntVect::Zero);


        exactCellF(error, levelDomain, dxLevel);

        for (dit.begin(); dit.ok(); ++dit)
          {
            FArrayBox& thisF = computedF[dit];
            FArrayBox& thisErr = error[dit];

            thisErr.minus(thisF,
                          levelGrids[dit],
                          0,0,1);
          }

        DisjointBoxLayout* fineGridPtr = NULL;
        int nRefFine = -1;
        Interval errInterval(0,0);

        Real L1Err = computeNorm(error, fineGridPtr,
                                 nRefFine, dxLevel[0],
                                 errInterval, 1);

        Real L2Err = computeNorm(error, fineGridPtr,
                                 nRefFine, dxLevel[0],
                                 errInterval, 2);

        Real maxErr = computeNorm(error, fineGridPtr,
                                  nRefFine, dxLevel[0],
                                  errInterval, 0);

        pout() << "   L1(err(FG/G)) = " << L1Err
               << "   L2(err(FG/G)) = " << L2Err
               << "   Max(err(FG/G)) = " << maxErr << endl;

        if (level > 0)
          {
            Real L1rate, L2rate, maxRate;
            L1rate = log(Abs(L1DivideErrCrse/L1Err))/log(2.0);
            L2rate = log(Abs(L2DivideErrCrse/L2Err))/log(2.0);
            maxRate = log(Abs(maxDivideErrCrse/maxErr))/log(2.0);
            pout() << "         rate:  L1 = " << L1rate
                   << "   L2 = " << L2rate
                   << "   Max = " << maxRate << endl;

            // check convergence rates for pass/fail
            Real baselineRate = 3.8;
            // it can take a little while to get into
            // the asymptotic regime
            if (level == 1) baselineRate = 3.5;
;
            bool fail = false;
            if ((L1rate < baselineRate) && (L1Err > precision))
              {
                fail = true;
              }
            if ((L2rate < baselineRate) && (L2Err > precision))
              {
                fail = true;
              }
            // a bit more relaxed for this maxnorm...
            if ((maxRate < baselineRate-0.25) && (maxErr > precision))
              {
                fail = true;
              }

            if (fail)
              {
                returnCode += 1;
                pout() << "CC division convergence test FAILS"
                       << endl;
              }
          }
        L1DivideErrCrse = L1Err;
        L2DivideErrCrse = L2Err;
        maxDivideErrCrse = maxErr;

      }




      if (level < nRefVect.size())
        {
          levelDomain.refine(nRefVect[level]);
          dxLevel /= nRefVect[level];
          for (int i=0; i<boxes.size(); i++)
            {
              boxes[i].refine(nRefVect[level]);
            }
        }
    } // end loop over levels

  return returnCode;

}



