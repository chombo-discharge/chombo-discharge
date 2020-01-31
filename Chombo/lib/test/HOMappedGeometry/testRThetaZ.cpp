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
#include "newMappedGridIO.H"
using std::cout;
using std::endl;
using std::ofstream;

#include "parstream.H"
#include "ParmParse.H"
#include "RThetaZCS.H"
#include "BRMeshRefine.H"
#include "BoxIterator.H"
#include "computeNorm.H"
#include "FABView.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testRThetaZ" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = false ;
static bool writePlotFiles = true;
static bool writeN11 = false;


#ifdef CH_USE_DOUBLE
static Real precision = 1.0e-15;
#else
static Real precision = 1.0e-7;
#endif

enum probtypeEnum
{
  constantF = 0,
  linearRZ,
  linearRZsinTheta,
  quadraticRZTheta,
  quadraticRZsinTheta,
  num_probtype
};

//int probtype = constantF;
//int probtype = linearRZ;
//int probtype = linearRZsinTheta;
//int probtype = quadraticRZTheta;
int probtype = quadraticRZsinTheta;
//int probtype = cubicF;
//int probtype = quarticF;
//int probtype = quinticF;

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

void initF(LevelData<FluxBox>& a_F,
           const ProblemDomain& a_levelDomain,
           const RealVect& a_dxLevel,
           const NewCoordSys& a_coordsys)
{
  //Real Pi = 4.0*atan(1.0);
  LevelData<FluxBox> Ferror(a_F.getBoxes(), a_F.nComp(),
                            a_F.ghostVect());

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
          FArrayBox& thisFErrDir = Ferror[dit][dir];

          // this is the really slow way to do this
          BoxIterator bit(thisFdir.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc = iv*a_dxLevel + offset;
              RealVect realLoc = a_coordsys.realCoord(loc);
              //first shot at things --  set up a constant field
              if (probtype == constantF)
                {
                  D_TERM(
                         thisFdir(iv, 0) = 1;,
                         thisFdir(iv,1) = 1;,
                         thisFdir(iv,2) = 1;)
                }
              else if (probtype == linearRZ)
                {
                  D_TERM(
                         thisFdir(iv,0) = realLoc[0];,
                         thisFdir(iv,1) = realLoc[1];,
                         thisFdir(iv,2) = realLoc[2];)
                }
              else if (probtype == linearRZsinTheta)
                {
                  Real uTheta = sin(loc[1]);
                  D_TERM(
                         thisFdir(iv,0) = realLoc[0] - sin(loc[1])*uTheta;,
                         thisFdir(iv,1) = realLoc[1] + cos(loc[1])*uTheta;,
                         thisFdir(iv,2) = realLoc[2];)
                }
              else if (probtype == quadraticRZTheta)
                {
                  //Real uTheta = loc[0]*loc[1]*(2*Pi - loc[1]);
                  Real uTheta = 0.0;
                  Real uR = loc[0]*loc[0];
                  //Real uR = 0.0;

                  //bool useFaceAvgVariables = true;
                  //if (useFaceAvgVariables) {
                  // face-average variables
                  if (dir == 0)
                    {
                      // r- faces
                      D_TERM(
                             thisFdir(iv,0) = uR*(sin(loc[1]+0.5*a_dxLevel[1]) - sin(loc[1]-0.5*a_dxLevel[1]))/a_dxLevel[1];,
                             thisFdir(iv,1) = uR*(cos(loc[1]-0.5*a_dxLevel[1]) - cos(loc[1]+0.5*a_dxLevel[1]))/a_dxLevel[1];,
                             thisFdir(iv,2) = loc[2]*loc[2] + a_dxLevel[2]*a_dxLevel[2]/12.0;)
                        }
                  else if (dir == 1)
                    {
                      // theta- faces
                      D_TERM(
                             thisFdir(iv,0) = (uR+a_dxLevel[0]*a_dxLevel[0]/12.0)*cos(loc[1]);,
                             thisFdir(iv,1) = (uR+a_dxLevel[0]*a_dxLevel[0]/12.0)*sin(loc[1]);,
                             thisFdir(iv,2) = loc[2]*loc[2] + a_dxLevel[2]*a_dxLevel[2]/12.0;)
                        }
#if CH_SPACEDIM >= 2
                  else if (dir == 2)
                    {
                      // z-faces
                      // this is simpler from the corner rather than the face center
                      RealVect cornerLoc = a_dxLevel*iv;
                      D_TERM(
                             Real r = cornerLoc[0];,
                             Real theta = cornerLoc[1];,
                             ;
                             //Real z = cornerLoc[2];
                             )

                      D_TERM(
                             Real dr = a_dxLevel[0];,
                             Real dtheta = a_dxLevel[1];,
                             ;
                             //Real dz = a_dxLevel[2];
                             )

                      D_TERM(
                             thisFdir(iv,0) = 0.5*(4.*r*r*r + 6*r*r*dr +4*r*dr*dr + dr*dr*dr)*(sin(theta+dtheta) - sin(theta))/(dtheta*(2*r+dr));,
                             thisFdir(iv,1) = 0.5*(4.*r*r*r + 6*r*r*dr +4*r*dr*dr + dr*dr*dr)*(cos(theta) - cos(theta+dtheta))/(dtheta*(2*r+dr));,
                             thisFdir(iv,2) = loc[2]*loc[2];)
                        }
#endif
                  //}
                  //else
                  //{

                  // compute 4th order approximation to face-averaged variables

                  D_TERM(thisFErrDir(iv,0) = thisFdir(iv,0);,
                         thisFErrDir(iv,1) = thisFdir(iv,1);,
                         thisFErrDir(iv,2) = thisFdir(iv,2);)

                    // point values
                    D_TERM(
                           thisFdir(iv,0) = uR*cos(loc[1]) - sin(loc[1])*uTheta;,
                           thisFdir(iv,1) = uR*sin(loc[1]) + cos(loc[1])*uTheta;,
                           thisFdir(iv,2) = loc[2]*loc[2];)

                    bool fourthOrderApprox = true;

                  if (fourthOrderApprox)
                    {
                      // now increment with corrections
                      for (int tanDir=0; tanDir<SpaceDim; tanDir++)
                        {
                          if (tanDir != dir)
                            {
                              RealVect plusLoc = loc;
                              plusLoc[tanDir] += a_dxLevel[tanDir];
                              RealVect minusLoc = loc;
                              minusLoc[tanDir] -= a_dxLevel[tanDir];
                              Real uRPlus = plusLoc[0]*plusLoc[0];
                              //Real uRPlus = 0.0;
                              //Real uThetaPlus = plusLoc[0]*plusLoc[1]*(2*Pi - plusLoc[1]);
                              Real uThetaPlus = 0.0;
#if CH_SPACEDIM == 3
                              Real uZPlus = plusLoc[2]*plusLoc[2];
#endif
                              Real uRMinus = minusLoc[0]*minusLoc[0];
                              //Real uRMinus = 0.0;
                              //Real uThetaMinus = minusLoc[0]*minusLoc[1]*(2*Pi - minusLoc[1]);
                              Real uThetaMinus = 0.0;
#if CH_SPACEDIM == 3
                              Real uZMinus = minusLoc[2]*minusLoc[2];
#endif
                              RealVect Fplus, Fminus;
                              D_TERM(
                                     Fplus[0] = uRPlus*cos(plusLoc[1]) - sin(plusLoc[1])*uThetaPlus;,
                                     Fplus[1] = uRPlus*sin(plusLoc[1]) + cos(plusLoc[1])*uThetaPlus;,
                                     Fplus[2] = uZPlus;)

                              D_TERM(
                                     Fminus[0] = uRMinus*cos(minusLoc[1]) - sin(minusLoc[1])*uThetaMinus;,
                                     Fminus[1] = uRMinus*sin(minusLoc[1]) + cos(minusLoc[1])*uThetaMinus;,
                                     Fminus[2] = uZMinus;)

                              RealVect secondDerivTerm;
                              D_TERM(// note undivided differences, since we'd just multiply by dx^2 anyway
                                     secondDerivTerm[0] = (Fplus[0] + Fminus[0] - 2.0*thisFdir(iv,0));,
                                     secondDerivTerm[1] = (Fplus[1] + Fminus[1] - 2.0*thisFdir(iv,1));,
                                     secondDerivTerm[2] = (Fplus[2] + Fminus[2] - 2.0*thisFdir(iv,2));)

                              RealVect exactSecondDeriv;
                              if (tanDir == 0)
                                {
                                  D_TERM(
                                         exactSecondDeriv[0] = 2*cos(loc[1]);,
                                         exactSecondDeriv[1] = 2*sin(loc[1]);,
                                         exactSecondDeriv[2] = 0.0;)
                                    }
                              else if (tanDir == 1)
                                {
                                  D_TERM(
                                         exactSecondDeriv[0] = -loc[0]*loc[0]*cos(loc[1]);,
                                         exactSecondDeriv[1] = -loc[0]*loc[0]*sin(loc[1]);,
                                         exactSecondDeriv[2] = 0.0;)
                                    }
                              else if (tanDir == 2)
                                {
                                  D_TERM(
                                         exactSecondDeriv[0] = 0.0;,
                                         exactSecondDeriv[1] = 0.0;,
                                         exactSecondDeriv[2] = 2.0;)
                                    }
                              exactSecondDeriv *= a_dxLevel[tanDir]*a_dxLevel[tanDir];

                              //RealVect secondDerivErr = exactSecondDeriv - secondDerivTerm;
#if 0
                              D_TERM(thisFdir(iv,0) += exactSecondDeriv[0]/24.0;,
                                     thisFdir(iv,1) += exactSecondDeriv[1]/24.0;,
                                     thisFdir(iv,2) += exactSecondDeriv[2]/24.0;)

#endif
                              D_TERM(thisFdir(iv,0) += secondDerivTerm[0]/24.0;,
                                     thisFdir(iv,1) += secondDerivTerm[1]/24.0;,
                                     thisFdir(iv,2) += secondDerivTerm[2]/24.0;)


                            } // end if tangential direction
                        } // end loop over tangential directions
                    } // end if using 4th-order corrections

                    //} // end numerically computing 4th order face values
                    } // end if quadraticRZTheta

              else if (probtype == quadraticRZsinTheta)
                {
                  Real uTheta = loc[0]*sin(loc[1]);
                  //Real uTheta = 0.0;
                  Real uR = loc[0]*loc[0];
                  //Real uR = 0.0;

                  // this is simpler from the corner than from the face center
                  RealVect cornerLoc = a_dxLevel*iv;

                  //bool useFaceAvgVariables = true;
                  //if (useFaceAvgVariables) {
                  // face-average variables
                  if (dir == 0)
                    {
                      // r- faces
                      D_TERM(
                             thisFdir(iv,0) = uR*(sin(loc[1]+0.5*a_dxLevel[1]) - sin(loc[1]-0.5*a_dxLevel[1]))/a_dxLevel[1] - 0.5*loc[0]*(a_dxLevel[1] - 0.5*(sin(2.*(cornerLoc[1]+a_dxLevel[1])) - sin(2.*cornerLoc[1])))/a_dxLevel[1];,
                             thisFdir(iv,1) = uR*(cos(loc[1]-0.5*a_dxLevel[1]) - cos(loc[1]+0.5*a_dxLevel[1]))/a_dxLevel[1] - 0.25*loc[0]*(cos(2.*(cornerLoc[1]+a_dxLevel[1])) - cos(2.*cornerLoc[1]))/(a_dxLevel[1]);,
                             thisFdir(iv,2) = loc[2]*loc[2] + a_dxLevel[2]*a_dxLevel[2]/12.0;)
                        }
                  else if (dir == 1)
                    {
                      // theta- faces
                      D_TERM(
                             thisFdir(iv,0) = (uR+a_dxLevel[0]*a_dxLevel[0]/12.0)*cos(loc[1]) -(loc[0]*sin(loc[1])*sin(loc[1]));,
                             thisFdir(iv,1) = (uR+a_dxLevel[0]*a_dxLevel[0]/12.0)*sin(loc[1]) +(loc[0]*sin(loc[1])*cos(loc[1]));,
                             thisFdir(iv,2) = loc[2]*loc[2] + a_dxLevel[2]*a_dxLevel[2]/12.0;)
                        }
#if CH_SPACEDIM >= 2
                  else if (dir == 2)
                    {
                      // z-faces
                      // this is simpler from the corner rather than the face center
                      //RealVect cornerLoc = a_dxLevel*iv;
                      D_TERM(
                             Real r = cornerLoc[0];,
                             Real theta = cornerLoc[1];,
                             ;
                             //Real z = cornerLoc[2];
                             )

                      D_TERM(
                             Real dr = a_dxLevel[0];,
                             Real dtheta = a_dxLevel[1];,
                             ;
                             //Real dz = a_dxLevel[2];
                             )

                      D_TERM(
                             thisFdir(iv,0) = 0.5*(4.*r*r*r + 6*r*r*dr +4*r*dr*dr + dr*dr*dr)*(sin(theta+dtheta) - sin(theta))/(dtheta*(2*r+dr))
                                             -(r*r +r*dr +dr*dr/3.)*(a_dxLevel[1] - 0.5*(sin(2.*(cornerLoc[1]+a_dxLevel[1])) - sin(2.*cornerLoc[1])))/((2.*r +dr)*dtheta);,
                             thisFdir(iv,1) = 0.5*(4.*r*r*r + 6*r*r*dr +4*r*dr*dr + dr*dr*dr)*(cos(theta) - cos(theta+dtheta))/(dtheta*(2*r+dr))
                                             -0.25*(r*r +r*dr +dr*dr/3.)*(cos(2.*(cornerLoc[1]+a_dxLevel[1])) - cos(2.*cornerLoc[1]))/((r+0.5*dr)*dtheta);,
                             thisFdir(iv,2) = loc[2]*loc[2];)
                        }
#endif
                  //}
                  //else
                  //{

                  // compute 4th order approximation to face-averaged variables

                  D_TERM(thisFErrDir(iv,0) = thisFdir(iv,0);,
                         thisFErrDir(iv,1) = thisFdir(iv,1);,
                         thisFErrDir(iv,2) = thisFdir(iv,2);)

                    // point values
                    D_TERM(
                           thisFdir(iv,0) = uR*cos(loc[1]) - sin(loc[1])*uTheta;,
                           thisFdir(iv,1) = uR*sin(loc[1]) + cos(loc[1])*uTheta;,
                           thisFdir(iv,2) = loc[2]*loc[2];)

                    bool fourthOrderApprox = true;

                  if (fourthOrderApprox)
                    {
                      RealVect lapTerm = RealVect::Zero;
                      // now increment with corrections
                      for (int tanDir=0; tanDir<SpaceDim; tanDir++)
                        {
#if CH_SPACEDIM >= 2
                          if (tanDir != dir)
                            {
                              RealVect plusLoc = loc;
                              plusLoc[tanDir] += a_dxLevel[tanDir];
                              RealVect minusLoc = loc;
                              minusLoc[tanDir] -= a_dxLevel[tanDir];
                              D_TERM(
                                     Real uRPlus = plusLoc[0]*plusLoc[0];
                                     Real uRMinus = minusLoc[0]*minusLoc[0];,
                                     //Real uRMinus = 0.0;
                                     //Real uRPlus = 0.0;
                                     Real uThetaPlus = plusLoc[0]*sin(plusLoc[1]);
                                     Real uThetaMinus = minusLoc[0]*sin(minusLoc[1]);,
                                     //Real uThetaMinus = 0.0;
                                     //Real uThetaPlus = 0.0;
                                     Real uZPlus = plusLoc[2]*plusLoc[2];
                                     Real uZMinus = minusLoc[2]*minusLoc[2];)
                              RealVect Fplus, Fminus;
                              D_TERM(
                                     Fplus[0] = uRPlus*cos(plusLoc[1]) - sin(plusLoc[1])*uThetaPlus;,
                                     Fplus[1] = uRPlus*sin(plusLoc[1]) + cos(plusLoc[1])*uThetaPlus;,
                                     Fplus[2] = uZPlus;)

                              D_TERM(
                                     Fminus[0] = uRMinus*cos(minusLoc[1]) - sin(minusLoc[1])*uThetaMinus;,
                                     Fminus[1] = uRMinus*sin(minusLoc[1]) + cos(minusLoc[1])*uThetaMinus;,
                                     Fminus[2] = uZMinus;)

                              RealVect secondDerivTerm;
                              D_TERM(// note undivided differences, since we'd just multiply by dx^2 anyway
                                     secondDerivTerm[0] = (Fplus[0] + Fminus[0] - 2.0*thisFdir(iv,0));,
                                     secondDerivTerm[1] = (Fplus[1] + Fminus[1] - 2.0*thisFdir(iv,1));,
                                     secondDerivTerm[2] = (Fplus[2] + Fminus[2] - 2.0*thisFdir(iv,2));)

                              RealVect exactSecondDeriv;
                              if (tanDir == 0)
                                {
                                  D_TERM(
                                         exactSecondDeriv[0] = 2*cos(loc[1]);,
                                         exactSecondDeriv[1] = 2*sin(loc[1]);,
                                         exactSecondDeriv[2] = 0.0;)
                                    }
                              else if (tanDir == 1)
                                {
                                  D_TERM(
                                         exactSecondDeriv[0] = -loc[0]*loc[0]*cos(loc[1]) -2.*loc[0]*(cos(loc[1])*cos(loc[1]) - sin(loc[1])*sin(loc[1]));,
                                         exactSecondDeriv[1] = -loc[0]*loc[0]*sin(loc[1])-4.*loc[0]*sin(loc[1])*cos(loc[1]);,
                                         exactSecondDeriv[2] = 0.0;)
                                    }
                              else if (tanDir == 2)
                                {
                                  D_TERM(
                                         exactSecondDeriv[0] = 0.0;,
                                         exactSecondDeriv[1] = 0.0;,
                                         exactSecondDeriv[2] = 2.0;)
                                    }
                              exactSecondDeriv *= a_dxLevel[tanDir]*a_dxLevel[tanDir];

                            //  RealVect secondDerivErr = exactSecondDeriv - secondDerivTerm;

                              // if we're using analytic second-derivatives here
#if 0
                              D_TERM(lapTerm[0] += exactSecondDeriv[0];,
                                     lapTerm[1] += exactSecondDeriv[1];,
                                     lapTerm[2] += exactSecondDeriv[2];)
#endif

                                //#if 0
                              D_TERM(lapTerm[0] += secondDerivTerm[0];,
                                     lapTerm[1] += secondDerivTerm[1];,
                                     lapTerm[2] += secondDerivTerm[2];)
                                //#endif

                            } // end if tangential direction
#endif
                        } // end loop over tangential directions

                      // now add in 4th-order correction
                      D_TERM(thisFdir(iv,0) += lapTerm[0]/24.0;,
                             thisFdir(iv,1) += lapTerm[1]/24.0;,
                             thisFdir(iv,2) += lapTerm[2]/24.0;)

                    } // end if using 4th-order corrections

                } // end if quadraticRZsinTheta
              else
                {
                  // bad probtype
                  pout() << "Bad probtype == " << probtype << endl;
                  MayDay::Error("Don't know what to do");
                }

            } // end loop over face cells
        } // end loop over dir
    } // end loop over boxes

  // convert to 4th order avgs
  //a_coordsys->fourthOrderAverage(a_F);

  // if verbose, compute errors and convergence rates for F
  if (verbose)
    {
      for (dit.begin(); dit.ok(); ++dit)
        {
          for (int dir=0; dir<SpaceDim; dir++)
            {
              FArrayBox& thisFErrDir = Ferror[dit][dir];
              FArrayBox& thisFdir = a_F[dit][dir];

              BoxIterator bit(thisFdir.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();

                  D_TERM(thisFErrDir(iv,0) -= thisFdir(iv,0);,
                         thisFErrDir(iv,1) -= thisFdir(iv,1);,
                         thisFErrDir(iv,2) -= thisFdir(iv,2);)
                    }

              RealVect L1Err, L2Err, MaxErr;
              for (int comp = 0; comp<SpaceDim; comp++)
                {
                  L1Err[comp] = thisFErrDir.norm(1, comp, 1);
                  L2Err[comp] = thisFErrDir.norm(2, comp, 1);
                  MaxErr[comp] = thisFErrDir.norm(0, comp,1);
                }
              Real cellVol =  D_TERM(a_dxLevel[0],*a_dxLevel[1],*a_dxLevel[2]);
              L1Err *= cellVol;
              L2Err *= sqrt(cellVol);
              pout() << "faceDir = " << dir << ": L1(err) = " << L1Err;
              if (oldfaceFerrL1[dir][0] != 0)
                {
                  // compute convergence rate
                  RealVect rate;
                  for (int comp=0; comp<SpaceDim; comp++)
                    {
                      rate[comp] = log(oldfaceFerrL1[dir][comp]/L1Err[comp])/log(2.0);
                    }
                  pout() << ", rate = " << rate;
                }
              pout() << endl;
              pout() << "             L2(err) = " << L2Err;
              if (oldfaceFerrL2[dir][0] != 0)
                {
                  // compute convergence rate
                  RealVect rate;
                  for (int comp=0; comp<SpaceDim; comp++)
                    {
                      rate[comp] = log(oldfaceFerrL2[dir][comp]/L2Err[comp])/log(2.0);
                    }
                  pout() << ", rate = " << rate;
                }

              pout() << endl;
              pout() << "            Max(err) = " << MaxErr;
              if (oldfaceFerrMax[dir][0] != 0)
                {
                  // compute convergence rate
                  RealVect rate;
                  for (int comp=0; comp<SpaceDim; comp++)
                    {
                      rate[comp] = log(oldfaceFerrMax[dir][comp]/MaxErr[comp])/log(2.0);
                    }
                  pout() << ", rate = " << rate;
                }

              pout() << endl;

              oldfaceFerrL1[dir] = L1Err;
              oldfaceFerrL2[dir] = L2Err;
              oldfaceFerrMax[dir] = MaxErr;


            } //end loop over face directions
        } // end loop over boxes
    } // end if verbose for computing errors and convergence rates
}



void exactDivF(LevelData<FArrayBox>& a_divF,
               const ProblemDomain& a_levelDomain,
               const RealVect& a_dxLevel)
{
  Real Pi = 4.0*atan(1.0);

  // set offest (location of cell center relative
  // to dx*intVect
  RealVect offset = a_dxLevel;
  offset *= 0.5;

  //Real cellVol = D_TERM(a_dxLevel[0],*a_dxLevel[1],*a_dxLevel[2]);

  DataIterator dit = a_divF.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& thisDiv = a_divF[dit];

      // this is the really slow way to do this
      BoxIterator bit(thisDiv.box());
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          RealVect loc = iv*a_dxLevel + offset;

          if (probtype == constantF)
            {
              thisDiv(iv, 0) = 0.0;
            }
          else if (probtype == linearRZ)
            {
              thisDiv(iv,0) = 3.0;
            }
          else if (probtype == linearRZsinTheta)
            {
              Real cosTheta = cos(loc[1]);
              thisDiv(iv,0) = 3.0+ cosTheta/loc[0];
            }
          else if (probtype == quadraticRZTheta)
            {
              // use cell volume
              //thisDiv(iv,0) = D_TERM(3*loc[0], +(2*Pi-2*loc[1]), +2*loc[2] );
              // first remove offset (reference from lower left corner)
              loc -= offset;
              thisDiv(iv,0) = D_TERM((3*loc[0]*loc[0]+3.0*loc[0]*a_dxLevel[0] + a_dxLevel[0]*a_dxLevel[0])/(loc[0]+0.5*a_dxLevel[0]),
                                     +(2.0*Pi-2.0*loc[1]-a_dxLevel[1]),
                                     +(2*loc[2] + a_dxLevel[2]));


              thisDiv(iv,0) = D_TERM((3*loc[0]*loc[0]+3.0*loc[0]*a_dxLevel[0] + a_dxLevel[0]*a_dxLevel[0])/(loc[0]+0.5*a_dxLevel[0]),
                                     +0.0,
                                     +(2*loc[2] + a_dxLevel[2]));

              /*
              thisDiv(iv,0) = D_TERM(0.0,
                                     +0.0,
                                     +(2*loc[2] + a_dxLevel[2]));
              */

              //thisDiv(iv,0) = D_TERM(0, +(2*Pi-2*loc[1]), +2*loc[2] );
              //thisDiv(iv,0) = D_TERM(3*loc[0], +0, +2*loc[2] );


            }
          else if (probtype == quadraticRZsinTheta)
            {
              // use cell volume
              //thisDiv(iv,0) = D_TERM(3*loc[0], +cos(loc[1]), +2*loc[2] );
              // first remove offset (reference from lower left corner)
              loc -= offset;
              thisDiv(iv,0) = D_TERM((3*loc[0]*loc[0]+3.0*loc[0]*a_dxLevel[0] + a_dxLevel[0]*a_dxLevel[0])/(loc[0]+0.5*a_dxLevel[0]),
                                     +(sin(loc[1]+a_dxLevel[1]) - sin(loc[1]))/a_dxLevel[1],
                                     +(2*loc[2] + a_dxLevel[2]));


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


void
writeMeshToFile(char* a_filestr,
                const NewCoordSys& a_coordSys,
                const DisjointBoxLayout& a_grids,
                const ProblemDomain& a_domain,
                const RealVect& a_dx)
{

  ofstream os(a_filestr, ios::out);

  if (os.fail())
    {
      pout() << "cannot open grid output file " << a_filestr << endl;
      MayDay::Error();
    }


  DataIterator dit = a_grids.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& gridBox = a_grids[dit];

      IntVect loVect = gridBox.smallEnd();
      IntVect hiVect = gridBox.bigEnd();
      // account for cell centering by adding the Unit vector
      hiVect += IntVect::Unit;

      int ilo = loVect[0];
      int ihi = hiVect[0];
      int jlo = loVect[1];
      int jhi = hiVect[1];

      // r-direction
      IntVect iv(loVect);
      for (int i=ilo; i< ihi; i++)
        {
          iv[0] = i;

          for (int j=jlo; j<=jhi; j++)
            {
              iv[1] = j;

              RealVect loc(iv);
              loc *= a_dx;

              RealVect x = a_coordSys.realCoord(loc);

              os << D_TERM(x[0] << "  ",
                           << x[1] << "  ",
                           << "   ")  << endl;
            } // end loop over j

          // add a blank line between coordinate lines
          os << endl;
        } // end loop over i
      // end x-direction

      // theta-direction
      for (int j=jlo; j<=jhi; j++)
        {
          iv[1] = j;

          for (int i=ilo; i< ihi; i++)
            {
              iv[0] = i;

              RealVect loc(iv);
              loc *= a_dx;

              RealVect x = a_coordSys.realCoord(loc);

              os << D_TERM(x[0] << "  ",
                           << x[1] << "  ",
                           << "   ")  << endl;
            } // end loop over i

          // add a blank line between coordinate lines
          os << endl;
        } // end loop over j
      // end y-direction

    } // end loop over boxes
  os.close();

} // end if we're writing out the mesh



int
testRThetaZ();

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testRThetaZ();

  if ( status == 0 )
    cout << indent << "RThetaZCoordSys test" << " passed." << endl ;
  else
    cout << indent << "RThetaZCoordSys test" << " failed with return code " << status << endl ;



  if ( status == 0 )
    cout << indent << pgmname << " passed." << endl ;
  else
    cout << indent << pgmname << " failed with return code " << status << endl ;


#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize ();
#endif
  return status ;
}

int
testRThetaZ()
{
  int returnCode = 0;
#if CH_USE_DOUBLE
#if CH_SPACEDIM >= 2
  // error tolerance
  Real eps = 1.0e-7;

    int numCells = 8;
  // int numCells = 4;

  if (probtype == constantF)
    {
      pout() << "constant F test problem" << endl;
    }
  else if (probtype == linearRZ)
    {
      pout() << "linear test problem" << endl;
    }
  else if (probtype == linearRZsinTheta)
    {
      pout() << "linear test problem with sin(theta) uTheta " << endl;
    }
  else if (probtype == quadraticRZTheta)
    {
      pout() << "quadratic test problem" << endl;
    }
  else if (probtype == quadraticRZsinTheta)
    {
      pout() << "sin theta, quadratic r,z test problem" << endl;
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

  // want to avoid singularity at the origin, so shift r coordinates
  // outward by, say, 5
  int offset = 5;
  loVect[0] += offset;
  hiVect[0] += offset;
  Box domainBox(loVect, hiVect);
  ProblemDomain baseDomain(domainBox);

  // always periodic in the theta direction
  baseDomain.setPeriodic(1, true);
  //  int numLevels = 5;
  int numLevels = 4;
  if (SpaceDim == 3) numLevels = 3;

  Vector<int> nRefVect(numLevels-1, 2);

  Real dx = 1.0/numCells;
  RealVect dxVect(D_DECL(dx, dx, dx));
  // "dtheta" should be 2Pi/numCells
  Real Pi = 4.0*atan(1.0);
  dxVect[1] = 2.0*Pi/numCells;

  RealVect stretch = RealVect::Unit;

  // set interval for tests
  Interval volInterval(0,SpaceDim-1);
  //Interval volInterval(0,0);
  pout() << "using interval of (" << volInterval.begin()
         << "," << volInterval.end()
         << ") for volume computation" << endl;


  // simplest possible test
  int maxBoxSize = 16;
  ProblemDomain levelDomain = baseDomain;

  // set numbers of ghost cells
  int baseNGhost = 1;
  int numNghost = baseNGhost+3;
  int numCellVolGhost = baseNGhost + 2;
  int numJinvGhost = baseNGhost +1;
  int NinvJGhost = baseNGhost;

  IntVect ghostVect = IntVect::Unit;

  Real L1Crse, L2Crse, maxCrse;
  Real L1volCrse, L2volCrse,maxVolCrse;

  Real maxNErr[SpaceDim][SpaceDim][SpaceDim];
  Real maxNErrCrse[SpaceDim][SpaceDim][SpaceDim];

  RealVect MaxJinverseCrse;
  RealVect MaxNJinverseCrse[SpaceDim][SpaceDim];

  RealVect dxLevel = dxVect;
  for (int level = 0; level<numLevels; level++)
    {

      //Vector<Box> boxes(1, domainBox);
      Vector<Box> boxes;
      domainSplit(domainBox, boxes, maxBoxSize);
      Vector<int> procAssign(boxes.size(),0);


      DisjointBoxLayout levelGrids(boxes, procAssign, levelDomain);

      RThetaZCS levelCoordSys(dxVect, stretch);
      int numNcomp = levelCoordSys.getNumN();
      LevelData<FluxBox> metricTerms(levelGrids, numNcomp,
                                     numNghost*ghostVect);


      {
        // if desired, write mesh to a file
        bool writeMesh = false;
        if (writeMesh)
          {
            char filestr[80];
            RealVect plainDx = dxLevel;
            sprintf(filestr, "grid%d.out", level);
            writeMeshToFile(filestr, levelCoordSys, levelGrids,
                            levelDomain, plainDx);
          }


        // invertibility test
        {
          Real maxErr = 0.0;
          DataIterator dit = levelGrids.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              BoxIterator bit(levelGrids[dit]);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  RealVect Xi = iv*dxLevel;
                  RealVect x = levelCoordSys.realCoord(Xi);
                  RealVect mapped_x = levelCoordSys.mappedCoord(x);
                  RealVect error = Xi - mapped_x;
                  error *= error;
                  Real magErr = sqrt(error.sum());
                  if (magErr > maxErr)
                    {
                      maxErr = magErr;
                    }
                } // end loop over points in this box
            } // end loop over grids

          if (maxErr > eps)
            {
              // fail
              returnCode += 10;
              pout() << "Invertibility test FAILS"
                     << " -- Max(error) = " << maxErr << endl;
            }
        } // end invertibility test


        // check face metric terms
        {

          LevelData<FluxBox> metricErrors(levelGrids, metricTerms.nComp(),
                                          metricTerms.ghostVect());

          // initialize errors to zero
          for (int faceDir=0; faceDir<SpaceDim; faceDir++)
            {
              for (int dir=0; dir<SpaceDim; dir++)
                  {
                    D_TERM(maxNErr[faceDir][dir][0] = 0.0;,
                           maxNErr[faceDir][dir][1] = 0.0;,
                           maxNErr[faceDir][dir][2] = 0.0;)
                      }
            }


          DataIterator dit = metricTerms.dataIterator();

          // compute N and print N_11
          {

            cout << "levelDomain = " << levelDomain << endl;

            int compN11 = levelCoordSys.getNcomponent(0,0);
            for (dit.begin(); dit.ok(); ++dit)
              {
                FluxBox& thisN = metricTerms[dit];
                levelCoordSys.getN(thisN, thisN.box());

                if (writeN11)
                  {

                    for (int dir=0; dir<SpaceDim; dir++)
                      {
                        const FArrayBox & thisNDir = thisN[dir];

                        cout << "N(1,1) in direction " << dir << " is:"
                             << endl;

                        BoxIterator bit = thisNDir.box();
                        for (bit.begin(); bit.ok(); ++bit)
                          {
                            IntVect iv = bit();

                            cout << iv[0] << " " << iv[1] << " " << thisNDir(iv,compN11) << endl;
                          }
                      }
                  } // end if we're writing N11
              }
          }

          for (dit.begin(); dit.ok(); ++dit)
            {
              const FluxBox& thisMetricTerms = metricTerms[dit];
              FluxBox& thisError = metricErrors[dit];
              thisError.setVal(0.0);

              for (int faceDir=0; faceDir<SpaceDim; faceDir++)
                {
                  const FArrayBox& thisMetricDir = thisMetricTerms[faceDir];
                  FArrayBox& thisErrorDir = thisError[faceDir];
                  thisErrorDir.setVal(0.0);

                  // slow, but hopefully not too bad...
                  BoxIterator bit = thisErrorDir.box();
                  for (bit.begin(); bit.ok(); ++bit)
                    {
                      IntVect iv = bit();
                      RealVect rThetaZ = iv*dxLevel;

                      if (SpaceDim == 2)
                        {
                          // r-faces
                          if (faceDir == 0)
                            {
                              // s = 0, d = 0
                              thisErrorDir(iv,0) = rThetaZ[0]*(sin(rThetaZ[1]+dxLevel[1]) - sin(rThetaZ[1]));
                              // s = 0, d = 1
                              thisErrorDir(iv,1) = -rThetaZ[0]*(cos(rThetaZ[1]+dxLevel[1]) - cos(rThetaZ[1]));
                              // s = 1, d = 0
                              thisErrorDir(iv,2) = (cos(rThetaZ[1]+dxLevel[1]) - cos(rThetaZ[1]));
                              //thisErrorDir(iv,2) = -sin(rThetaZ[1] + 0.5*dxLevel[1]);
                              // s = 1, d = 1
                              thisErrorDir(iv,3) = (sin(rThetaZ[1]+dxLevel[1]) - sin(rThetaZ[1]));
                              //thisErrorDir(iv,3) = cos(rThetaZ[1] + 0.5*dxLevel[1]);
                            }
                          else if (faceDir == 1)
                            {
                              thisErrorDir(iv,0) = 0.5*dxLevel[0]*(2.0*rThetaZ[0]+dxLevel[0])*cos(rThetaZ[1]);
                              //thisErrorDir(iv,0) = (rThetaZ[0]+0.5*dxLevel[0])*cos(rThetaZ[1]);
                              thisErrorDir(iv,1) = 0.5*dxLevel[0]*(2.0*rThetaZ[0]+dxLevel[0])*sin(rThetaZ[1]);
                              //thisErrorDir(iv,1) = (rThetaZ[0]+0.5*dxLevel[0])*sin(rThetaZ[1]);
                              thisErrorDir(iv,2) = -dxLevel[0]*sin(rThetaZ[1]);
                              thisErrorDir(iv,3) = dxLevel[0]*cos(rThetaZ[1]);
                            }
                        } // end if 2d
                      else if (SpaceDim == 3)
                        {
                          if (faceDir == 0)
                            {
                              thisErrorDir(iv,0) = dxLevel[2]*rThetaZ[0]*(sin(rThetaZ[1]+dxLevel[1]) - sin(rThetaZ[1]));
                              thisErrorDir(iv,1) = -dxLevel[2]*rThetaZ[0]*(cos(rThetaZ[1]+dxLevel[1]) - cos(rThetaZ[1]));
                              thisErrorDir(iv,2) = 0.0;

                              thisErrorDir(iv,3) = dxLevel[2]*(cos(rThetaZ[1]+dxLevel[1]) - cos(rThetaZ[1]));
                              thisErrorDir(iv,4) = dxLevel[2]*(sin(rThetaZ[1]+dxLevel[1]) - sin(rThetaZ[1]));
                              thisErrorDir(iv,5) = 0.0;

                              thisErrorDir(iv,6) = 0.0;
                              thisErrorDir(iv,7) = 0.0;
                              thisErrorDir(iv,8) = dxLevel[1]*dxLevel[2]*rThetaZ[0];
                            }
                          else if (faceDir == 1)
                            {
                              thisErrorDir(iv,0)  =0.5*dxLevel[0]*(2.0*rThetaZ[0]+dxLevel[0])*dxLevel[2]*cos(rThetaZ[1]);
                              thisErrorDir(iv,1) = 0.5*dxLevel[0]*(2.0*rThetaZ[0]+dxLevel[0])*dxLevel[2]*sin(rThetaZ[1]);
                              thisErrorDir(iv,2) = 0.0;

                              thisErrorDir(iv,3) = -dxLevel[0]*dxLevel[2]*sin(rThetaZ[1]);
                              thisErrorDir(iv,4) = dxLevel[0]*dxLevel[2]*cos(rThetaZ[1]);
                              thisErrorDir(iv,5) = 0.0;

                              thisErrorDir(iv,6) = 0.0;
                              thisErrorDir(iv,7) = 0.0;
                              thisErrorDir(iv,8) = 0.5*dxLevel[0]*(2.0*rThetaZ[0]+dxLevel[0])*dxLevel[2];
                            }
                          else if (faceDir == 2)
                            {
                              thisErrorDir(iv,0)  =0.5*dxLevel[0]*(2.0*rThetaZ[0]+dxLevel[0])*(sin(rThetaZ[1] + dxLevel[1]) -sin(rThetaZ[1]));
                              thisErrorDir(iv,1) = -0.5*dxLevel[0]*(2.0*rThetaZ[0]+dxLevel[0])*(cos(rThetaZ[1] + dxLevel[1]) -cos(rThetaZ[1]));
                              thisErrorDir(iv,2) = 0.0;


                              thisErrorDir(iv,3) = dxLevel[0]*(cos(rThetaZ[1] + dxLevel[1]) -cos(rThetaZ[1]));
                              thisErrorDir(iv,4) = dxLevel[0]*(sin(rThetaZ[1] + dxLevel[1]) -sin(rThetaZ[1]));
                              thisErrorDir(iv,5) = 0.0;

                              thisErrorDir(iv,6) = 0.0;
                              thisErrorDir(iv,7) = 0.0;
                              thisErrorDir(iv,8) = dxLevel[0]*dxLevel[1]*(2*rThetaZ[0] + dxLevel[0])/2;
                            }
                        } // end if 3D
                    } // end loop over points

                  thisErrorDir -= thisMetricDir;

                  if (verbose)
                    pout() << "FaceDir = " << faceDir << endl;

                  // for now, just check normal comps
                  for (int dir = 0; dir<SpaceDim; dir++)
                    {
                      for (int comp=0; comp<SpaceDim; comp++)
                        {
                          int Ncomp = levelCoordSys.getNcomponent(comp,dir);
                          Real thisMin = thisErrorDir.min(Ncomp);
                          Real thisMax = thisErrorDir.max(Ncomp);
                          if (verbose)
                            {
                              pout () << "comp " << comp
                                      << " thisMax = " << thisMax
                                      << endl;

                              pout () << "comp " << comp
                                      << " thisMin = " << thisMin
                                      << endl;
                            }

                          if (abs(thisMin) > maxNErr[faceDir][dir][comp]) maxNErr[faceDir][dir][comp] = abs(thisMin);
                          if (abs(thisMax) > maxNErr[faceDir][dir][comp]) maxNErr[faceDir][dir][comp] = abs(thisMax);
                        } // end loop over components (s)
                    } // end loop over direction (d)
                } // end loop over face directions
            } // end loop over grids

          bool fail = false;
          for (int faceDir =0; faceDir<SpaceDim; faceDir++)
            {
              pout() << "FaceDir = " << faceDir << ":" << endl;

              for (int dDir=0; dDir<SpaceDim; dDir++)
                {
                  Real rate[SpaceDim];

                  pout () << "level " << level << ":";
                  pout()  << "   Max(err(FaceMetric[" << dDir
                          << "] = " << D_TERM(maxNErr[faceDir][dDir][0],
                                              << "    " << maxNErr[faceDir][dDir][1],
                                              << "    " << maxNErr[faceDir][dDir][2]) << endl;
                  if (level > 0)
                    {
                      D_TERM(rate[0]=log(Abs(maxNErrCrse[faceDir][dDir][0]/maxNErr[faceDir][dDir][0]))/log(2.0);,
                             rate[1]=log(Abs(maxNErrCrse[faceDir][dDir][1]/maxNErr[faceDir][dDir][1]))/log(2.0);,
                             rate[2]=log(Abs(maxNErrCrse[faceDir][dDir][2]/maxNErr[faceDir][dDir][2]))/log(2.0);)

                        pout() << "                            rate = "
                               << D_TERM(rate[0],
                                         << "    " << rate[1],
                                         << "    " << rate[2])
                               << endl;

                      Real baselineRate = 3.75 + SpaceDim-1;

                      for (int n=0; n<SpaceDim; n++)
                        {
                          if ((rate[n] < baselineRate)  && (maxNErr[faceDir][dDir][n] > precision))
                            fail = true;
                        }
                    } // end if level > 0

                  D_TERM(maxNErrCrse[faceDir][dDir][0] = maxNErr[faceDir][dDir][0];,
                         maxNErrCrse[faceDir][dDir][1] = maxNErr[faceDir][dDir][1];,
                         maxNErrCrse[faceDir][dDir][2] = maxNErr[faceDir][dDir][2];)

                    } // end loop over d directions on this face
            } // end loop over face directions

          if (fail)
            {
              // fail
              returnCode += 100;
              pout() << "Face metric terms test FAILS" << endl;
            }


        } // end test of face metric terms


        // test cell volumes
        {
          LevelData<FArrayBox> cellVol(levelGrids, 1,
                                       numCellVolGhost*IntVect::Unit);
          LevelData<FArrayBox> cellVolErr(levelGrids, 1,
                                          numCellVolGhost*IntVect::Unit);

          Real baseVol = 0.5*D_TERM(dxLevel[0],*dxLevel[1],*dxLevel[2]);
          DataIterator dit=cellVol.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              FArrayBox& thisCellVol = cellVol[dit];
              levelCoordSys.cellVol(thisCellVol, metricTerms[dit],
                                    thisCellVol.box());

              FArrayBox& thisCellVolErr = cellVolErr[dit];


              BoxIterator bit(thisCellVolErr.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  Real rad = iv[0]*dxLevel[0];
                  Real volume = (2.0*rad + dxLevel[0])*baseVol;
                  thisCellVolErr(iv,0) = volume;
                }
              cellVolErr[dit]-=thisCellVol;

            }

          DisjointBoxLayout* fineGridPtr = NULL;
          int nRefFine = -1;
          Interval divInterval(0,0);

          Real L1Err = computeNorm(cellVolErr, fineGridPtr,
                                   nRefFine, dxLevel[0],
                                   divInterval, 1);

          Real L2Err = computeNorm(cellVolErr, fineGridPtr,
                                   nRefFine, dxLevel[0],
                                   divInterval, 2);

          Real maxErr = computeNorm(cellVolErr, fineGridPtr,
                                    nRefFine, dxLevel[0],
                                    divInterval, 0);

          pout () << "level " << level << ":" << endl;
          pout()  << "   err(cellVol): L1 = " << L1Err
                  << "   L2 = " << L2Err
                  << "   Max = " << maxErr << endl;

          if (level > 0)
            {
              Real L1Volrate, L2Volrate, maxVolRate;
              L1Volrate = log(Abs(L1volCrse/L1Err))/log(2.0);
              L2Volrate = log(Abs(L2volCrse/L2Err))/log(2.0);
              maxVolRate = log(Abs(maxVolCrse/maxErr))/log(2.0);
              pout() << "         rate:  L1 = " << L1Volrate
                     << "   L2 = " << L2Volrate
                     << "   Max = " << maxVolRate << endl;

              // test convergence rate, and fail if it isn't good enough
              bool fail = false;
              Real baselineRate = 3.75 + SpaceDim;

              if ((L1Volrate < baselineRate) && (L1Err > precision) )
                {
                  fail = true;
                }
              if ((L2Volrate < baselineRate) && (L2Err > precision) )
                {
                  fail = true;
                }
              if ((maxVolRate < baselineRate) && (maxErr > precision) )
                {
                  fail = true;
                }

              if (fail)
                {
                  pout() << "Cell volume convergence test FAILS!" << endl;
                  returnCode += 1000;
                }
            }
          L1volCrse = L1Err;
          L2volCrse = L2Err;
          maxVolCrse = maxErr;


        }

        {
          // check face-centered 1/J and N/J terms
          LevelData<FluxBox> Jinverse(levelGrids, 1,
                                      numJinvGhost*IntVect::Unit);

          LevelData<FluxBox> NJinverse(levelGrids, metricTerms.nComp(),
                                       NinvJGhost*IntVect::Unit);

          LevelData<FluxBox> JInverseErrors(levelGrids, Jinverse.nComp(),
                                            Jinverse.ghostVect());

          LevelData<FluxBox> NJInverseErrors(levelGrids, NJinverse.nComp(),
                                             NJinverse.ghostVect());


          RealVect maxJinverseErr(RealVect::Zero);
          RealVect maxErrNJinverse[SpaceDim][SpaceDim];
          for (int dir=0; dir<SpaceDim; dir++)
            {
              D_TERM(
                     maxErrNJinverse[dir][0] = RealVect::Zero;,
                     maxErrNJinverse[dir][1] = RealVect::Zero;,
                     maxErrNJinverse[dir][2] = RealVect::Zero;)
                }

          DataIterator dit = Jinverse.dataIterator();
          for (dit.begin(); dit.ok(); ++dit)
            {
              FluxBox& thisJinverse = Jinverse[dit];
              levelCoordSys.getAvgJinverse(thisJinverse,
                                           thisJinverse.box());

              FluxBox& thisError = JInverseErrors[dit];
              thisError.setVal(0.0);

              FluxBox& thisNJinverse = NJinverse[dit];
              const FluxBox& thisMetricTerms = metricTerms[dit];
              levelCoordSys.computeNJinverse(thisNJinverse,
                                              thisJinverse,
                                              thisMetricTerms,
                                              thisNJinverse.box());

              FluxBox& thisNJerr = NJInverseErrors[dit];

              for (int faceDir=0; faceDir<SpaceDim; faceDir++)
                {
                  const FArrayBox& thisJinverseDir = thisJinverse[faceDir];

                  FArrayBox& thisErrorDir = thisError[faceDir];
                  // exact value should be integral of 1/r over the face
                  BoxIterator bit(thisErrorDir.box());
                  for (bit.begin(); bit.ok(); ++bit)
                    {
                      IntVect iv = bit();
                      // this is r at the lower corner in cell-centered
                      // directions
                      Real r = iv[0]*dxLevel[0];
                      Real exactVal = 1.0/r;
                      if (faceDir != 0)
                          {
                            exactVal = log((r+dxLevel[0])/r)/dxLevel[0];
                          }

                      thisErrorDir(iv, 0) = exactVal;
                    } // end loop over cells to compute exact solution
                  thisErrorDir -= thisJinverseDir;

                  if (verbose)
                    pout() << "Testing 1/J, N/J: FaceDir = " << faceDir << endl;

                  // don't expect correct values at highest and lowest
                  // faces in FaceDir, so don't count those in the norms
                  Box normBox(thisErrorDir.box());
                  normBox.grow(faceDir, -1);

                  Real thisMin = thisErrorDir.min(normBox);
                  Real thisMax = thisErrorDir.max(normBox);
                  if (verbose)
                    {
                      pout () << "faceDir " << faceDir
                              << " thisMax(1/J) = " << thisMax
                              << endl;

                      pout () << "faceDir " << faceDir
                              << " thisMin(1/J) = " << thisMin
                              << endl;
                    }

                  if (abs(thisMin) > maxJinverseErr[faceDir])
                    maxJinverseErr[faceDir] = abs(thisMin);
                  if (abs(thisMax) > maxJinverseErr[faceDir])
                    maxJinverseErr[faceDir] = abs(thisMax);


                  const FArrayBox& thisNJinverseDir = thisNJinverse[faceDir];
                  const FArrayBox& thisNdir = thisMetricTerms[faceDir];
                  FArrayBox& thisNJerrDir = thisNJerr[faceDir];

                  Box NJintersectBox(thisNdir.box());
                  NJintersectBox &= thisNJinverseDir.box();

                  // now compute errors for N/J
                  thisNJerrDir.copy(thisNdir, NJintersectBox);

                  BoxIterator bit2(thisNJerrDir.box());
                  for (bit2.begin(); bit2.ok(); ++bit2)
                    {
                      IntVect iv = bit2();
                      // this is r at the lower corner in cell-centered
                      // directions
                      RealVect rThetaZ = iv*dxLevel;
                      //Real exactJinverse = 1.0/rThetaZ[0];

                      if (SpaceDim == 2)
                        {
                          // r-faces
                          if (faceDir == 0)
                            {
                              // s = 0, d = 0
                              thisNJerrDir(iv,0) = sin(rThetaZ[1]+dxLevel[1]) - sin(rThetaZ[1]);
                              // s = 1, d = 0
                              thisNJerrDir(iv,1) = -1*(cos(rThetaZ[1]+dxLevel[1]) - cos(rThetaZ[1]));
                              // s = 0, d = 1
                              thisNJerrDir(iv,2) = (cos(rThetaZ[1]+dxLevel[1]) - cos(rThetaZ[1]))/rThetaZ[0];
                              // s = 1, d = 1
                              thisNJerrDir(iv,3) = (sin(rThetaZ[1]+dxLevel[1]) - sin(rThetaZ[1]))/rThetaZ[0];
                            } // end if r-face
                          else if (faceDir == 1)
                            {
                              // s = 0, d = 0
                              thisNJerrDir(iv,0) =  dxLevel[0]*cos(rThetaZ[1]);
                              // s = 1, d = 0
                              thisNJerrDir(iv,1) = dxLevel[0]*sin(rThetaZ[1]);
                              // s = 0, d = 1
                              thisNJerrDir(iv,2) = -(log(rThetaZ[0] + dxLevel[0])-log(rThetaZ[0]))*sin(rThetaZ[1]);
                              // s = 1, d = 1
                              thisNJerrDir(iv,3) =  (log(rThetaZ[0] + dxLevel[0])-log(rThetaZ[0]))*cos(rThetaZ[1]);
                            } // end if theta-face
                        } // end if 2d
                      else if (SpaceDim == 3)
                        {
                          if (faceDir == 0)
                            {
                              thisNJerrDir(iv,0) = dxLevel[2]*(sin(rThetaZ[1]+dxLevel[1]) - sin(rThetaZ[1]));
                              thisNJerrDir(iv,1) = -dxLevel[2]*(cos(rThetaZ[1]+dxLevel[1]) - cos(rThetaZ[1]));
                              thisNJerrDir(iv,2) = 0.0;

                              thisNJerrDir(iv,3) = dxLevel[2]*(cos(rThetaZ[1]+dxLevel[1]) - cos(rThetaZ[1]))/rThetaZ[0];
                              thisNJerrDir(iv,4) = dxLevel[2]*(sin(rThetaZ[1]+dxLevel[1]) - sin(rThetaZ[1]))/rThetaZ[0];
                              thisNJerrDir(iv,5) = 0.0;

                              thisNJerrDir(iv,6) = 0.0;
                              thisNJerrDir(iv,7) = 0.0;
                              thisNJerrDir(iv,8) = dxLevel[1]*dxLevel[2];
                            }
                          else if (faceDir == 1)
                            {
                              thisNJerrDir(iv,0)  =dxLevel[0]*dxLevel[2]*cos(rThetaZ[1]);
                              thisNJerrDir(iv,1) = dxLevel[0]*dxLevel[2]*sin(rThetaZ[1]);
                              thisNJerrDir(iv,2) = 0.0;

                              thisNJerrDir(iv,3) = -(log(rThetaZ[0]+dxLevel[0])-log(rThetaZ[0]))*dxLevel[2]*sin(rThetaZ[1]);
                              thisNJerrDir(iv,4) =  (log(rThetaZ[0]+dxLevel[0])-log(rThetaZ[0]))*dxLevel[2]*cos(rThetaZ[1]);
                              thisNJerrDir(iv,5) = 0.0;

                              thisNJerrDir(iv,6) = 0.0;
                              thisNJerrDir(iv,7) = 0.0;
                              thisNJerrDir(iv,8) = dxLevel[0]*dxLevel[2];
                            }
                          else if (faceDir == 2)
                            {
                              thisNJerrDir(iv,0) = dxLevel[0]*(sin(rThetaZ[1] + dxLevel[1]) -sin(rThetaZ[1]));
                              thisNJerrDir(iv,1) =-dxLevel[0]*(cos(rThetaZ[1] + dxLevel[1]) -cos(rThetaZ[1]));
                              thisNJerrDir(iv,2) = 0.0;


                              thisNJerrDir(iv,3) = (log(rThetaZ[0]+dxLevel[0])-log(rThetaZ[0]))*(cos(rThetaZ[1] + dxLevel[1]) -cos(rThetaZ[1]));
                              thisNJerrDir(iv,4) = (log(rThetaZ[0]+dxLevel[0])-log(rThetaZ[0]))*(sin(rThetaZ[1] + dxLevel[1]) -sin(rThetaZ[1]));
                              thisNJerrDir(iv,5) = 0.0;

                              thisNJerrDir(iv,6) = 0.0;
                              thisNJerrDir(iv,7) = 0.0;
                              thisNJerrDir(iv,8) = dxLevel[0]*dxLevel[1];
                            }
                        }

                    } // end loop over faces to compute exact NJinverse

                  thisNJerrDir -= thisNJinverseDir;

                  //for (int comp=0; comp<thisNJinverseDir.nComp(); comp++)
                  for (int sComp=0; sComp<SpaceDim; sComp++)
                    {
                      for (int dDir=0; dDir<SpaceDim; dDir++)
                        {
                          int thisComp = levelCoordSys.getNcomponent(sComp,
                                                                     dDir);
                          Real thisMin = thisNJerrDir.min(thisComp);
                          Real thisMax = thisNJerrDir.max(thisComp);
                          if (verbose)
                            {
                              pout () << "sComp, dDir = " << sComp << ", " << dDir
                                      << " thisMax(N/J) = " << thisMax
                                      << endl;

                              pout () << "sComp, dDir = " << sComp << ", " << dDir
                                      << " thisMin(N/J) = " << thisMin
                                      << endl;
                            }

                          if (abs(thisMin) > maxErrNJinverse[faceDir][dDir][sComp])
                            maxErrNJinverse[faceDir][dDir][sComp] = abs(thisMin);
                          if (abs(thisMax) > maxErrNJinverse[faceDir][dDir][sComp])
                            maxErrNJinverse[faceDir][dDir][sComp] = abs(thisMax);
                        } // end loop over dDir
                    } // end loop over sComp
                } // end loop over face directions
            } // end loop over grids

          pout() << endl;
          pout() << "   max(JinverseError) = " << maxJinverseErr;
          if (level > 0)
            {
              RealVect rate;
              bool fail = false;
              Real baselineRate = 3.75;
              for (int dir=0; dir<SpaceDim; dir++)
                {
                  rate[dir]=log(Abs(MaxJinverseCrse[dir]/maxJinverseErr[dir]))/log(2.0);
                  if ((rate[dir] < baselineRate) && (maxJinverseErr[dir] > precision))
                    {
                      fail = true;
                    }
                }
              pout() << ",  rate = " << rate;
              if (fail)
                {
                  returnCode += 1000000;
                  pout() << endl;
                  pout() << "Jinverse convergence test FAILS!" << endl;
                }
            }
          MaxJinverseCrse = maxJinverseErr;

          // now check rate for N/J
          for (int faceDir=0; faceDir<SpaceDim; faceDir++)
            {
              pout() << endl;
              pout() << "FaceDir = " << faceDir << ":" << endl;
              for (int dDir =0; dDir<SpaceDim; ++dDir)
                {
                  pout() << "   max(N/J Error)[" << dDir <<
                    "] = " << maxErrNJinverse[faceDir][dDir];
                  if (level > 0)
                    {
                      bool fail = false;
                      Real baselineRate = 3.75;
                      if ((level <  3  ) && (faceDir == 2) )
                        {
                          // 3d z-direction seems to take a bit longer to get into the
                          // asymptotic regime...
                          baselineRate = 3.35;
                        }
                      RealVect rate;
                      for (int comp=0; comp<SpaceDim; comp++)
                        {
                          rate[comp]=
                            log(Abs(MaxNJinverseCrse[faceDir][dDir][comp]/maxErrNJinverse[faceDir][dDir][comp]))/log(2.0);
                          if ((rate[comp] < baselineRate) && (maxErrNJinverse[faceDir][dDir][comp] > precision))
                            {
                              fail = true;
                              returnCode += 10000;
                              pout() << endl;
                              pout() << "NJinverse test FAILS"
                                     << " -- Max(error) for NJinverse terms = " ;
                              for (int j=0; j< SpaceDim; j++)
                                {
                                  pout() << maxErrNJinverse[faceDir][j] << "   ";
                                }
                              pout() << endl;
                            }

                        }
                      pout() << ",  rate = " << rate;
                    } // end if level > 0
                  pout() << endl;
                } // end loop over dDir

            } // end loop over faceDirections

          for (int faceDir=0; faceDir<SpaceDim; faceDir++)
            {
              for (int dDir = 0; dDir<SpaceDim; dDir++)
                {
                  MaxNJinverseCrse[faceDir][dDir] = maxErrNJinverse[faceDir][dDir];
                }
            }
        }   // end test of Jinverse

        pout() << endl;
        pout() << endl;

        // now test divergence operator
        LevelData<FluxBox> F(levelGrids, SpaceDim, ghostVect);

        initF(F, levelDomain, dxLevel, levelCoordSys);

        LevelData<FArrayBox> divF(levelGrids, 1, IntVect::Zero);

        DataIterator dit = divF.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            Interval divInterval(0, divF.nComp()-1);
            levelCoordSys.computeDivergence(divF[dit], F[dit],
                                            metricTerms[dit],
                                            levelGrids[dit],
                                            divInterval);
          }

        // now evaluate error
        {
          LevelData<FArrayBox> divError(levelGrids, 1);
          exactDivF(divError, levelDomain, dxLevel);

          Real baseVol = 0.5*D_TERM(dxLevel[0],*dxLevel[1],*dxLevel[2]);

          for (dit.begin(); dit.ok(); ++dit)
            {
              FArrayBox& thisDiv = divF[dit];
              BoxIterator bit(thisDiv.box());
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  Real rad = iv[0]*dxLevel[0];
                  Real volume = (2.0*rad + dxLevel[0])*baseVol;
                  thisDiv(iv,0) /=volume;
                }
              divError[dit].minus(divF[dit]);
            }

          DisjointBoxLayout* fineGridPtr = NULL;
          int nRefFine = -1;
          Interval divInterval(0,0);

          Real L1Err = computeNorm(divError, fineGridPtr,
                                   nRefFine, dxLevel[0],
                                   divInterval, 1);

          Real L2Err = computeNorm(divError, fineGridPtr,
                                   nRefFine, dxLevel[0],
                                   divInterval, 2);

          Real maxErr = computeNorm(divError, fineGridPtr,
                                    nRefFine, dxLevel[0],
                                    divInterval, 0);

          pout() << "   L1(err(div) = " << L1Err
                 << "   L2(err(div) = " << L2Err
                 << "   Max(err(div) = " << maxErr << endl;

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
              Real baselineRate = 3.6;
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
                  returnCode += 100000;
                  pout() << "divergence convergence rate test FAILS" << endl;
                }
            }
          L1Crse = L1Err;
          L2Crse = L2Err;
          maxCrse = maxErr;

          // this is a good place to write a test plot file
          // with the divergence and the error in the divergence
          if (writePlotFiles)
            {
#ifdef CH_USE_HDF5
              LevelData<FArrayBox> plotData(levelGrids, 2);
              divF.copyTo(divF.interval(),
                          plotData, divF.interval());
              Interval errorInterval(1,1);
              divError.copyTo(divError.interval(),
                              plotData, errorInterval);

              string fileRoot("RThetaZData.");
              // use unigrid version, since we've only got one level here
                WriteMappedUGHDF5(fileRoot, levelGrids,
                                  plotData, levelCoordSys,
                                  levelDomain.domainBox(),0);
#endif
              } // end if we're writing plot files
        }

        if (level < nRefVect.size())
          {
            domainBox.refine(nRefVect[level]);
            levelDomain.refine(nRefVect[level]);
            dxLevel /= nRefVect[level];
            dxVect /= nRefVect[level];
          }
      } // end if we returnd a valid CoordSys
    } // end loop over levels
#endif
#endif
  return returnCode;
}
