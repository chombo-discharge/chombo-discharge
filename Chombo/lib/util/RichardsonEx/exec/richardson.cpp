#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::ifstream;
using std::ios;

// #define HALEM_PROC_SPEED
#ifdef HALEM_PROC_SPEED
#include <cstdio>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#endif

#ifdef CH_MPI
#include "CH_Attach.H"
#endif

#include "FABView.H"
#include <fstream>

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "memusage.H"
#include "Vector.H"

#include "DisjointBoxLayout.H"
#include "BoxLayout.H"
#include "LevelData.H"
#include "PatchOpIF_F.H"
#include "DebugOut.H"
#include "CoarseAverage.H"
#include "computeNorm.H"
#include "LevelOP.H"
#include "VesselOP.H"

// Debugging stuff
//#include <fenv.h>

#include "UsingNamespace.H"

#ifdef CH_Linux
// Should be undefined by default
#define TRAP_FPE
#undef  TRAP_FPE
#endif

int main(int a_argc, char* a_argv[])
{
//feenableexcept(FE_INVALID | FE_OVERFLOW);
#ifdef CH_MPI
  // Start MPI
  MPI_Init(&a_argc,&a_argv);
#ifdef CH_AIX
  H5dont_atexit();
#endif
  // setChomboMPIErrorHandler();
#endif

  int rank, number_procs;

#ifdef CH_MPI
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  MPI_Comm_size(Chombo_MPI::comm, &number_procs);
  pout() << "Error: Only set up for run in serial" << endl;
  return 0;
#else
  rank = 0;
  number_procs = 1;
#endif

/// Only set up for serial running at the moment

  char* inFile=NULL;

  if (a_argc > 1)
    {
      inFile = a_argv[1];
    }
  else
    {
      pout() << "No input file specified" << endl;
      return -1;
    }

  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

  ParmParse ppre("re");


// This is the refinement between different solution grids
  int refRatio = 0;
  ppre.get("ref_ratio",refRatio);

// Expected order of accuracy
  int orderAccuracy = 1;
  ppre.get("order_accuracy",orderAccuracy);

  int compt = 0;
  ppre.get("comp",compt);

  int verbosity = 0;
  ppre.get("verbosity",verbosity);

  bool useAMR = false;
  ppre.get("amr_bool",useAMR);

  int numLevels = 0;
  ppre.query("num_levels",numLevels);

  std::string singlefileFolder;
  ppre.query("file_folder_single",singlefileFolder);
  singlefileFolder = singlefileFolder + "/";

  std::string amrfileFolder;
  if(useAMR)
    {
      ppre.query("file_folder_amr",amrfileFolder);
      CH_assert(numLevels > 0);
      amrfileFolder = amrfileFolder + "/";
    }

  // Set the ouput directory to the singlefileFolder by default
  std::string outputFile;
  outputFile = singlefileFolder + "/output";
  ppre.query("output_file",outputFile);
  outputFile = outputFile + "/";

  bool latexOutput;
  ppre.query("latex_output",latexOutput);
  std::string latexFileName;
  if(latexOutput)
    {
      ppre.get("latex_file",latexFileName);
    }

  LevelOP postfineLE(compt,orderAccuracy,verbosity);
  LevelOP coarseLE(compt,orderAccuracy,verbosity);
  LevelOP mediumLE(compt,orderAccuracy,verbosity);
  LevelOP fineLE(compt,orderAccuracy,verbosity);
  LevelOP precoarseLE(compt,orderAccuracy,verbosity);

  postfineLE.setSolRefRate(refRatio);
  fineLE.setSolRefRate(refRatio);
  mediumLE.setSolRefRate(refRatio*2.);
  coarseLE.setSolRefRate(refRatio*4.);
  precoarseLE.setSolRefRate(refRatio*8.);

  bool precoarse = false;
  // Create a character string that contains the component number
  char varch[21];
  sprintf(varch, "%d", compt);

  {
// Error check for the defines    
    int ier = 0;
// The first set of solutions do not use AMR
    bool wAMR = false;
//****************************************//
// Collect data from the coarse solution
    std::string coarseFile;
    ppre.query("coarse_file_single",coarseFile);

    coarseFile = singlefileFolder + coarseFile;

    ier += coarseLE.define(coarseFile,
                           0,
                           wAMR);
  
//****************************************//
// Collect data from the medium solution
    std::string mediumFile;
    ppre.query("medium_file_single",mediumFile);

    mediumFile = singlefileFolder + mediumFile;

    ier += mediumLE.define(mediumFile,
                           0,
                           wAMR);

//****************************************//
// Collect data from the fine solution
    std::string fineFile;
    ppre.query("fine_file_single",fineFile);

    fineFile = singlefileFolder + fineFile;

    ier += fineLE.define(fineFile,
                         0,
                         wAMR);

//****************************************//
// Collect data from the post-fine solution
    std::string postfineFile;
    ppre.query("postfine_file_single",postfineFile);

    postfineFile = singlefileFolder + postfineFile;
    
    ier += postfineLE.define(postfineFile,
                             0,
                             wAMR);
    
//****************************************//
// Collect data from the post-fine solution
    std::string precoarseFile;
    ppre.query("precoarse_file_single",precoarseFile);
    if(!precoarseFile.empty())
      {
        precoarseFile = singlefileFolder + precoarseFile;
        precoarseLE.define(precoarseFile,
                           0,
                           wAMR);
        precoarse = true;
      }
    if(ier > 0)
      {
        pout() << "Error reading data from files" << endl;
        return 1;
      }
  }
      
  VesselOP coarseVessel;
  VesselOP mediumVessel;
  VesselOP fineVessel;

  if(useAMR)
    {
// These solution do use AMR
//****************************************//
// Collect data from the coarse and refined coarse solution
      std::string coarseFile;
      ppre.query("coarse_file_amr",coarseFile);

      coarseFile = amrfileFolder + coarseFile;

      coarseVessel.define(coarseFile,
                          numLevels,
                          refRatio,
                          2,
                          compt,
                          orderAccuracy,
                          verbosity);

      std::string mediumFile;
      ppre.query("medium_file_amr",mediumFile);

      mediumFile = amrfileFolder + mediumFile;

      mediumVessel.define(mediumFile,
                          numLevels,
                          refRatio,
                          1,
                          compt,
                          orderAccuracy,
                          verbosity);

      std::string fineFile;
      ppre.query("fine_file_amr",fineFile);

      fineFile = amrfileFolder + fineFile;

      fineVessel.define(fineFile,
                        numLevels,
                        refRatio,
                        0,
                        compt,
                        orderAccuracy,
                        verbosity);
    }

  if(verbosity > 3)
    {
      pout() << "Coarse dx: " << coarseLE.gridDX() << endl;
      pout() << "Medium dx: " << mediumLE.gridDX() << endl;
      pout() << "Fine dx: " << fineLE.gridDX() << endl;
      pout() << "Post fine dx: " << postfineLE.gridDX() << endl;
    }
  

// Define the interpolated level data in each LE object

  coarseLE.defineCoarseAverage(postfineLE.DBL(),
                               postfineLE.LD());
  mediumLE.defineCoarseAverage(postfineLE.DBL(),
                               postfineLE.LD());
  fineLE.defineCoarseAverage(postfineLE.DBL(),
                             postfineLE.LD());
  if(precoarse)
    {
      precoarseLE.defineCoarseAverage(postfineLE.DBL(),
                                      postfineLE.LD());
    }

// Use Richardson Extrapolation (RE) to find a solution that is one order
// of magnitude higher than the order of accuracy used in the solution

// These solutions will be stored in the coarser of the two meshes used
// to find the "analytical" solution

  coarseLE.solRE();
  mediumLE.solRE();
  fineLE.solRE();
  if(precoarse)
    {
      precoarseLE.solRE();
    }


  /// Open up the file to place the convergence rates in
  std::ofstream convrate;

  std::string convrateFile = "convratevar";

  convrateFile = outputFile + convrateFile + varch + ".dat";

  char* conv_cstr = new char[convrateFile.length()+1];

  std::strcpy(conv_cstr,convrateFile.c_str());

  convrate.open(conv_cstr);

  delete[] conv_cstr;

  int precoarsenumcells = 0;
  if(precoarse)
    {
      precoarsenumcells = precoarseLE.boxdomLD().size(0);
    }
  int coarsenumcells = coarseLE.boxdomLD().size(0);
  int mediumnumcells = mediumLE.boxdomLD().size(0);
  int finenumcells = fineLE.boxdomLD().size(0);

  Real totalarea = 1.;
  for(int dir = 0; dir != SpaceDim; ++dir)
    {
      totalarea *= coarseLE.gridDX(dir)*coarseLE.boxdomLD().size(dir);
    }

  {
    if(verbosity > 3)
      {
        pout() << "Creating DAT files" << endl;
      }

    {
      /// Solve for the error on the single grids
      coarseLE.differenceSolutionSolve();
      mediumLE.differenceSolutionSolve();
      fineLE.differenceSolutionSolve();
      if(precoarse)
        {
          precoarseLE.differenceSolutionSolve();
        }
      
      char file_extension[50];
      sprintf(file_extension, "%dd.hdf5", SpaceDim);

        
      /// Write the HDF5 files of the error for the single grids
      {
        std::string difffileName;
        
        difffileName = "singlecoarsediff.";

        difffileName = outputFile + difffileName + file_extension;

        coarseLE.writeData(difffileName);
      }
      {
        std::string difffileName;
        
        difffileName = "singlemediumdiff.";

        difffileName = outputFile + difffileName + file_extension;

        mediumLE.writeData(difffileName);
      }
      {
        std::string difffileName;
        
        difffileName = "singlefinediff.";

        difffileName = outputFile + difffileName + file_extension;

        fineLE.writeData(difffileName);
      }
      if(precoarse)
        {
          std::string difffileName;
        
          difffileName = "singleprecoarsediff.";

          difffileName = outputFile + difffileName + file_extension;

          precoarseLE.writeData(difffileName);
        }

      /// Write the error files for the single grids
      std::ofstream slerrorDat;

      std::string slerrorFile = "slerrorvar";
      
      slerrorFile = outputFile + slerrorFile + varch + ".dat";
        
      char* sl_cstr = new char[slerrorFile.length()+1];
        
      std::strcpy(sl_cstr,slerrorFile.c_str());
        
      slerrorDat.open(sl_cstr);

      delete[] sl_cstr;

      slerrorDat << "N Linfty L1 L2 " << endl;

      Interval UInterval(compt,compt);

      Vector<LevelData<FArrayBox>* > phi(1);

      Vector<int> refVec(1);
      refVec[0] = 1;
      // lnorms is a vector that lists the norms in order of
      // 0: precoarse L inf norm
      // 1: precoarse L 1 norm
      // 2: precoarse L 2 norm
      // 3: coarse L inf norm
      // etc
      Vector<Real> lnorms;
      int numNorms = 3;
      int numCol = 2;
      if(precoarse)
        {
          numCol++;
        }
      if(precoarse)
        {
          phi[0] = &precoarseLE.diffLD();
          Real inDX0 = precoarseLE.gridDX(0);
          slerrorDat << precoarsenumcells;
          for(int i = 0; i != 3; ++i)
            {
              Real normVal = computeNorm(phi,
                                         refVec,
                                         inDX0,
                                         UInterval,
                                         i);
              if(i != 0)
                {
                  normVal *= pow(totalarea,-1./i);
                }
              lnorms.push_back(normVal);
              slerrorDat << " " << lnorms.back();
            }
          slerrorDat << endl;
        }

      phi[0] = &coarseLE.diffLD();
      Real inDX1 = coarseLE.gridDX(0);
      slerrorDat << coarsenumcells;
      for(int i = 0; i != 3; ++i)
        {
          Real normVal = computeNorm(phi,
                                     refVec,
                                     inDX1,
                                     UInterval,
                                     i);
          if(i != 0)
            {
              normVal *= pow(totalarea,-1./i);
            }
          lnorms.push_back(normVal);
          slerrorDat << " " << lnorms.back();
        }
      slerrorDat << endl;

      phi[0] = &mediumLE.diffLD();
      Real inDX2 = mediumLE.gridDX(0);
      slerrorDat << mediumnumcells;
      for(int i = 0; i != 3; ++i)
        {
          Real normVal = computeNorm(phi,
                                     refVec,
                                     inDX2,
                                     UInterval,
                                     i);
          if(i != 0)
            {
              normVal *= pow(totalarea,-1./i);
            }
          lnorms.push_back(normVal);
          slerrorDat << " " << lnorms.back();
        }
      slerrorDat << endl;

      phi[0] = &fineLE.diffLD();
      Real inDX3 = fineLE.gridDX(0);
      slerrorDat << finenumcells;
      for(int i = 0; i != 3; ++i)
        {
          Real normVal = computeNorm(phi,
                                     refVec,
                                     inDX3,
                                     UInterval,
                                     i);
          if(i != 0)
            {
              normVal *= pow(totalarea,-1./i);
            }
          lnorms.push_back(normVal);
          slerrorDat << " " << lnorms.back();
        }
      slerrorDat << endl;

      slerrorDat.close();
      convrate << "Variable: " << compt << endl;
      pout() << "Variable: " << compt << endl;
      convrate << "Single Level: " << endl;
      pout() << "Single Level: " << endl;
      for(int i = 0; i != numNorms; ++i)
        {
          pout() << "L" << i << " ";
          convrate << "L" << i << " ";
          for(int col = 0; col != numCol; ++col)
            {
              const int indx1 = i + numNorms*col;
              const int indx2 = i + numNorms*(col + 1);
              const Real LV1 = lnorms[indx1];
              const Real LV2 = lnorms[indx2];
              const Real Nval = log(LV1/LV2)/(log(refRatio));
              pout() << Nval << " ";
              convrate << Nval << " ";
            }
          pout() << endl;
          convrate << endl;
        }
      if(latexOutput)
        {
          Vector<string> normNames(3);
          normNames[0] = "$L_{\\infty}$";
          normNames[1] = "$L_1$";
          normNames[2] = "$L_2$";
          /// Write latex formatted data
          std::ofstream latexDat;

          std::string latexForm;
      
          latexForm = outputFile + latexFileName + varch + ".tex";
        
          char* sl_cstr = new char[latexForm.length()+1];
        
          std::strcpy(sl_cstr,latexForm.c_str());
        
          latexDat.open(sl_cstr);

          delete[] sl_cstr;
          latexDat.precision(3);
          latexDat << std::scientific;
          int normIndx = 0;
          latexDat << " & $L_{\\infty}$ & " << lnorms[normIndx];
          for(int col = 0; col != numCol; ++col)
            {
              latexDat << " & ";
              const int indx1 = normIndx + numNorms*col;
              const int indx2 = normIndx + numNorms*(col + 1);
              const Real LV1 = lnorms[indx1];
              const Real LV2 = lnorms[indx2];
              const Real Nval = log(LV1/LV2)/(log(refRatio));
              latexDat << "\\textcolor{blue}{" << std::fixed << Nval << "} & "
                       << std::scientific << LV2;
            }
          latexDat << " \\\\" << endl;
          normIndx = 1;
          latexDat << " & $L_1$ & " << lnorms[normIndx];
          for(int col = 0; col != numCol; ++col)
            {
              latexDat << " & ";
              const int indx1 = normIndx + numNorms*col;
              const int indx2 = normIndx + numNorms*(col + 1);
              const Real LV1 = lnorms[indx1];
              const Real LV2 = lnorms[indx2];
              const Real Nval = log(LV1/LV2)/(log(refRatio));
              latexDat << "\\textcolor{blue}{" << std::fixed << Nval << "} & "
                       << std::scientific << LV2;
            }
          latexDat << " \\\\" << endl;
          normIndx = 2;
          latexDat << " & $L_2$ & " << lnorms[normIndx];
          for(int col = 0; col != numCol; ++col)
            {
              latexDat << " & ";
              const int indx1 = normIndx + numNorms*col;
              const int indx2 = normIndx + numNorms*(col + 1);
              const Real LV1 = lnorms[indx1];
              const Real LV2 = lnorms[indx2];
              const Real Nval = log(LV1/LV2)/(log(refRatio));
              latexDat << "\\textcolor{blue}{" << std::fixed << Nval << "} & "
                       << std::scientific << LV2;
            }
          latexDat << " \\\\" << endl;
          latexDat.close();
        }
    }

    /// Create the error HDF5 and dat files for the AMR grids
    if(useAMR)
      {
        int numNorms = 3;
        int numCol = 2;
        Interval UInterval(compt,compt);

        /// Calculate the error for the AMR grids

        coarseVessel.calcError(fineLE.DBL(),
                               fineLE.solLD());
        
        mediumVessel.calcError(fineLE.DBL(),
                               fineLE.solLD());
        
        fineVessel.calcError(fineLE.DBL(),
                             fineLE.solLD());
        
        /// Write the HDF5 files for AMR grids
        char file_extension[50];
        sprintf(file_extension, "%dd.hdf5", SpaceDim);

        {
          std::string difffileName;
        
          difffileName = "amrcoarsediff.";

          difffileName = outputFile + difffileName + file_extension;

          coarseVessel.writeData(difffileName);
        }


        {
          std::string difffileName;
        
          difffileName = "amrmediumdiff.";

          difffileName = outputFile + difffileName + file_extension;
          mediumVessel.writeData(difffileName);
        }
        {
          std::string difffileName;
        
          difffileName = "amrfinediff.";

          difffileName = outputFile + difffileName + file_extension;

          fineVessel.writeData(difffileName);
        }


        /// Compute the norms

        coarseVessel.ComputeNorms();

        mediumVessel.ComputeNorms();

        fineVessel.ComputeNorms();

        /// Output the norms and slopes

        std::ofstream amrerrorDat;

        std::string amrerrorFile = "amrerrorvar";
        
        amrerrorFile = outputFile + amrerrorFile + varch + ".dat";
        
        char* amr_cstr = new char[amrerrorFile.length()+1];
        
        std::strcpy(amr_cstr,amrerrorFile.c_str());
        
        amrerrorDat.open(amr_cstr);

        delete[] amr_cstr;
        Vector<Real> amrlnorms;

        amrerrorDat << "N Linfty L1 L2 " << endl;
        
        Vector<Real> coarseNorms = coarseVessel.returnNorms();
        for(int i = 0; i != 3; ++i)
          {
            amrlnorms.push_back(coarseNorms[i]);
            if(i != 0)
              {
                amrlnorms.back() *= pow(totalarea,-1./i);
              }
            amrerrorDat << " " << amrlnorms.back();
          }
        amrerrorDat << endl;
        
        Vector<Real> mediumNorms = mediumVessel.returnNorms();
        for(int i = 0; i != 3; ++i)
          {
            amrlnorms.push_back(mediumNorms[i]);
            if(i != 0)
              {
                amrlnorms.back() *= pow(totalarea,-1./i);
              }
            amrerrorDat << " " << amrlnorms.back();
          }
        amrerrorDat << endl;
        
        Vector<Real> fineNorms = fineVessel.returnNorms();
        for(int i = 0; i != 3; ++i)
          {
            amrlnorms.push_back(fineNorms[i]);
            if(i != 0)
              {
                amrlnorms.back() *= pow(totalarea,-1./i);
              }
            amrerrorDat << " " << amrlnorms.back();
          }
        amrerrorDat << endl;

        amrerrorDat.close();
        {
          pout() << "Variable: " << compt << endl;
          pout() << "Two Level: " << endl;

          for(int i = 0; i != numNorms; ++i)
            {
              pout() << "L" << i << " ";
              convrate << "L" << i << " ";
              for(int col = 0; col != numCol; ++col)
                {
                  const int indx1 = i + numNorms*col;
                  const int indx2 = i + numNorms*(col + 1);
                  const Real LV1 = amrlnorms[indx1];
                  const Real LV2 = amrlnorms[indx2];
                  const Real Nval = log(LV1/LV2)/(log(refRatio));
                  pout() << Nval << " ";
                  convrate << Nval << " ";
                }
              pout() << endl;
              convrate << endl;
            }
          if(latexOutput)
            {
              Vector<string> normNames(3);
              normNames[0] = "$L_{\\infty}$";
              normNames[1] = "$L_1$";
              normNames[2] = "$L_2$";
              /// Write latex formatted data
              std::ofstream latexDat;

              std::string latexForm;
      
              latexForm = outputFile + latexFileName + varch + "amr.tex";
        
              char* sl_cstr = new char[latexForm.length()+1];
        
              std::strcpy(sl_cstr,latexForm.c_str());
        
              latexDat.open(sl_cstr);

              delete[] sl_cstr;
              latexDat.precision(3);
              latexDat << std::scientific;
              int normIndx = 0;
              latexDat << " & $L_{\\infty}$ & " << amrlnorms[normIndx];
              for(int col = 0; col != numCol; ++col)
                {
                  latexDat << " & ";
                  const int indx1 = normIndx + numNorms*col;
                  const int indx2 = normIndx + numNorms*(col + 1);
                  const Real LV1 = amrlnorms[indx1];
                  const Real LV2 = amrlnorms[indx2];
                  const Real Nval = log(LV1/LV2)/(log(refRatio));
                  latexDat << "\\textcolor{blue}{" << std::fixed
                           << Nval << "} & " << std::scientific << LV2;
                }
              latexDat << " \\\\" << endl;
              normIndx = 1;
              latexDat << " & $L_1$ & " << amrlnorms[normIndx];
              for(int col = 0; col != numCol; ++col)
                {
                  latexDat << " & ";
                  const int indx1 = normIndx + numNorms*col;
                  const int indx2 = normIndx + numNorms*(col + 1);
                  const Real LV1 = amrlnorms[indx1];
                  const Real LV2 = amrlnorms[indx2];
                  const Real Nval = log(LV1/LV2)/(log(refRatio));
                  latexDat << "\\textcolor{blue}{" << std::fixed
                           << Nval << "} & " << std::scientific << LV2;
                }
              latexDat << " \\\\" << endl;
              normIndx = 2;
              latexDat << " & $L_2$ & " << amrlnorms[normIndx];
              for(int col = 0; col != numCol; ++col)
                {
                  latexDat << " & ";
                  const int indx1 = normIndx + numNorms*col;
                  const int indx2 = normIndx + numNorms*(col + 1);
                  const Real LV1 = amrlnorms[indx1];
                  const Real LV2 = amrlnorms[indx2];
                  const Real Nval = log(LV1/LV2)/(log(refRatio));
                  latexDat << "\\textcolor{blue}{" << std::fixed
                           << Nval << "} & " << std::scientific << LV2;
                }
              latexDat << " \\\\" << endl;
              latexDat.close();
            }
        }
      }
  }
  convrate.close();
  return 0;
#ifdef CH_MPI
  // Exit MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}


