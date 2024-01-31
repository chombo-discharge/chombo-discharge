/* chombo-discharge
 * Copyright © 2024 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DischargeIOImplem.H
  @brief  Implementation of CD_DischargeIO.H
  @author Robert Marskar
*/

#ifndef CD_DischargeIOImplem_H
#define CD_DischargeIOImplem_H

// Std includes
#include <hdf5.h>

// Our includes
#include <CD_DischargeIO.H>
#include <CD_NamespaceHeader.H>

#ifdef CH_USE_HDF5
template <size_t M, size_t N>
void
DischargeIO::writeH5Part(const std::string                               a_filename,
                         const ParticleContainer<GenericParticle<M, N>>& a_particles,
                         const std::vector<std::string>                  a_realVars,
                         const std::vector<std::string>                  a_vectVars,
                         const RealVect                                  a_shift) noexcept
{
  CH_TIME("DischargeIO::writeH5Part");

  CH_assert(a_realVars.size() == 0 || a_realVars.size() == M);
  CH_assert(a_vectVars.size() == 0 || a_vecVars.size() == N);

  std::vector<std::string> realVariables(M);
  std::vector<std::string> vectVariables(N);

  if (a_realVars.size() == M) {
    realVariables = a_realVars;

    for (int i = 0; i < M; i++) {
      if (realVariables[i] == "") {
        realVariables[i] = "real-" + std::to_string(i);
      }
    }
  }
  else {
    for (int i = 0; i < M; i++) {
      realVariables[i] = "real-" + std::to_string(i);
    }
  }

  if (a_vectVars.size() == N) {
    vectVariables = a_vectVars;

    for (int i = 0; i < N; i++) {
      if (vectVariables[i] == "") {
        vectVariables[i] = "vect-" + std::to_string(i);
      }
    }
  }
  else {
    for (int i = 0; i < N; i++) {
      vectVariables[i] = "vect-" + std::to_string(i);
    }
  }

  // Figure out the number of particles on each rank
  const unsigned long long numParticlesLocal  = a_particles.getNumberOfValidParticlesLocal();
  const unsigned long long numParticlesGlobal = a_particles.getNumberOfValidParticlesGlobal();

  std::vector<unsigned long long> particlesPerRank;
#ifdef CH_MPI
  particlesPerRank.resize(numProc(), 0ULL);

  std::vector<int> recv(numProc(), 1);
  std::vector<int> displ(numProc(), 0);
  for (int i = 0; i < numProc(); i++) {
    displ[i] = i;
  }
  MPI_Allgatherv(&numParticlesLocal,
                 1,
                 MPI_UNSIGNED_LONG_LONG,
                 &particlesPerRank[0],
                 &recv[0],
                 &displ[0],
                 MPI_UNSIGNED_LONG_LONG,
                 Chombo_MPI::comm);
#else
  particlesPerRank.resize(1);
  particlesPerRank[0] = numParticlesGlobal;
#endif

  // Set up file access and create the file.
  hid_t fileAccess = 0;
#ifdef CH_MPI
  fileAccess = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fileAccess, Chombo_MPI::comm, MPI_INFO_NULL);
#endif

  hid_t fileID = H5Fcreate(a_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fileAccess);
  H5Pclose(fileAccess);

  // Define the top group necessary for the H5Part file format
  hid_t grp = H5Gcreate2(fileID, "Step#0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Define dataspace dimensions.
  hsize_t dims[1];
  dims[0] = numParticlesGlobal;

  hid_t fileSpaceID = H5Screate_simple(1, dims, nullptr);

  // Memory space
  hsize_t memDims[1];
  memDims[0]       = numParticlesLocal;
  hid_t memSpaceID = H5Screate_simple(1, memDims, nullptr);

  // Set hyperslabs for file and memory
  hsize_t fileStart[1];
  hsize_t fileCount[1];

  hsize_t memStart[1];
  hsize_t memCount[1];

  memStart[0]  = 0;
  fileStart[0] = 0;
  for (int i = 0; i < procID(); i++) {
    fileStart[0] += particlesPerRank[i];
  }
  fileCount[0] = particlesPerRank[procID()];
  memCount[0]  = particlesPerRank[procID()];

  H5Sselect_hyperslab(fileSpaceID, H5S_SELECT_SET, fileStart, nullptr, fileCount, nullptr);
  H5Sselect_hyperslab(memSpaceID, H5S_SELECT_SET, memStart, nullptr, memCount, nullptr);

  // Create the ID and positional data sets
  hid_t datasetID = H5Dcreate2(grp, "id", H5T_NATIVE_ULLONG, fileSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t datasetX  = H5Dcreate2(grp, "x", H5T_NATIVE_DOUBLE, fileSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t datasetY  = H5Dcreate2(grp, "y", H5T_NATIVE_DOUBLE, fileSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#if CH_SPACEDIM == 3
  hid_t datasetZ = H5Dcreate2(grp, "z", H5T_NATIVE_DOUBLE, fileSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

  // Write ID data set
  std::vector<unsigned long long> id;
  std::vector<double>             x;
  std::vector<double>             y;
  std::vector<double>             z;

  for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
    for (DataIterator dit(a_particles.getGrids()[lvl]); dit.ok(); ++dit) {
      const List<GenericParticle<M, N>>& particles = a_particles[lvl][dit()].listItems();

      for (ListIterator<GenericParticle<M, N>> lit(particles); lit.ok(); ++lit) {
        id.push_back(procID());
        x.push_back(lit().position()[0] - a_shift[0]);
        y.push_back(lit().position()[1] - a_shift[0]);
#if CH_SPACEDIM == 3
        z.push_back(lit().position()[2] - a_shift[0]);
#endif
      }
    }
  }

  H5Dwrite(datasetID, H5Dget_type(datasetID), memSpaceID, fileSpaceID, H5P_DEFAULT, &id[0]);
  H5Dwrite(datasetX, H5Dget_type(datasetX), memSpaceID, fileSpaceID, H5P_DEFAULT, &x[0]);
  H5Dwrite(datasetY, H5Dget_type(datasetY), memSpaceID, fileSpaceID, H5P_DEFAULT, &y[0]);
#if CH_SPACEDIM == 3
  H5Dwrite(datasetZ, H5Dget_type(datasetY), memSpaceID, fileSpaceID, H5P_DEFAULT, &z[0]);
#endif

  id.resize(0);
  x.resize(0);
  y.resize(0);
  z.resize(0);

  // Close the ID and positional data sets
  H5Dclose(datasetID);
  H5Dclose(datasetX);
  H5Dclose(datasetY);
#if CH_SPACEDIM == 3
  H5Dclose(datasetZ);
#endif

  // Write the M real-variables
  for (int curVar = 0; curVar < M; curVar++) {
    hid_t dataset = H5Dcreate2(grp,
                               realVariables[curVar].c_str(),
                               H5T_NATIVE_DOUBLE,
                               fileSpaceID,
                               H5P_DEFAULT,
                               H5P_DEFAULT,
                               H5P_DEFAULT);

    std::vector<double> ds;

    for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
      for (DataIterator dit(a_particles.getGrids()[lvl]); dit.ok(); ++dit) {
        const List<GenericParticle<M, N>>& particles = a_particles[lvl][dit()].listItems();

        for (ListIterator<GenericParticle<M, N>> lit(particles); lit.ok(); ++lit) {
          ds.push_back(lit().getReals()[curVar]);
        }
      }
    }

    // Write and clsoe dataset
    H5Dwrite(dataset, H5Dget_type(dataset), memSpaceID, fileSpaceID, H5P_DEFAULT, &ds[0]);
    H5Dclose(dataset);
  }

  // Write the N vector variables
  for (int curVar = 0; curVar < N; curVar++) {
    for (int dir = 0; dir < SpaceDim; dir++) {
      std::string varName;
      if (dir == 0) {
        varName = vectVariables[curVar] + "-x";
      }
      else if (dir == 1) {
        varName = vectVariables[curVar] + "-y";
      }
      else if (dir == 2) {
        varName = vectVariables[curVar] + "-z";
      }

      hid_t dataset =
        H5Dcreate2(grp, varName.c_str(), H5T_NATIVE_DOUBLE, fileSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      std::vector<double> ds;

      for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
        for (DataIterator dit(a_particles.getGrids()[lvl]); dit.ok(); ++dit) {
          const List<GenericParticle<M, N>>& particles = a_particles[lvl][dit()].listItems();

          for (ListIterator<GenericParticle<M, N>> lit(particles); lit.ok(); ++lit) {
            const RealVect var = lit().getVects()[curVar];

            ds.push_back(var[dir]);
          }
        }
      }

      // Write and close dataset
      H5Dwrite(dataset, H5Dget_type(dataset), memSpaceID, fileSpaceID, H5P_DEFAULT, &ds[0]);
      H5Dclose(dataset);
    }
  }

  // Close top group and file
  H5Gclose(grp);
  H5Fclose(fileID);
}
#endif

#include <CD_NamespaceFooter.H>

#endif
