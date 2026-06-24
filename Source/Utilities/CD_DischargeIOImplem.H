/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_DischargeIOImplem.H
  @brief  Implementation of CD_DischargeIO.H
  @author Robert Marskar
*/

#ifndef CD_DISCHARGEIOIMPLEM_H
#define CD_DISCHARGEIOIMPLEM_H

// Std includes
#ifdef CH_USE_HDF5
#include <hdf5.h>
#endif

// Chombo includes
#include <ParticleIO.H>

// Our includes
#include <CD_DischargeIO.H>
#include <CD_NamespaceHeader.H>

template <typename P, typename Traits>
void
DischargeIO::writeH5Part(std::string                                                  a_filename,
                         const ParticleContainer<P, Traits>&                          a_particles,
                         std::vector<std::pair<std::string, H5PartScalar<P, Traits>>> a_scalarVars,
                         std::vector<std::pair<std::string, H5PartVect<P, Traits>>>   a_vectorVars,
                         RealVect                                                     a_shift,
                         Real                                                         a_time) noexcept
{
#ifdef CH_USE_HDF5
  CH_TIME("DischargeIO::writeH5Part(SoA)");

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
#ifdef CH_MPI
  H5Pclose(fileAccess);
#endif

  // Define the top group necessary for the H5Part file format
  hid_t grp = H5Gcreate2(fileID, "Step#0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Write the time attribute
  hid_t scal     = H5Screate(H5S_SCALAR);
  hid_t timeAttr = H5Acreate(grp, "time", H5T_NATIVE_DOUBLE, scal, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(timeAttr, H5T_NATIVE_DOUBLE, &a_time);
  H5Aclose(timeAttr);

  H5Sclose(scal);

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
  hid_t datasetID = H5Dcreate2(grp, "id", H5T_NATIVE_LLONG, fileSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t datasetX  = H5Dcreate2(grp, "x", H5T_NATIVE_DOUBLE, fileSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t datasetY  = H5Dcreate2(grp, "y", H5T_NATIVE_DOUBLE, fileSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#if CH_SPACEDIM == 3
  hid_t datasetZ = H5Dcreate2(grp, "z", H5T_NATIVE_DOUBLE, fileSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

  std::vector<long long> id;
  std::vector<double>    x;
  std::vector<double>    y;
  std::vector<double>    z;

  for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

      for (std::size_t i = 0; i < leaf.size(); i++) {
        const RealVect pos = leaf.position(i);

        id.push_back(static_cast<long long>(leaf.particleID(i)));
        x.push_back(pos[0] - a_shift[0]);
        y.push_back(pos[1] - a_shift[1]);
#if CH_SPACEDIM == 3
        z.push_back(pos[2] - a_shift[2]);
#endif
      }
    }
  }

  H5Dwrite(datasetID, H5Dget_type(datasetID), memSpaceID, fileSpaceID, H5P_DEFAULT, id.data());
  H5Dwrite(datasetX, H5Dget_type(datasetX), memSpaceID, fileSpaceID, H5P_DEFAULT, x.data());
  H5Dwrite(datasetY, H5Dget_type(datasetY), memSpaceID, fileSpaceID, H5P_DEFAULT, y.data());
#if CH_SPACEDIM == 3
  H5Dwrite(datasetZ, H5Dget_type(datasetZ), memSpaceID, fileSpaceID, H5P_DEFAULT, z.data());
#endif

  H5Dclose(datasetID);
  H5Dclose(datasetX);
  H5Dclose(datasetY);
#if CH_SPACEDIM == 3
  H5Dclose(datasetZ);
#endif

  // Write the scalar payload datasets (caller-named; e.g. weight + payload scalar columns).
  for (const auto& var : a_scalarVars) {
    hid_t dataset = H5Dcreate2(grp,
                               var.first.c_str(),
                               H5T_NATIVE_DOUBLE,
                               fileSpaceID,
                               H5P_DEFAULT,
                               H5P_DEFAULT,
                               H5P_DEFAULT);

    std::vector<double> ds;

    for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&              din  = dit[mybox];
        const ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

        for (std::size_t i = 0; i < leaf.size(); i++) {
          ds.push_back(var.second(leaf, i));
        }
      }
    }

    H5Dwrite(dataset, H5Dget_type(dataset), memSpaceID, fileSpaceID, H5P_DEFAULT, ds.data());
    H5Dclose(dataset);
  }

  // Write the vector payload datasets as name-x/name-y/name-z.
  for (const auto& var : a_vectorVars) {
    for (int dir = 0; dir < SpaceDim; dir++) {
      std::string varName = var.first;
      if (dir == 0) {
        varName += "-x";
      }
      else if (dir == 1) {
        varName += "-y";
      }
      else if (dir == 2) {
        varName += "-z";
      }

      hid_t dataset = H5Dcreate2(grp,
                                 varName.c_str(),
                                 H5T_NATIVE_DOUBLE,
                                 fileSpaceID,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT);

      std::vector<double> ds;

      for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
        const DataIterator&      dit = dbl.dataIterator();

        const int nbox = dit.size();

        for (int mybox = 0; mybox < nbox; mybox++) {
          const DataIndex&              din  = dit[mybox];
          const ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

          for (std::size_t i = 0; i < leaf.size(); i++) {
            ds.push_back(var.second(leaf, i)[dir]);
          }
        }
      }

      H5Dwrite(dataset, H5Dget_type(dataset), memSpaceID, fileSpaceID, H5P_DEFAULT, ds.data());
      H5Dclose(dataset);
    }
  }

  // Close top group and file
  H5Gclose(grp);
  H5Fclose(fileID);
#endif
}

#ifdef CH_USE_HDF5
template <typename P, typename Traits>
void
DischargeIO::writeCheckParticlesToHDF(HDF5Handle&                               a_handle,
                                      const LayoutData<ParticleSoA<P, Traits>>& a_particles,
                                      const std::string&                        a_dataType)
{
  CH_TIMERS("writeParticlesToHDF(SoA)");

  const BoxLayout& grids = a_particles.boxLayout();

  std::vector<unsigned long long> locParticlesPerBox(grids.size(), 0);
  unsigned long long              numLocalParticles = 0;

  for (DataIterator di(grids); di.ok(); ++di) {
    const size_t numItems = a_particles[di].size();
    numLocalParticles += (unsigned long long)numItems;
    locParticlesPerBox[grids.index(di())] = (unsigned long long)numItems;
  }

  std::vector<unsigned long long> particlesPerBox(grids.size());

#ifdef CH_MPI
  int result = MPI_Allreduce(&locParticlesPerBox[0],
                             &particlesPerBox[0],
                             locParticlesPerBox.size(),
                             MPI_UNSIGNED_LONG_LONG,
                             MPI_SUM,
                             Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("MPI communication error in ParticleIO (SoA)");
  }
#else
  for (int i = 0; i < grids.size(); i++) {
    particlesPerBox[i] = (unsigned long long)locParticlesPerBox[i];
  }
#endif

  write_hdf_part_header(a_handle, grids, particlesPerBox, a_dataType);

  // The SoA HDF5 record is position + weight + the payload h5 subset, all doubles.
  const size_t objSize  = ParticleSoA<P, Traits>::h5BytesPerParticle();
  const size_t numComps = objSize / sizeof(double);

  std::vector<unsigned long long> offsets;
  unsigned long long              totNumParticles = 0;
  for (int i = 0; i < grids.size(); i++) {
    offsets.push_back(numComps * totNumParticles);
    totNumParticles += particlesPerBox[i];
  }

  if (totNumParticles == 0) {
    return;
  }

  hsize_t     dataSize  = totNumParticles * numComps;
  hid_t       dataspace = H5Screate_simple(1, &dataSize, nullptr);
  hid_t       H5T_type  = H5T_NATIVE_DOUBLE;
  std::string dataname  = a_dataType + ":data";

#ifdef H516
  hid_t dataset = H5Dcreate(a_handle.groupID(), dataname.c_str(), H5T_type, dataspace, H5P_DEFAULT);
#else
  hid_t dataset = H5Dcreate2(a_handle.groupID(),
                             dataname.c_str(),
                             H5T_type,
                             dataspace,
                             H5P_DEFAULT,
                             H5P_DEFAULT,
                             H5P_DEFAULT);
#endif

  if (numLocalParticles > 0) {
    const size_t chunkSize = objSize * _CHUNK;
    char* const  chunkBase = new char[chunkSize];
    char*        chunk     = chunkBase;

    if (chunkBase == nullptr) {
      MayDay::Error("WritePart(SoA)::Error: new returned NULL pointer ");
    }

    for (DataIterator di(grids); di.ok(); ++di) {
      size_t offset = offsets[grids.index(di())];

      const ParticleSoA<P, Traits>& leaf = a_particles[di];
      const size_t                  n    = leaf.size();

      for (size_t ip = 0; ip < n; ip++) {
        leaf.h5LinearizeParticle((void*)chunk, ip);
        chunk += objSize;

        if ((ip + 1) % _CHUNK == 0) {
          chunk -= chunkSize;
          writeDataChunk(offset, dataspace, dataset, H5T_type, chunkSize / sizeof(double), chunk);
        }
      }

      const size_t nResidual = n % _CHUNK;
      if (nResidual > 0) {
        chunk -= (nResidual * objSize);
        writeDataChunk(offset, dataspace, dataset, H5T_type, nResidual * objSize / sizeof(double), chunk);
      }

      chunk = chunkBase;
    }

    delete[] chunkBase;
  }

  H5Sclose(dataspace);
  H5Dclose(dataset);
}
#endif

#ifdef CH_USE_HDF5
template <typename P, typename Traits>
void
DischargeIO::readCheckParticlesFromHDF(HDF5Handle&                         a_handle,
                                       LayoutData<ParticleSoA<P, Traits>>& a_particles,
                                       const std::string&                  a_dataType)
{
  CH_TIMERS("readParticlesFromHDF(SoA)");

  std::vector<unsigned long long> particlesPerBox;
  Vector<Box>                     grids;
  read_hdf_part_header(a_handle, grids, particlesPerBox, a_dataType, a_handle.getGroup());

  const BoxLayout& bl = a_particles.boxLayout();

  const size_t objSize  = ParticleSoA<P, Traits>::h5BytesPerParticle();
  const size_t numComps = objSize / sizeof(double);

  std::vector<unsigned long long> offsets;
  unsigned long long              totNumParticles = 0;
  for (int i = 0; i < grids.size(); i++) {
    offsets.push_back(numComps * totNumParticles);
    totNumParticles += particlesPerBox[i];
  }

  if (totNumParticles == 0) {
    return;
  }

  hid_t       H5T_type = H5T_NATIVE_DOUBLE;
  std::string dataname = a_dataType + ":data";

#ifdef H516
  hid_t dataset = H5Dopen(a_handle.groupID(), dataname.c_str());
#else
  hid_t dataset = H5Dopen2(a_handle.groupID(), dataname.c_str(), H5P_DEFAULT);
#endif

  hid_t dataspace = H5Dget_space(dataset);

  const size_t chunkSize = objSize * _CHUNK;
  char*        chunk     = new char[chunkSize];

  for (DataIterator di(bl); di.ok(); ++di) {
    size_t boxData = numComps * particlesPerBox[bl.index(di())];
    size_t offset  = offsets[bl.index(di())];

    ParticleSoA<P, Traits>& leaf = a_particles[di];

    size_t dataIn = 0;
    while (dataIn < boxData) {
      int size = ((boxData - dataIn) >= (chunkSize / sizeof(double))) ? chunkSize / sizeof(double) : boxData - dataIn;
      readDataChunk(offset, dataspace, dataset, H5T_type, size, chunk);

      for (int ip = 0; ip < size / (int)numComps; ip++) {
        leaf.h5DelinearizeAndAppend(chunk);
        chunk += objSize;
      }
      dataIn += size;
      chunk -= size * sizeof(double);
    }
  }

  delete[] chunk;
  chunk = nullptr;

  H5Sclose(dataspace);
  H5Dclose(dataset);
}
#endif

#include <CD_NamespaceFooter.H>

#endif
