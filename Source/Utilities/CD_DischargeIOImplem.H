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
#ifdef CH_USE_HDF5
#include <hdf5.h>
#endif

// Chombo includes
#include <ParticleIO.H>

// Our includes
#include <CD_DischargeIO.H>
#include <CD_NamespaceHeader.H>

template <size_t M, size_t N>
void
DischargeIO::writeH5Part(const std::string                               a_filename,
                         const ParticleContainer<GenericParticle<M, N>>& a_particles,
                         const std::vector<std::string>                  a_realVars,
                         const std::vector<std::string>                  a_vectVars,
                         const RealVect                                  a_shift,
                         const Real                                      a_time) noexcept
{
#ifdef CH_USE_HDF5
  CH_TIME("DischargeIO::writeH5Part");

  CH_assert(a_realVars.size() == 0 || a_realVars.size() == M);
  CH_assert(a_vectVars.size() == 0 || a_vectVars.size() == N);

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

  // Write the time attribute
  hid_t  scal = H5Screate(H5S_SCALAR);
  hid_t  attr = H5Acreate(grp, "time", H5T_NATIVE_DOUBLE, scal, H5P_DEFAULT, H5P_DEFAULT);
  herr_t err  = H5Awrite(attr, H5T_NATIVE_DOUBLE, &a_time);

  H5Sclose(scal);
  H5Aclose(attr);

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
  hid_t datasetID = H5Dcreate2(grp, "rank id", H5T_NATIVE_ULLONG, fileSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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

    const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const List<GenericParticle<M, N>>& particles = a_particles[lvl][din].listItems();

      for (ListIterator<GenericParticle<M, N>> lit(particles); lit.ok(); ++lit) {
        id.push_back(procID());
        x.push_back(lit().position()[0] - a_shift[0]);
        y.push_back(lit().position()[1] - a_shift[1]);
#if CH_SPACEDIM == 3
        z.push_back(lit().position()[2] - a_shift[2]);
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

      const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const List<GenericParticle<M, N>>& particles = a_particles[lvl][din].listItems();

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
          const DataIndex& din = dit[mybox];

          const List<GenericParticle<M, N>>& particles = a_particles[lvl][din].listItems();

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
#endif
}

#ifdef CH_USE_HDF5
template <class P>
void
DischargeIO::writeCheckParticlesToHDF(HDF5Handle&            a_handle,
                                      const ParticleData<P>& a_particles,
                                      const std::string&     a_dataType)
{
  // timer
  CH_TIMERS("writeParticlesToHDF");

  // grids
  const BoxLayout& grids = a_particles.getBoxes();

  // loc part per box
  vector<unsigned long long> locParticlesPerBox(grids.size(), 0);

  // loc and tot number of particles
  unsigned long long numLocalParticles = 0;

  // count number of particles per box
  for (DataIterator di(grids); di.ok(); ++di) {
    const size_t numItems = a_particles[di].numItems();
    numLocalParticles += (unsigned long long)numItems;
    locParticlesPerBox[grids.index(di())] = (unsigned long long)numItems;
  }

  vector<unsigned long long> particlesPerBox(grids.size());

#ifdef CH_MPI
  int result = MPI_Allreduce(&locParticlesPerBox[0],
                             &particlesPerBox[0],
                             locParticlesPerBox.size(),
                             MPI_UNSIGNED_LONG_LONG,
                             MPI_SUM,
                             Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("MPI communcation error in ParticleIO");
  }

#else

  for (int i = 0; i < grids.size(); i++) {
    particlesPerBox[i] = (unsigned long long)locParticlesPerBox[i];
  }

#endif

  write_hdf_part_header(a_handle, grids, particlesPerBox, a_dataType);

  // size of individual objects
  size_t objSize = P().H5size();

  // the number of components per particle
  size_t numComps = objSize / sizeof(Real);

  // compute offsets and tot number of particles; order matters
  vector<unsigned long long> offsets;
  unsigned long long         totNumParticles = 0;
  for (int i = 0; i < grids.size(); i++) {
    offsets.push_back(numComps * totNumParticles);
    totNumParticles += particlesPerBox[i];
  }

  // now store particle info
  if (totNumParticles == 0) {
    return;
  }

  CH_TIMER("writeParticlesToHDF::checkpoint", t_cp);
  CH_START(t_cp);

  // data size
  hsize_t dataSize = totNumParticles * numComps;

  // a single dataspace
  hid_t dataspace = H5Screate_simple(1, &dataSize, NULL);

  // create data type
  hid_t H5T_type = H5T_NATIVE_DOUBLE;

  //
  std::string dataname = a_dataType + ":data";

  // open dataset
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
    // allocate data buffer where to linearize the particle data
    const size_t chunkSize = objSize * _CHUNK;
    char*        chunk     = new char[chunkSize];

    if (chunk == NULL) {
      MayDay::Error("WritePart::Error: new returned NULL pointer ");
    }

    for (DataIterator di(grids); di.ok(); ++di) {
      size_t offset = offsets[grids.index(di())];

      // particle counter
      int ip = 0;

      // linearize particle data
      const List<P>& pList = a_particles[di].listItems();

      for (ListIterator<P> li(pList); li.ok(); ++li) {
        pList[li].H5linearOut((void*)chunk);
        chunk += objSize;

        if (++ip % _CHUNK == 0) {
          // rewind pointer position, write buffered data
          chunk -= chunkSize;
          writeDataChunk(offset, dataspace, dataset, H5T_type, chunkSize / sizeof(Real), chunk);
        }
      }

      // write our residual particles
      const size_t nResidual = pList.length() % _CHUNK;
      if (nResidual > 0) {
        // rewind pointer position and write buffered data
        chunk -= (nResidual * objSize);
        writeDataChunk(offset, dataspace, dataset, H5T_type, nResidual * objSize / sizeof(Real), chunk);
      }
    }

    // free buffer
    delete[] chunk;
    chunk = NULL;
  }

  // done
  H5Sclose(dataspace);
  H5Dclose(dataset);

  // stop timer
  CH_STOP(t_cp);
}
#endif

#ifdef CH_USE_HDF5
template <class P>
void
DischargeIO::readCheckParticlesFromHDF(HDF5Handle&        a_handle,
                                       ParticleData<P>&   a_particles,
                                       const std::string& a_dataType)
{
  // timer
  CH_TIMERS("readParticlesFromHDF");

  // number of particles on each box
  vector<unsigned long long> particlesPerBox;

  // grids
  Vector<Box> grids;
  read_hdf_part_header(a_handle, grids, particlesPerBox, a_dataType, a_handle.getGroup());

  // load balancing should be handled outside
  const BoxLayout& bl = a_particles.getBoxes();

  // size of individual objects
  size_t objSize = P().H5size();

  // the number of components per particle
  size_t numComps = objSize / sizeof(Real);

  // compute offsets and tot number of particles; order matters
  vector<unsigned long long> offsets;
  unsigned long long         totNumParticles = 0;
  for (int i = 0; i < grids.size(); i++) {
    offsets.push_back(numComps * totNumParticles);
    totNumParticles += particlesPerBox[i];
  }

  // now read in particle data
  if (totNumParticles == 0) {
    return;
  }

  CH_TIMER("readParticlesFromHDF::checkpoint", t_cp);
  CH_START(t_cp);

  // P object
  P p;

  // data type
  hid_t H5T_type = H5T_NATIVE_DOUBLE;

  //
  std::string dataname = a_dataType + ":data";

  // one single dataset
#ifdef H516
  hid_t dataset = H5Dopen(a_handle.groupID(), dataname.c_str());
#else
  hid_t dataset = H5Dopen2(a_handle.groupID(), dataname.c_str(), H5P_DEFAULT);
#endif

  hid_t dataspace = H5Dget_space(dataset);

  // read in data relative to this proc in small chuncks
  const size_t chunkSize = objSize * _CHUNK;

  // allocate data buffer to read in linearized particle data chunks
  char* chunk = new char[chunkSize];

  for (DataIterator di(bl); di.ok(); ++di) {
    size_t boxData = numComps * particlesPerBox[bl.index(di())];
    size_t offset  = offsets[bl.index(di())];

    List<P>& items = a_particles[di].listItems();

    size_t dataIn = 0;
    while (dataIn < boxData) {
      // read data chunk of the right size
      int size = ((boxData - dataIn) >= (chunkSize / sizeof(Real))) ? chunkSize / sizeof(Real) : boxData - dataIn;
      readDataChunk(offset, dataspace, dataset, H5T_type, size, chunk);

      // list to store the particle data from the HDF5 file
      for (int ip = 0; ip < size / numComps; ip++) {
        p.H5linearIn(chunk);
        chunk += objSize;

        // assign to data holder
        items.add(p);
      }
      dataIn += size;
      chunk -= size * sizeof(Real);
    }
  }

  // free alllocated memory
  delete[] chunk;
  chunk = NULL;

  // done with this component
  H5Sclose(dataspace);
  H5Dclose(dataset);

  // stop timer
  CH_STOP(t_cp);
}
#endif

#include <CD_NamespaceFooter.H>

#endif
