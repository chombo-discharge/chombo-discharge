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

#ifdef CH_USE_HDF5
namespace DischargeIO {

  /**
    @brief Detector (false case) for an optional Traits::h5PartColumns descriptor tuple.
    @tparam Traits Column descriptor for the payload.
  */
  template <typename Traits, typename = void>
  struct HasH5PartColumns : std::false_type
  {
  };

  /**
    @brief Detector (true case) for an optional Traits::h5PartColumns descriptor tuple.
    @tparam Traits Column descriptor declaring h5PartColumns.
  */
  template <typename Traits>
  struct HasH5PartColumns<Traits, detail::Void<decltype(Traits::h5PartColumns)>> : std::true_type
  {
  };

  /**
    @brief Write one named scalar dataset (one value per particle) into an open H5Part step group.
    @details Gathers the value for every valid particle across all AMR levels (in the same per-level,
    per-box, per-particle order as the id/position datasets) and writes it as an H5T_NATIVE_DOUBLE dataset.
    @param[in] a_grp        Open H5Part step group.
    @param[in] a_fileSpace  File dataspace (global particle count, hyperslab pre-selected).
    @param[in] a_memSpace   Memory dataspace (local particle count, hyperslab pre-selected).
    @param[in] a_particles  SoA particle container.
    @param[in] a_name       Dataset name.
    @param[in] a_gather     Callable (const ParticleSoA<P,Traits>&, std::size_t) -> double.
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
    @tparam Gather Gather-callable type.
  */
  template <typename P, typename Traits, typename Gather>
  inline void
  writeH5PartScalarDataset(hid_t                               a_grp,
                           hid_t                               a_fileSpace,
                           hid_t                               a_memSpace,
                           const ParticleContainer<P, Traits>& a_particles,
                           const std::string&                  a_name,
                           Gather                              a_gather) noexcept
  {
    hid_t dataset = H5Dcreate2(a_grp,
                               a_name.c_str(),
                               H5T_NATIVE_DOUBLE,
                               a_fileSpace,
                               H5P_DEFAULT,
                               H5P_DEFAULT,
                               H5P_DEFAULT);

    std::vector<double> ds;

    for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl  = a_particles.getGrids()[lvl];
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();

      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&              din  = dit[mybox];
        const ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

        for (std::size_t i = 0; i < leaf.size(); i++) {
          ds.push_back(a_gather(leaf, i));
        }
      }
    }

    H5Dwrite(dataset, H5Dget_type(dataset), a_memSpace, a_fileSpace, H5P_DEFAULT, ds.data());
    H5Dclose(dataset);
  }

  /**
    @brief Write the SpaceDim component datasets (name-x/name-y/name-z) of one H5Part::Vector descriptor.
    @param[in] a_grp        Open H5Part step group.
    @param[in] a_fileSpace  File dataspace.
    @param[in] a_memSpace   Memory dataspace.
    @param[in] a_particles  SoA particle container.
    @param[in] a_dirs       Index sequence [0, SpaceDim) of the spatial components.
    @tparam I      Index of the descriptor within Traits::h5PartColumns.
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
    @tparam Dirs   Spatial component indices [0, SpaceDim).
  */
  template <std::size_t I, typename P, typename Traits, std::size_t... Dirs>
  inline void
  writeH5PartVectorColumn(hid_t                                         a_grp,
                          hid_t                                         a_fileSpace,
                          hid_t                                         a_memSpace,
                          const ParticleContainer<P, Traits>&           a_particles,
                          [[maybe_unused]] std::index_sequence<Dirs...> a_dirs) noexcept
  {
    constexpr auto    desc      = std::get<I>(Traits::h5PartColumns);
    const char* const suffix[3] = {"-x", "-y", "-z"};

    (writeH5PartScalarDataset(a_grp,
                              a_fileSpace,
                              a_memSpace,
                              a_particles,
                              std::string(desc.name) + suffix[Dirs],
                              [](const ParticleSoA<P, Traits>& a_leaf, std::size_t a_i) -> double {
                                constexpr auto member = std::get<Dirs>(std::get<I>(Traits::h5PartColumns).members);
                                return static_cast<double>(a_leaf.template get<member>(a_i));
                              }),
     ...);
  }

  /**
    @brief Write one payload descriptor of Traits::h5PartColumns (scalar -> one dataset, vector ->
           SpaceDim datasets).
    @param[in] a_grp        Open H5Part step group.
    @param[in] a_fileSpace  File dataspace.
    @param[in] a_memSpace   Memory dataspace.
    @param[in] a_particles  SoA particle container.
    @tparam I      Index of the descriptor within Traits::h5PartColumns.
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
  */
  template <std::size_t I, typename P, typename Traits>
  inline void
  writeH5PartColumn(hid_t                               a_grp,
                    hid_t                               a_fileSpace,
                    hid_t                               a_memSpace,
                    const ParticleContainer<P, Traits>& a_particles) noexcept
  {
    constexpr auto desc = std::get<I>(Traits::h5PartColumns);

    if constexpr (decltype(desc)::s_isVector) {
      writeH5PartVectorColumn<I>(a_grp, a_fileSpace, a_memSpace, a_particles, std::make_index_sequence<SpaceDim>{});
    }
    else {
      writeH5PartScalarDataset(a_grp,
                               a_fileSpace,
                               a_memSpace,
                               a_particles,
                               desc.name,
                               [](const ParticleSoA<P, Traits>& a_leaf, std::size_t a_i) -> double {
                                 constexpr auto member = std::get<I>(Traits::h5PartColumns).member;
                                 return static_cast<double>(a_leaf.template get<member>(a_i));
                               });
    }
  }

  /**
    @brief Write every payload descriptor in Traits::h5PartColumns.
    @param[in] a_grp        Open H5Part step group.
    @param[in] a_fileSpace  File dataspace.
    @param[in] a_memSpace   Memory dataspace.
    @param[in] a_particles  SoA particle container.
    @param[in] a_indices    Index sequence [0, size(h5PartColumns)).
    @tparam P      Payload type of the container.
    @tparam Traits Column descriptor for the payload.
    @tparam I      Descriptor indices [0, size(h5PartColumns)).
  */
  template <typename P, typename Traits, std::size_t... I>
  inline void
  writeH5PartColumns(hid_t                                      a_grp,
                     hid_t                                      a_fileSpace,
                     hid_t                                      a_memSpace,
                     const ParticleContainer<P, Traits>&        a_particles,
                     [[maybe_unused]] std::index_sequence<I...> a_indices) noexcept
  {
    (writeH5PartColumn<I>(a_grp, a_fileSpace, a_memSpace, a_particles), ...);
  }

} // namespace DischargeIO
#endif

template <typename P, typename Traits>
void
DischargeIO::writeH5Part(std::string                         a_filename,
                         const ParticleContainer<P, Traits>& a_particles,
                         RealVect                            a_shift,
                         Real                                a_time) noexcept
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

  // Weight is a container-owned column and is always written.
  writeH5PartScalarDataset(grp,
                           fileSpaceID,
                           memSpaceID,
                           a_particles,
                           "weight",
                           [](const ParticleSoA<P, Traits>& a_leaf, std::size_t a_i) -> double {
                             return a_leaf.weight(a_i);
                           });

  // Payload datasets, derived declaratively from the optional Traits::h5PartColumns descriptor tuple.
  if constexpr (HasH5PartColumns<Traits>::value) {
    using H5PartTuple = std::remove_cv_t<std::remove_reference_t<decltype(Traits::h5PartColumns)>>;
    writeH5PartColumns(grp,
                       fileSpaceID,
                       memSpaceID,
                       a_particles,
                       std::make_index_sequence<std::tuple_size<H5PartTuple>::value>{});
  }

  // Close dataspaces, top group and file
  H5Sclose(fileSpaceID);
  H5Sclose(memSpaceID);
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

  // serial: accumulates into the running counter numLocalParticles
  const DataIterator& dit  = grids.dataIterator();
  const int           nbox = dit.size();
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din      = dit[mybox];
    const size_t     numItems = a_particles[din].size();
    numLocalParticles += (unsigned long long)numItems;
    locParticlesPerBox[grids.index(din)] = (unsigned long long)numItems;
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

  // The SoA HDF5 record (position + weight + the payload h5 subset) is stored as raw bytes
  // (H5T_NATIVE_CHAR). Payload columns may be ParticleReal (float), so the per-particle byte size need
  // NOT be a multiple of sizeof(double); a double-typed dataset would truncate the extent
  // (numComps = objSize/sizeof(double)) and overflow the hyperslab on write/read. Byte-based storage is
  // exact for any payload and any SpaceDim.
  const size_t objSize = ParticleSoA<P, Traits>::h5BytesPerParticle();

  std::vector<unsigned long long> offsets;
  unsigned long long              totNumParticles = 0;
  for (int i = 0; i < grids.size(); i++) {
    offsets.push_back(objSize * totNumParticles);
    totNumParticles += particlesPerBox[i];
  }

  if (totNumParticles == 0) {
    return;
  }

  hsize_t     dataSize  = totNumParticles * objSize;
  hid_t       dataspace = H5Screate_simple(1, &dataSize, nullptr);
  hid_t       H5T_type  = H5T_NATIVE_CHAR;
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

    // serial: writes into the shared HDF5 dataset at running per-box offsets and shares the chunk buffer
    const DataIterator& dit  = grids.dataIterator();
    const int           nbox = dit.size();

    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din    = dit[mybox];
      size_t           offset = offsets[grids.index(din)];

      const ParticleSoA<P, Traits>& leaf = a_particles[din];
      const size_t                  n    = leaf.size();

      for (size_t ip = 0; ip < n; ip++) {
        leaf.h5LinearizeParticle((void*)chunk, ip);
        chunk += objSize;

        if ((ip + 1) % _CHUNK == 0) {
          chunk -= chunkSize;
          writeDataChunk(offset, dataspace, dataset, H5T_type, chunkSize, chunk);
        }
      }

      const size_t nResidual = n % _CHUNK;
      if (nResidual > 0) {
        chunk -= (nResidual * objSize);
        writeDataChunk(offset, dataspace, dataset, H5T_type, nResidual * objSize, chunk);
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

  // Byte-based storage (matches writeCheckParticlesToHDF): the per-particle record may not be a multiple
  // of sizeof(double) (ParticleReal payload columns), so everything is counted in bytes (H5T_NATIVE_CHAR).
  const size_t objSize = ParticleSoA<P, Traits>::h5BytesPerParticle();

  std::vector<unsigned long long> offsets;
  unsigned long long              totNumParticles = 0;
  for (int i = 0; i < grids.size(); i++) {
    offsets.push_back(objSize * totNumParticles);
    totNumParticles += particlesPerBox[i];
  }

  if (totNumParticles == 0) {
    return;
  }

  hid_t       H5T_type = H5T_NATIVE_CHAR;
  std::string dataname = a_dataType + ":data";

#ifdef H516
  hid_t dataset = H5Dopen(a_handle.groupID(), dataname.c_str());
#else
  hid_t dataset = H5Dopen2(a_handle.groupID(), dataname.c_str(), H5P_DEFAULT);
#endif

  hid_t dataspace = H5Dget_space(dataset);

  const size_t chunkSize = objSize * _CHUNK;
  char*        chunk     = new char[chunkSize];

  // serial: reads from the shared HDF5 dataset at per-box offsets and shares the chunk buffer
  const DataIterator& dit  = bl.dataIterator();
  const int           nbox = dit.size();
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din     = dit[mybox];
    size_t           boxData = objSize * particlesPerBox[bl.index(din)]; // bytes
    size_t           offset  = offsets[bl.index(din)];                   // bytes

    ParticleSoA<P, Traits>& leaf = a_particles[din];

    size_t dataIn = 0;
    while (dataIn < boxData) {
      int size = ((boxData - dataIn) >= chunkSize) ? chunkSize : boxData - dataIn; // bytes
      readDataChunk(offset, dataspace, dataset, H5T_type, size, chunk);

      for (int ip = 0; ip < size / (int)objSize; ip++) {
        leaf.h5DelinearizeAndAppend(chunk);
        chunk += objSize;
      }
      dataIn += size;
      chunk -= size;
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
