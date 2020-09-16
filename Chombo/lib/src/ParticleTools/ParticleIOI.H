#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// functions for I/O of particle data

#ifndef _PARTICLEIOI_H_
#define _PARTICLEIOI_H_

#include "ParticleIO.H"

#include "NamespaceHeader.H"

// save data in small chunks
#define _CHUNK (1024*1024)

#ifdef CH_USE_HDF5

template <class T>
void write_vect_to_header(HDF5Handle&         a_handle,
                          const vector<T>&    a_vect,
                          const hid_t&        H5T_type,
                          const std::string&  a_dataname)
{
  // store vect info
  hsize_t length = a_vect.size();
  hid_t dataspace= H5Screate_simple(1, &length, NULL);

#ifdef H516
  hid_t dataset  = H5Dcreate(a_handle.groupID(), a_dataname.c_str(),
                             H5T_type, dataspace,
                             H5P_DEFAULT);
#else
  hid_t dataset  = H5Dcreate2(a_handle.groupID(), a_dataname.c_str(),
			      H5T_type, dataspace,
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif  

  CH_assert(dataspace >= 0);
  CH_assert(dataset >= 0);
  if (procID() == 0)
    {
      hid_t memdataspace = H5Screate_simple(1, &length, NULL);
      CH_assert(memdataspace >= 0);
      int err = H5Dwrite(dataset, H5T_type, memdataspace,
                         dataspace, H5P_DEFAULT, &a_vect[0]);
      CH_assert(err >= 0);
      H5Sclose(memdataspace);
    }
  H5Sclose(dataspace);
  H5Dclose(dataset);
}

// read vect to header
template <class T>
void read_vect_from_header(HDF5Handle&        a_handle,
                          vector<T>&          a_vect,
                          const hid_t&        H5T_type,
                          const std::string&  a_dataname)
{
  // read in offset info
  // hsize_t length = sizeof(T) * a_vect.size();
  hsize_t length = a_vect.size();

  hid_t dataspace= H5Screate_simple(1, &length, NULL);
#ifdef H516
  hid_t dataset = H5Dopen(a_handle.groupID(), a_dataname.c_str());
#else
  hid_t dataset = H5Dopen2(a_handle.groupID(), a_dataname.c_str(), H5P_DEFAULT);
#endif 
  CH_assert(dataspace >= 0);
  CH_assert(dataset >= 0);

  hid_t memdataspace = H5Screate_simple(1, &length, NULL);
  CH_assert(memdataspace >= 0);
  int err = H5Dread(dataset, H5T_type, memdataspace,
                    dataspace, H5P_DEFAULT, &a_vect[0]);
  CH_assert(err >= 0);

  H5Sclose(memdataspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
}

// write whole hierarchy of particle data
template<class P> void
writeParticlesToHDF(HDF5Handle&            a_handle,
                    const ParticleData<P>& a_particles,
                    const std::string&     a_dataType)
{
  // timer
  CH_TIMERS("writeParticlesToHDF");

  // grids
  const BoxLayout& grids = a_particles.getBoxes();

  // loc part per box
  vector<unsigned long long> locParticlesPerBox(grids.size(),0);

  // loc and tot number of particles
  unsigned long long numLocalParticles = 0;

  // count number of particles per box
  for (DataIterator di(grids); di.ok(); ++di)
    {
      const size_t numItems=a_particles[di].numItems();
      numLocalParticles += (unsigned long long)numItems;
      locParticlesPerBox[grids.index(di())] = (unsigned long long)numItems;
    }

  vector<unsigned long long> particlesPerBox(grids.size());

#ifdef CH_MPI
  int result = MPI_Allreduce(&locParticlesPerBox[0], &particlesPerBox[0], locParticlesPerBox.size(),
                             MPI_UNSIGNED_LONG_LONG, MPI_SUM, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    {
      MayDay::Error("MPI communcation error in ParticleIO");
    }

#else

  for (int i=0; i<grids.size(); i++)
    {
      particlesPerBox[i]=(unsigned long long)locParticlesPerBox[i];
    }

#endif

  write_hdf_part_header(a_handle, grids, particlesPerBox, a_dataType);

  // size of individual objects
  size_t objSize = P().size();

  // the number of components per particle
  size_t numComps = objSize / sizeof(Real);

  // compute offsets and tot number of particles; order matters
  vector<unsigned long long> offsets;
  unsigned long long totNumParticles = 0;
  for (int i=0; i<grids.size(); i++)
    {
      offsets.push_back(numComps * totNumParticles);
      totNumParticles += particlesPerBox[i];
    }

  // now store particle info
  if (totNumParticles==0)
    {
      return;
    }

  CH_TIMER("writeParticlesToHDF::checkpoint",t_cp);
  CH_START(t_cp);

  // data size
  hsize_t dataSize = totNumParticles * numComps;

  // a single dataspace
  hid_t dataspace = H5Screate_simple(1, &dataSize, NULL);

  // create data type
  hid_t H5T_type = H5T_NATIVE_DOUBLE;

  //
  std::string dataname = a_dataType+":data";

  // open dataset
#ifdef H516
  hid_t dataset  = H5Dcreate(a_handle.groupID(), dataname.c_str(),
                             H5T_type, dataspace,
                             H5P_DEFAULT);
#else
  hid_t dataset  = H5Dcreate2(a_handle.groupID(), dataname.c_str(),
			      H5T_type, dataspace,
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif  

  if (numLocalParticles > 0)
    {
      // allocate data buffer where to linearize the particle data
      const size_t chunkSize = objSize*_CHUNK;
      char* chunk = new char[chunkSize];

      if (chunk == NULL) {
        MayDay::Error("WritePart::Error: new returned NULL pointer ");
      }

      for (DataIterator di(grids); di.ok(); ++di)
        {
          size_t offset = offsets[ grids.index(di()) ];

          // particle counter
          int ip = 0;

          // linearize particle data
          const List<P>& pList = a_particles[di].listItems();

          for (ListIterator<P> li(pList); li.ok(); ++li)
            {
              pList[li].linearOut((void*)chunk);
              chunk += objSize;

              if (++ip % _CHUNK == 0)
                {
                  // rewind pointer position, write buffered data
                  chunk -= chunkSize;
                  writeDataChunk(offset,dataspace,dataset,H5T_type,
                                 chunkSize / sizeof(Real),chunk);
                }
            }

          // write our residual particles
          const size_t nResidual = pList.length()%_CHUNK;
          if (nResidual>0)
            {
              // rewind pointer position and write buffered data
              chunk -= (nResidual*objSize);
              writeDataChunk(offset,dataspace,dataset,H5T_type,
                             nResidual*objSize / sizeof(Real),chunk);
            }
        }

      // free buffer
      delete[] chunk; chunk=NULL;
    }

  // done
  H5Sclose(dataspace);
  H5Dclose(dataset);

  // stop timer
  CH_STOP(t_cp);
}

// read whole set of particle data
template<class P>
void readParticlesFromHDF(HDF5Handle&        a_handle,
                          ParticleData<P>&   a_particles,
                          const std::string& a_dataType)
{
  // timer
  CH_TIMERS("readParticlesFromHDF");

  // number of particles on each box
  vector<unsigned long long> particlesPerBox;

  // grids
  Vector<Box> grids;
  read_hdf_part_header(a_handle,
                       grids,
                       particlesPerBox,
                       a_dataType,
                       a_handle.getGroup());

   // load balancing should be handled outside
   const BoxLayout& bl = a_particles.getBoxes();

  // size of individual objects
  size_t objSize = P().size();

  // the number of components per particle
  size_t numComps = objSize / sizeof(Real);

  // compute offsets and tot number of particles; order matters
  vector<unsigned long long> offsets;
  unsigned long long totNumParticles = 0;
  for (int i=0; i<grids.size(); i++)
    {
      offsets.push_back(numComps * totNumParticles);
      totNumParticles += particlesPerBox[i];
    }

  // now read in particle data
  if (totNumParticles==0)
    {
      return;
    }

  CH_TIMER("readParticlesFromHDF::checkpoint",t_cp);
  CH_START(t_cp);

  // P object
  P p;

  // data type
  hid_t H5T_type = H5T_NATIVE_DOUBLE;

  //
  std::string dataname = a_dataType+":data";

  // one single dataset
#ifdef H516
  hid_t dataset = H5Dopen(a_handle.groupID(), dataname.c_str());
#else
  hid_t dataset = H5Dopen2(a_handle.groupID(), dataname.c_str(), H5P_DEFAULT);
#endif 

  hid_t dataspace= H5Dget_space(dataset);

  // read in data relative to this proc in small chuncks
  const size_t chunkSize = objSize*_CHUNK;

  // allocate data buffer to read in linearized particle data chunks
  char* chunk = new char[chunkSize];

  for (DataIterator di(bl); di.ok(); ++di)
    {
      size_t boxData= numComps * particlesPerBox[ bl.index(di()) ];
      size_t offset = offsets[ bl.index(di()) ];

      List<P>& items= a_particles[di].listItems();

      size_t dataIn = 0;
      while (dataIn < boxData)
        {
          // read data chunk of the right size
          int size = ((boxData-dataIn) >= (chunkSize / sizeof(Real))) ? chunkSize / sizeof(Real) : boxData - dataIn;
          readDataChunk(offset, dataspace, dataset, H5T_type, size, chunk);

          // list to store the particle data from the HDF5 file
          for (int ip=0; ip<size/numComps; ip++)
            {
              p.linearIn(chunk);
              chunk += objSize;

              // assign to data holder
              items.add(p);
            }
          dataIn += size;
          chunk  -= size*sizeof(Real);
        }
    }

  // free alllocated memory
  delete[] chunk; chunk=NULL;

  // done with this component
  H5Sclose(dataspace);
  H5Dclose(dataset);

  // stop timer
  CH_STOP(t_cp);
}

// write whole hierarchy of particle data
template<class P> void
writeParticlesToHDF(HDF5Handle&            a_handle,
                    const Vector<P>&       a_particles,
		    const Box&             a_domain,
                    const std::string&     a_dataType)
{
  // timer
  CH_TIMERS("writeParticlesToHDF");

  // We make a fake BoxLayout with only one Box.
  Vector<Box> boxes;
  boxes.push_back(a_domain);
  Vector<int> procIDs;
  procIDs.push_back(0);
  BoxLayout grids(boxes, procIDs);

  unsigned long long totNumParticles = a_particles.size();

  // In this case, there is only one box
  vector<unsigned long long> particlesPerBox;
  particlesPerBox.push_back(totNumParticles);

  write_hdf_part_header(a_handle, grids, particlesPerBox, a_dataType);

  // size of individual objects
  size_t objSize = P().size();

  // the number of components per particle
  size_t numComps = objSize / sizeof(Real);

    // don't do anything if there are no particles
  if (totNumParticles==0)
    {
      return;
    }

  // data size
  hsize_t dataSize = totNumParticles * numComps;

  // a single dataspace
  hid_t dataspace = H5Screate_simple(1, &dataSize, NULL);

  // create data type
  hid_t H5T_type = H5T_NATIVE_DOUBLE;

  // 
  std::string dataname = a_dataType+":data";

  // open dataset

#ifdef H516
  hid_t dataset  = H5Dcreate(a_handle.groupID(), dataname.c_str(),
                             H5T_type, dataspace,
                             H5P_DEFAULT);
#else
  hid_t dataset  = H5Dcreate2(a_handle.groupID(), dataname.c_str(),
			      H5T_type, dataspace,
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif  

  // allocate data buffer in which to linearize the particle data
  const size_t bufferSize = objSize*totNumParticles;
  char* buffer = new char[bufferSize];

  // are we good to go?
  if (buffer == NULL) {
    MayDay::Error("WritePart::Error: new returned NULL pointer ");
  }
  
  // remember, only one box
  size_t offset = 0;
  
  // write out the particles to the buffer
  for (int ip = 0; ip < totNumParticles; ip++)
    {
      a_particles[ip].linearOut((void*)buffer);
      buffer += objSize;
    }

  // rewind the buffer pointer
  buffer -= bufferSize;

  // write the buffer to the HDF5
  writeDataChunk(offset, dataspace, dataset, H5T_type,
		 bufferSize / sizeof(Real), buffer);
  
  // free buffer
  delete[] buffer; 
  buffer=NULL;
  
  // done
  H5Sclose(dataspace);
  H5Dclose(dataset);
}

#endif // HDF5

#include "NamespaceFooter.H"

#endif
