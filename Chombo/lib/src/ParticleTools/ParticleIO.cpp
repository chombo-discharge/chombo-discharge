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

#include "ParticleIO.H"

#include "NamespaceHeader.H"

#ifdef CH_USE_HDF5

// write chunk of data and upgrade offset
void writeDataChunk(size_t&             offset,
                    const hid_t&        dataspace,
                    const hid_t&        dataset,
                    const hid_t&        H5T_type,
                    const unsigned long dataLength,
                    const void* const   data)
{
  // timer
  CH_TIME("writeDataChunk");

  hsize_t hyperslabSize(dataLength);
  hsize_t hyperslabOffset(offset);
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
                      &hyperslabOffset, NULL,
                      &hyperslabSize, NULL);

  hid_t hyperslabSpace = H5Screate_simple(1, &hyperslabSize, NULL);
  int status = H5Dwrite(dataset, H5T_type,
                        hyperslabSpace, dataspace,
                        H5P_DEFAULT, data);

  // close
  H5Sclose(hyperslabSpace);

  // upgrade offset
  offset += hyperslabSize;

  if (status < 0)
    {
      MayDay::Error("WriteDataChunk::Error: H5Dwrite returned negative value");
    }
}

// write chunk of data and upgrade offset
void readDataChunk(size_t&             offset,
                   const hid_t&        dataspace,
                   const hid_t&        dataset,
                   const hid_t&        H5T_type,
                   const unsigned long dataLength,
                   void* const         data,
                   const size_t        stride,
                   const size_t        block)
{
  // timer
  CH_TIME("readDataChunk");

  hsize_t hsize(dataLength);
  hsize_t hoffset(offset);
  hsize_t hstride(stride);
  hsize_t hblock(block);

  herr_t herr = H5Sselect_hyperslab(dataspace,
                                    H5S_SELECT_SET,
                                    &hoffset,
                                    &hstride,
                                    &hsize,
                                    &hblock);
  if (herr < 0)
    {
      MayDay::Error("ReadDataChunk::Error: H5Sselect_hyperslab returned negative value");
    }

  hsize_t hcsize(hsize*hblock);
  //  hsize_t hcsize(hsize*(hstride+hblock));
  hid_t hspace = H5Screate_simple(1, &hcsize, NULL);

  herr = H5Dread(dataset,
                 H5T_type,
                 hspace,
                 dataspace,
                 H5P_DEFAULT,
                 data);

  // close
  H5Sclose(hspace);

  // upgrade offset
  offset += hsize*stride;

  if (herr < 0)
    {
      MayDay::Error("ReadDataChunk::Error: H5Dread returned negative value");
    }
}

// write header
void write_hdf_part_header(HDF5Handle&                       a_handle,
                           const BoxLayout&                  a_grids,
                           const vector<unsigned long long>& a_partPerBox,
                           const std::string&                a_dataType)
{
  // save part boxes
  const int gridStatus = write(a_handle,a_grids,a_dataType+":boxes");
  CH_assert(gridStatus == 0);

  // store part-per-box info
  write_vect_to_header(a_handle,a_partPerBox,H5T_NATIVE_ULLONG,a_dataType+":offsets");

  // store pids info
  vector<int> pids;
  for (LayoutIterator li=a_grids.layoutIterator(); li.ok(); ++li)
    {
      pids.push_back(a_grids.procID(li()));
    }
  write_vect_to_header(a_handle,pids,H5T_NATIVE_INT,a_dataType+":processors");
}

//
void read_hdf_part_header(HDF5Handle&                 a_handle,
                          Vector<Box>&                a_grids,
                          vector<unsigned long long>& a_partPerBox,
                          const std::string&          a_dataType,
                          const std::string&          a_path)
{
  // hdf5 group label for output file
  a_handle.setGroup(a_path);

  // Get the grids
  const int gridStatus = read(a_handle,a_grids,a_dataType+":boxes");

  if (gridStatus != 0)
    {
      MayDay::Error("compare_hdf_part: file A does not contain a Vector<Box>");
    }

  // read in offset info
  a_partPerBox.resize(a_grids.size(), 0);
  read_vect_from_header(a_handle,a_partPerBox,H5T_NATIVE_ULLONG,a_dataType+":offsets");
}

#endif // HDF5

#include "NamespaceFooter.H"
