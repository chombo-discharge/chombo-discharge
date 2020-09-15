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

#ifndef _PARTICLEIO_H_
#define _PARTICLEIO_H_

#include <list>
using std::list;

#include "REAL.H"
#include "Vector.H"
#include "LayoutIterator.H"
#include "CH_HDF5.H"
#include "SPMD.H"
#include "ListBox.H"
#include "ParticleData.H"

#include "NamespaceHeader.H"

#ifdef CH_USE_HDF5

/// Write chunk of data and upgrade offset. Not meant to be a user-facing function.
void writeDataChunk(size_t&             offset,
                    const hid_t&        dataspace,
                    const hid_t&        dataset,
                    const hid_t&        H5T_type,
                    const unsigned long dataLength,
                    const void* const   data);

/// Read chunk of data and upgrade offset. Not meant to be a user-facing function.
void readDataChunk(size_t&             offset,
                   const hid_t&        dataspace,
                   const hid_t&        dataset,
                   const hid_t&        H5T_type,
                   const unsigned long dataLength,
                   void* const         data,
                   const size_t        stride=1,
                   const size_t        block=1);

/// Write the particle header to the HDF5 file. Not meant to be a user-facing function.
void write_hdf_part_header(HDF5Handle&                       a_handle,
                           const BoxLayout&                  a_grids,
                           const vector<unsigned long long>& a_partPerBox,
                           const std::string&                a_dataType);

/// Read the particle header from the HDF5 file. Not meant to be a user-facing function.
void read_hdf_part_header(HDF5Handle&                 a_handle,
                          Vector<Box>&                a_grids,
                          vector<unsigned long long>& a_particlesPerBox,
                          const std::string&          a_dataType,
                          const std::string&          a_path);

/// A helper function that writes a vector of 'T's to the HDF5 file under a_dataname
template <class T>
void write_vect_to_header(HDF5Handle&         a_handle,
                          const vector<T>&    a_vect,
                          const hid_t&        H5T_type,
                          const std::string&  a_dataname);

/// A helper function that reads a vector of 'T's to the HDF5 file under a_dataname
template <class T>
void read_vect_from_header(HDF5Handle&        a_handle,
			   vector<T>&          a_vect,
			   const hid_t&        H5T_type,
			   const std::string&  a_dataname);

/// Write the particle data contained in a_particles to the HDF5 file described by a_handle.
/// This is the function that users should actually call. The string a_dataType gives a way
/// to label the particles, in case multiple particle types are present.
template<class P> void
writeParticlesToHDF(HDF5Handle&            a_handle,
                    const ParticleData<P>& a_particles,
                    const std::string&     a_dataType);

/// Read the particle data contained in the HDF5 file described by a_handle into a_particles.
/// This is the function that users should actually call. The string a_dataType gives a way
/// to label the particles, in case multiple particle types are present.
template<class P>
void readParticlesFromHDF(HDF5Handle&        a_handle,
                          ParticleData<P>&   a_particles,
                          const std::string& a_dataType);

/// Like the above, but takes a Vector<P> instead of a ParticleData.
template<class P> void
writeParticlesToHDF(HDF5Handle&            a_handle,
                    const ParticleData<P>& a_particles,
		    const Box&             a_domain,
                    const std::string&     a_dataType);

#endif // HDF5

#include "NamespaceFooter.H"

// implementation
#include "ParticleIOI.H"

#endif
