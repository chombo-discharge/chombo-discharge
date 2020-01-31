#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// bvs  06/06/06

#ifndef _MAPPEDDOMAINIO_H_
#define _MAPPEDDOMAINIO_H_


#include <iostream>
#include "Box.H"
#include "Misc.H"
#include "SPACE.H"
#include <map>
#include <string>

#include "NamespaceHeader.H"

class MappedDomain;

// ---------------------------------------------------------
/// Class for reading and writing MappedDomain objects.
class MappedDomainIO
{

public:

  ///
  MappedDomainIO();

  ///
  void writeUnmappedDomainPlotfile(const MappedDomain& a_domain,
                                   int a_coarseningFactor) const ;

  ///
  void writeDomainPlotfile(const MappedDomain& a_domain,
                           int a_coarseningFactor) const ;

  ///
  void read(MappedDomain& a_domain) const;


  std::string filename;

};

// ---------------------------------------------------------

class HDF5Handle;
class DisjointBoxLayout;

void writeMap(HDF5Handle& a_handle,
              int a_level,
              const DisjointBoxLayout& a_layout,
              const MappedDomain& a_map);



#include "NamespaceFooter.H"

#endif // include guard
