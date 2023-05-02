/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DischargeIO.cpp
  @brief  Implementation of DischargeIO.H
  @author Robert Marskar
*/

// Std includes
#include <sstream>

// Chombo includes
#include <CH_HDF5.H>
#include <EBAMRIO.H>
#include <PolyGeom.H>

// Our includes
#include <CD_DischargeIO.H>
#include <CD_NamespaceHeader.H>

std::string
DischargeIO::numberFmt(const long long n, char sep) noexcept
{
  CH_TIME("DischargeIO::numberFmt(long long, char)");

  std::stringstream fmt;
  fmt << n;
  string s = fmt.str();
  s.reserve(s.length() + s.length() / 3);

  for (int i = 0, j = 3 - s.length() % 3; i < s.length(); ++i, ++j)
    if (i != 0 && j % 3 == 0)
      s.insert(i++, 1, sep);

  return s;
}

Vector<std::string>
DischargeIO::numberFmt(const Vector<long long> a_numbers, char a_sep) noexcept
{
  CH_TIME("DischargeIO::numberFmt(Vector<long long>, char)");

  Vector<std::string> ret(a_numbers.size());

  for (int i = 0; i < a_numbers.size(); i++) {
    ret[i] = numberFmt(a_numbers[i], a_sep) + " ";
  }

  return ret;
}

#ifdef CH_USE_HDF5
void
DischargeIO::writeEBHDF5Header(HDF5Handle&                a_handleH5,
                               const int                  a_numLevels,
                               const RealVect&            a_probLo,
                               const Vector<std::string>& a_variableNames) noexcept
{
  CH_TIME("DischargeIO::writeEBHDF5Header");

  CH_assert(a_handleH5.isOpen());
  CH_assert(a_numLevels >= 0);

  const int numInputVars      = a_variableNames.size();
  const int indexVolFrac      = numInputVars;
  const int indexBoundaryArea = indexVolFrac + 1;
  const int indexAreaFrac     = indexBoundaryArea + 1;
  const int indexNormal       = indexAreaFrac + 2 * SpaceDim;
  const int indexDist         = indexNormal + SpaceDim;
  const int numCompTotal      = indexDist + 1;

  // Now create a vector of all the variable names. This is the user input variables plus the EB-related variables
  // for doing the EB reconstruction.
  Vector<std::string> variableNamesHDF5(numCompTotal);

  const std::string volFracName("fraction-0");
  const std::string boundaryAreaName("boundaryArea-0");

  Vector<std::string> areaName(6);
  areaName[0] = "xAreafractionLo-0";
  areaName[1] = "xAreafractionHi-0";
  areaName[2] = "yAreafractionLo-0";
  areaName[3] = "yAreafractionHi-0";
  areaName[4] = "zAreafractionLo-0";
  areaName[5] = "zAreafractionHi-0";

  Vector<std::string> normName(3);
  normName[0] = "xnormal-0";
  normName[1] = "ynormal-0";
  normName[2] = "znormal-0";

  const std::string distName("distance-0");

  // Start appending names.
  for (int i = 0; i < a_variableNames.size(); i++) {
    variableNamesHDF5[i] = a_variableNames[i];
  }

  variableNamesHDF5[indexVolFrac]      = volFracName;
  variableNamesHDF5[indexBoundaryArea] = boundaryAreaName;

  for (int i = 0; i < 2 * SpaceDim; i++) {
    variableNamesHDF5[indexAreaFrac + i] = areaName[i];
  }

  for (int i = 0; i < SpaceDim; i++) {
    variableNamesHDF5[indexNormal + i] = normName[i];
  }

  variableNamesHDF5[indexDist] = distName;

  // Write the header to file.
  HDF5HeaderData header;

  header.m_string["filetype"]    = "VanillaAMRFileType";
  header.m_int["num_levels"]     = a_numLevels;
  header.m_int["num_components"] = numCompTotal;
#if 0 // Uncommenting this because VisIt doesn't know what to do with it.
  header.m_realvect["prob_lo"]   = a_probLo;
#endif

  for (int comp = 0; comp < numCompTotal; comp++) {
    char labelString[100];
    sprintf(labelString, "component_%d", comp);

    std::string label(labelString);

    header.m_string[label] = variableNamesHDF5[comp];
  }

  header.writeToFile(a_handleH5);
}
#endif

#ifdef CH_USE_HDF5
void
DischargeIO::writeEBHDF5Level(HDF5Handle&                 a_handleH5,
                              const LevelData<EBCellFAB>& a_outputData,
                              const ProblemDomain         a_domain,
                              const Real                  a_dx,
                              const Real                  a_dt,
                              const Real                  a_time,
                              const int                   a_level,
                              const int                   a_refRatio,
                              const int                   a_numGhost) noexcept
{
  CH_TIMERS("DischargeIO::writeEBHDF5Level");
  CH_TIMER("DischargeIO::writeEBHDF5Level::alloc", t1);
  CH_TIMER("DischargeIO::writeEBHDF5Level::copy_vars", t2);
  CH_TIMER("DischargeIO::writeEBHDF5Level::average_multicells", t3);
  CH_TIMER("DischargeIO::writeEBHDF5Level::set_default_data", t4);
  CH_TIMER("DischargeIO::writeEBHDF5Level::set_eb_moments", t5);
  CH_TIMER("DischargeIO::writeEBHDF5Level::set_ghosts", t6);

  CH_assert(a_refRatio > 0);
  CH_assert(a_numGhost >= 0);

  const int numInputVars      = a_outputData.nComp();
  const int indexVolFrac      = numInputVars;
  const int indexBoundaryArea = indexVolFrac + 1;
  const int indexAreaFrac     = indexBoundaryArea + 1;
  const int indexNormal       = indexAreaFrac + 2 * SpaceDim;
  const int indexDist         = indexNormal + SpaceDim;
  const int numCompTotal      = indexDist + 1;

  const DisjointBoxLayout& dbl = a_outputData.disjointBoxLayout();

  CH_START(t1);
  LevelData<FArrayBox> levelData(dbl, numCompTotal, a_numGhost * IntVect::Unit);
  CH_STOP(t1);

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    FArrayBox&       levelFAB      = levelData[dit()];
    const EBCellFAB& outputData    = a_outputData[dit()];
    const FArrayBox& outputDataReg = outputData.getFArrayBox();
    const EBISBox&   ebisbox       = outputData.getEBISBox();
    const Box        outputBox     = levelFAB.box() & a_domain;

    // Copy the regular data. This will also do single-valued cut-cells. Later, we average the
    // contents in multi-valued cells.
    CH_START(t2);
    levelFAB.setVal(0.);
    levelFAB.copy(outputDataReg, 0, 0, numInputVars);
    CH_STOP(t2);

    // Run through the multi-valued cells and set the single-valued output-data to be the sum of the multi-valued data. I just
    // don't know of a different way of doing this that would also play well with HDF5..
    CH_START(t3);
    const IntVectSet multiCells = ebisbox.getMultiCells(outputBox);
    for (IVSIterator ivsIt(multiCells); ivsIt.ok(); ++ivsIt) {
      const IntVect iv = ivsIt();

      const Vector<VolIndex>& multiVofs = ebisbox.getVoFs(iv);
      const int               numVoFs   = ebisbox.numVoFs(iv);

      CH_assert(numVoFs > 0);

      const Real invNumVoFs = 1. / (Real(numVoFs));

      for (int comp = 0; comp < numInputVars; ++comp) {
        levelFAB(iv, comp) = 0.0;

        for (int ivof = 0; ivof < numVoFs; ivof++) {
          levelFAB(iv, comp) += outputData(multiVofs[ivof], comp);
        }

        levelFAB(iv, comp) *= invNumVoFs;
      }
    }
    CH_STOP(t3);

    CH_START(t4);
    // Set default volume fraction, area fraction, normal, and distance of EB from corner.
    levelFAB.setVal(1.0, indexVolFrac);
    levelFAB.setVal(0.0, indexBoundaryArea);
    levelFAB.setVal(0.0, indexDist);
    for (int i = 0; i < 2 * SpaceDim; i++) {
      levelFAB.setVal(1.0, indexAreaFrac + i);
    }

    for (int i = 0; i < SpaceDim; i++) {
      levelFAB.setVal(0.0, indexNormal + i);
    }
    CH_STOP(t4);

    // Set EB moment data for each cell.
    CH_START(t5);
    for (BoxIterator bit(outputBox); bit.ok(); ++bit) {
      const IntVect& iv = bit();

      const bool isCovered   = ebisbox.isCovered(iv);
      const bool isRegular   = ebisbox.isRegular(iv);
      const bool isIrregular = !isCovered && !isRegular;

      // In covered cells the volume and area fractions must be zero.
      if (isCovered) {
        levelFAB(iv, indexVolFrac) = 0.0;

        for (int i = 0; i < 2 * SpaceDim; i++) {
          levelFAB(iv, indexAreaFrac + i) = 0.0;
        }
      }
      else if (isIrregular) {
        const Vector<VolIndex> vofs = ebisbox.getVoFs(iv);

        // Take the last vof and put it in the data. For multi-cells this doesn't really make a whole lot of sense...
        const VolIndex& vof0 = vofs[vofs.size() - 1];

        Real     volFrac   = ebisbox.volFrac(vof0);
        Real     bndryArea = ebisbox.bndryArea(vof0);
        RealVect normal    = ebisbox.normal(vof0);

        if (bndryArea == 0.0) {
          if (volFrac > 0.5) {
            volFrac = 1.0;
          }
          else {
            volFrac = 0.0;
          }

          normal = RealVect::Zero;
        }

        // set volume fraction, EB boundary area, and face area fractions.
        levelFAB(iv, indexVolFrac)      = volFrac;
        levelFAB(iv, indexBoundaryArea) = bndryArea;

        for (int dir = 0; dir < SpaceDim; dir++) {
          const Vector<FaceIndex> facesLo = ebisbox.getFaces(vof0, dir, Side::Lo);
          const Vector<FaceIndex> facesHi = ebisbox.getFaces(vof0, dir, Side::Hi);

          if (facesLo.size() == 0) {
            levelFAB(iv, indexAreaFrac + 2 * dir) = 0.0;
          }
          else {
            levelFAB(iv, indexAreaFrac + 2 * dir) = ebisbox.areaFrac(facesLo[0]);
          }

          if (facesHi.size() == 0) {
            levelFAB(iv, indexAreaFrac + 2 * dir + 1) = 0.0;
          }
          else {
            levelFAB(iv, indexAreaFrac + 2 * dir + 1) = ebisbox.areaFrac(facesHi[0]);
          }
        }

        // Set the EB normal vector.
        for (int dir = 0; dir < SpaceDim; dir++) {
          levelFAB(iv, indexNormal + dir) = normal[dir];
        }

        // Set distance unless the length of the normal is zero.
        if (normal.vectorLength() > 0.0) {
          levelFAB(iv, indexDist) = -PolyGeom::computeAlpha(volFrac, normal) * a_dx;
        }
      }
    }
    CH_STOP(t5);

    // At this point we want to fill one layer of ghost cells OUTSIDE the domain.
    CH_START(t6);
    if (a_numGhost > 0) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        for (SideIterator sit; sit.ok(); ++sit) {
          const IntVect shift = sign(flip(sit())) * BASISV(dir); // => +1 for Lo side and -1 for high side.

          // Get the layer of cells immediately outside this box.
          const Box domainBox   = a_domain.domainBox();
          const Box validBox    = outputBox & domainBox;
          const Box boundaryBox = adjCellBox(validBox, dir, sit(), 1);

          if (!(domainBox.contains(boundaryBox))) {
            for (BoxIterator bit(boundaryBox); bit.ok(); ++bit) {
              const IntVect iv = bit();

              for (int comp = 0; comp < numCompTotal; comp++) {
                levelFAB(iv, comp) = levelFAB(iv + shift, comp);
              }
            }
          }
        }
      }
    }
    CH_STOP(t6);
  }

  const int success = writeLevel(a_handleH5,
                                 a_level,
                                 levelData,
                                 a_dx,
                                 a_dt,
                                 a_time,
                                 a_domain.domainBox(),
                                 a_refRatio,
                                 a_numGhost * IntVect::Unit,
                                 Interval(0, numCompTotal - 1));

  if (success != 0) {
    MayDay::Error("DischargeIO::writeEBHDF5 -- error in writeLevel");
  }
}
#endif

#ifdef CH_USE_HDF5
void
DischargeIO::writeEBHDF5(const std::string&                   a_filename,
                         const Vector<std::string>&           a_variableNames,
                         const Vector<DisjointBoxLayout>&     a_grids,
                         const Vector<LevelData<EBCellFAB>*>& a_data,
                         const Vector<ProblemDomain>&         a_domains,
                         const Vector<Real>                   a_dx,
                         const Vector<int>                    a_refinementRatios,
                         const Real                           a_dt,
                         const Real                           a_time,
                         const RealVect                       a_probLo,
                         const int                            a_numLevels,
                         const int                            a_numGhost) noexcept
{
  CH_TIME("DischargeIO::writeEBHDF5");

  // Basic assertions to make sure the input makes sense.
  CH_assert(a_numLevels > 0);
  CH_assert(a_numGhost >= 0);
  CH_assert(a_grids.size() >= a_numLevels);
  CH_assert(a_data.size() >= a_numLevels);
  CH_assert(a_refinementRatios.size() >= a_numLevels - 1);
  CH_assert(a_variableNames.size() >= numInputVars);

  // Indices for where we store the Chombo stuff. This is for storing the volume fraction, EB boundary area,
  // face area fractions etc. These things are written AFTER the user input variables.
  const int numInputVars      = a_data[0]->nComp();
  const int indexVolFrac      = numInputVars;
  const int indexBoundaryArea = indexVolFrac + 1;
  const int indexAreaFrac     = indexBoundaryArea + 1;
  const int indexNormal       = indexAreaFrac + 2 * SpaceDim;
  const int indexDist         = indexNormal + SpaceDim;
  const int numCompTotal      = indexDist + 1;

  // Write header.
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  HDF5Handle handle(a_filename.c_str(), HDF5Handle::CREATE);

  // write header data
  DischargeIO::writeEBHDF5Header(handle, a_numLevels, a_probLo, a_variableNames);

  // Write levfle-by-level
  for (int lvl = 0; lvl < a_numLevels; lvl++) {
    const int refRat = (lvl < a_numLevels) ? a_refinementRatios[lvl] : 1;
    DischargeIO::writeEBHDF5Level(handle,
                                  *a_data[lvl],
                                  a_domains[lvl],
                                  a_dx[lvl],
                                  a_dt,
                                  a_time,
                                  lvl,
                                  refRat,
                                  a_numGhost);
  }

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  handle.close();
}
#endif

#ifdef CH_USE_HDF5
void
DischargeIO::writeEBHDF5(const EBAMRCellData& a_data, const std::string& a_file)
{
  CH_TIME("DischargeIO::writeEBHDF5(debug)");

  Vector<LevelData<EBCellFAB>*> raw(a_data.size());

  for (int lvl = 0; lvl < a_data.size(); lvl++) {
    raw[lvl] = &(*a_data[lvl]);
  }

  writeEBAMRname(&raw, a_file.c_str());
}
#endif

#include <CD_NamespaceFooter.H>
