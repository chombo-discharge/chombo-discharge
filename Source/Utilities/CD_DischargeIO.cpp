/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DischargeIO.cpp
  @brief  Implementation of DischargeIO.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_HDF5.H>
#include <EBAMRIO.H>
#include <PolyGeom.H>

// Our includes
#include <CD_DischargeIO.H>
#include <CD_NamespaceHeader.H>

#ifdef CH_USE_HDF5
void DischargeIO::writeEBHDF5(const std::string&                     a_filename,
			      const Vector<std::string>&             a_variableNames,			      
			      const Vector<DisjointBoxLayout>&       a_grids,
			      const Vector<LevelData<EBCellFAB>* > & a_data,
			      const Vector<ProblemDomain>&           a_domains,
			      const Vector<Real>                     a_dx,
			      const Vector<int>                      a_refinementRatios,
			      const Real                             a_dt,
			      const Real                             a_time,
			      const RealVect                         a_probLo,
			      const int                              a_numLevels,
			      const int                              a_numGhost) {
  CH_TIME("DischargeIO::writeEBHDF5");

  const int numInputVars  = a_data[0]->nComp();  

  // Basic assertions to make sure the input makes sense. 
  CH_assert(a_numLevels               >  0              );
  CH_assert(a_numGhost                >= 0              );
  CH_assert(a_grids.           size() >= a_numLevels    );
  CH_assert(a_data.            size() >= a_numLevels    );
  CH_assert(a_refinementRatios.size() >= a_numLevels - 1);
  CH_assert(a_variableNames.   size() >= numInputVars   );

  // Indices for where we store the Chombo stuff. This is for storing the volume fraction, EB boundary area,
  // face area fractions etc. These things are written AFTER the user input variables. 
  const int indexVolFrac      = numInputVars;
  const int indexBoundaryArea = indexVolFrac+1;
  const int indexAreaFrac     = indexBoundaryArea+1;
  const int indexNormal       = indexAreaFrac+2*SpaceDim;
  const int indexDist         = indexNormal+SpaceDim;

  // Total number of variables -- the sum of the user input variables and the EB variables. 
  const int numCompTotal      = indexDist + 1;

  // Now create a vector of all the variable names. This is the user input variables plus the EB-related variables
  // for doing the EB reconstruction. 
  Vector<std::string> variableNamesHDF5(numCompTotal);

  const std::string volFracName     ("fraction-0"    );
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
  for (int i = 0; i < a_variableNames.size(); i++){
    variableNamesHDF5[i] = a_variableNames[i];
  }  

  variableNamesHDF5[indexVolFrac     ] = volFracName;
  variableNamesHDF5[indexBoundaryArea] = boundaryAreaName;

  for (int i = 0; i < 2*SpaceDim; i++) {
    variableNamesHDF5[indexAreaFrac+i] = areaName[i];
  }

  for (int i = 0; i < SpaceDim; i++){
    variableNamesHDF5[indexNormal+i] = normName[i];
  }

  variableNamesHDF5[indexDist] = distName;

  // Done with variables. Now set up some storage for Chombo data on each level. 
  Vector<LevelData<FArrayBox>* > chomboData(a_numLevels, nullptr);

  // set things up for each level
  for (int lvl = 0; lvl < a_numLevels; lvl++){
    const DisjointBoxLayout& dbl  =  a_grids[lvl];

    // This is the domain and grid resolution on this level.
    const Real          dx     = a_dx     [lvl];
    const ProblemDomain domain = a_domains[lvl];

    // Allocate the Chombo FArrayBox data. 
    chomboData[lvl] = new LevelData<FArrayBox> (dbl, numCompTotal, a_numGhost*IntVect::Unit);
      
    LevelData<FArrayBox>& fabData = *chomboData[lvl];

    // Go through all the grids on this level
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      // Handle to single-valued data.
      FArrayBox& currentFab = fabData[dit()];

      // Grid information and user input variables.
      const EBCellFAB& data    = (*a_data[lvl])[dit()];
      const EBISBox&   ebisbox = data.getEBISBox();
      const Box        box     = currentFab.box() & domain;

      // Copy the regular data. This will also do single-valued cut-cells. 
      constexpr int srcStartComp = 0;
      constexpr int dstStartComp = 0;
      
      currentFab.setVal(0.);      
      currentFab.copy(data.getSingleValuedFAB(), srcStartComp, dstStartComp, numInputVars);

      // Run through the multi-valued cells and set the single-valued output-data to be the sum of the multi-valued data. I just
      // don't know of a different way of doing this that would also play well with HDF5..
      const IntVectSet multiCells = ebisbox.getMultiCells(currentFab.box());
      for (IVSIterator ivsIt(multiCells); ivsIt.ok(); ++ivsIt){
	const IntVect iv = ivsIt();

	const Vector<VolIndex>& multiVofs = ebisbox.getVoFs(iv);
	const int               numVoFs   = ebisbox.numVoFs(iv);

	const Real invNumVoFs = 1./(Real(numVoFs));

	// Set to average
	for (int comp = 0; comp < numInputVars; ++comp){
	  currentFab(iv, comp) = 0.0;

	  for (int ivof = 0; ivof < numVoFs; ivof++){
	    currentFab(iv, comp) += data(multiVofs[ivof], comp);
	  }

	  currentFab(iv, comp) *= invNumVoFs;
	}
      }

      // Set default volume fraction, area fraction, normal, and distance of EB from corner. 
      currentFab.setVal(1.0, indexVolFrac     );
      currentFab.setVal(0.0, indexBoundaryArea);
      currentFab.setVal(0.0, indexDist        );
	  
      for (int i = 0; i < 2*SpaceDim; i++){
	currentFab.setVal(1.0, indexAreaFrac+i);
      }

      for (int i = 0; i < SpaceDim; i++){
	currentFab.setVal(0.0, indexNormal+i);
      }

      // Run through the grid box and set the grid information. 
      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect& iv = bit();

	const bool isCovered   = ebisbox.isCovered(iv);
	const bool isRegular   = ebisbox.isRegular(iv);
	const bool isIrregular = !isCovered && !isRegular;

	// In covered cells the volume and area fractions must be zero. 
	if (isCovered){
	  currentFab(iv, indexVolFrac) = 0.0;

	  for (int i = 0; i < 2*SpaceDim; i++){
	    currentFab(iv, indexAreaFrac+i) = 0.0;
	  }
	}
	else if(isIrregular){
	  const Vector<VolIndex> vofs = ebisbox.getVoFs(iv);

	  // Take the last vof and put it in the data. For multi-cells this doesn't really make a whole lot of sense...
	  const VolIndex& vof0 = vofs[vofs.size()-1];
	  
	  Real     volFrac   = ebisbox.volFrac  (vof0);
	  Real     bndryArea = ebisbox.bndryArea(vof0);	  
	  RealVect normal    = ebisbox.normal   (vof0);
	  
	  if (bndryArea == 0.0){
	    if (volFrac > 0.5){
	      volFrac = 1.0;
	    }
	    else{
	      volFrac = 0.0;
	    }

	    normal = RealVect::Zero;
	  }

	  // set volume fraction, EB boundary area, and face area fractions. 
	  currentFab(iv, indexVolFrac     ) = volFrac  ;
	  currentFab(iv, indexBoundaryArea) = bndryArea;

	  for (int dir = 0; dir < SpaceDim; dir++){
	    const Vector<FaceIndex> facesLo = ebisbox.getFaces(vof0, dir, Side::Lo);
	    const Vector<FaceIndex> facesHi = ebisbox.getFaces(vof0, dir, Side::Hi);	    
	    
	    if (facesLo.size() == 0){
	      currentFab(iv, indexAreaFrac + 2*dir) = 0.0;
	    }
	    else{
	      currentFab(iv, indexAreaFrac + 2*dir) = ebisbox.areaFrac(facesLo[0]);
	    }

	    if (facesHi.size() == 0){
	      currentFab(iv, indexAreaFrac + 2*dir + 1) = 0.0;
	    }
	    else{
	      currentFab(iv, indexAreaFrac + 2*dir + 1) = ebisbox.areaFrac(facesHi[0]);
	    }
	  }

	  // Set the EB normal vector. 
	  for (int dir = 0; dir < SpaceDim; dir++){
	    currentFab(iv, indexNormal + dir) = normal[dir];
	  }

	  // Set distance unless the length of the normal is zero.
	  if (normal.vectorLength() > 0.0){
	    currentFab(iv, indexDist) = -PolyGeom::computeAlpha(volFrac,normal) * dx;
	  }
	}
      }

      // At this point we want to fill one layer of ghost cells OUTSIDE the domain.
      if(a_numGhost > 0){
	for (int dir = 0; dir < SpaceDim; dir++){
	  for (SideIterator sit; sit.ok(); ++sit){
	    const IntVect shift = sign(flip(sit())) * BASISV(dir); // => +1 for Lo side and -1 for high side. 

	    // Get the layer of cells immediately outside this box.
	    const Box domainBox   = a_domains[lvl].domainBox();	    
	    const Box validBox    = box & domainBox;
	    const Box boundaryBox = adjCellBox(validBox, dir, sit(), 1);
	    
	    if(!(domainBox.contains(boundaryBox))){
	      for (BoxIterator bit(boundaryBox); bit.ok(); ++bit){
		const IntVect iv = bit();

		for (int comp = 0; comp < numCompTotal; comp++){
		  currentFab(iv, comp) = currentFab(iv + shift, comp);
		}
	      }
	    }
	  }
	}
      }
    } // End grid patch loop.
  } // End of level loop.

  // Now write the data to HDF5.
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  HDF5Handle handle(a_filename.c_str(), HDF5Handle::CREATE);

  // Write the header to file. 
  HDF5HeaderData header;

  header.m_string  ["filetype"]       = "VanillaAMRFileType";
  header.m_int     ["num_levels"]     = a_numLevels;
  header.m_int     ["num_components"] = numCompTotal;
#if 0 // Uncommenting this because although VisIt uses the attribute, slice operators tend to break. 
  header.m_realvect["prob_lo"]        = a_probLo;
#endif

  for (int comp = 0; comp < numCompTotal; comp++){
    char labelString[100];
    sprintf(labelString, "component_%d", comp);

    std::string label(labelString);

    header.m_string[label] = variableNamesHDF5[comp];
  }

  header.writeToFile(handle);

  // Go through each grid level and write it to file. 
  for (int lvl = 0; lvl < a_numLevels; lvl++){
    int refRatio = 1;

    if(lvl != a_numLevels -1){
      refRatio = a_refinementRatios[lvl];
    }

    const int success = writeLevel(handle,
				   lvl,
				   *chomboData[lvl],
				   a_dx[lvl],
				   a_dt,
				   a_time,
				   a_domains[lvl].domainBox(),
				   refRatio,
				   a_numGhost*IntVect::Unit,
				   Interval(0, numCompTotal-1));

    if(success != 0){
      MayDay::Error("DischargeIO::writeEBHDF5 -- error in writeLevel");
    }
  }


#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif  
  handle.close();

  // Clean up memory.
  for (int lvl = 0; lvl < a_numLevels; lvl++) {
    delete chomboData[lvl];
  }  
}
#endif

#ifdef CH_USE_HDF5
void DischargeIO::writeEBHDF5(const EBAMRCellData& a_data,
			      const std::string&   a_file) {
  CH_TIME("DischargeIO::writeEBHDF5(debug)");

  Vector<LevelData<EBCellFAB>* > raw(a_data.size());

  for (int lvl = 0; lvl < a_data.size(); lvl++){
    raw[lvl] = &(*a_data[lvl]);
  }

  writeEBAMRname(&raw, a_file.c_str());
}
#endif

#include <CD_NamespaceFooter.H>
