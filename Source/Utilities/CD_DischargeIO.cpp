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
#include <EBAMRIO.H>
#include <AMRIO.H>
#include <CH_HDF5.H>
#include <PolyGeom.H>

// Our includes
#include <CD_DischargeIO.H>
#include <CD_NamespaceHeader.H>

void DischargeIO::writeEBHDF5(const std::string&                     a_filename,
			      const Vector<DisjointBoxLayout>&       a_grids,
			      const Vector<LevelData<EBCellFAB>* > & a_data,
			      const Vector<std::string>&             a_variableNames,
			      const ProblemDomain&                   a_coarsestDomain,
			      const Real&                            a_coarsestDx,
			      const Real&                            a_dt,
			      const Real&                            a_time,
			      const RealVect&                        a_probLo,			      
			      const Vector<int>&                     a_refinementRatios,
			      const int&                             a_numLevels,
			      const IntVect                          a_ghostVect) {
  CH_TIME("DischargeIO::writeEBHDF5");

  const int numInputVars  = a_data[0]->nComp();  

  // Basic assertions to make sure the input makes sense. 
  CH_assert(a_numLevels               >  0              );
  CH_assert(a_grids.           size() >= a_numLevels    );
  CH_assert(a_data.            size() >= a_numLevels    );
  CH_assert(a_refinementRatios.size() >= a_numLevels - 1);
  CH_assert(a_variablesNames.  size() >= numInputVars   );

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

    // This is the grid resolution on this level.
    Real dx = a_coarsestDx;
    for (int i = 0; i < lvl; i++){
      dx *= a_refinementRatios[lvl];
    }

    // Allocate the Chombo FArrayBox data. 
    chomboData[lvl] = new LevelData<FArrayBox> (dbl, numCompTotal, a_ghostVect);
      
    LevelData<FArrayBox>& fabData = *chomboData[lvl];

    // Go through all the grids on this level
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      // Handle to single-valued data.
      FArrayBox& currentFab = fabData[dit()];

      // Grid information and user input variables.
      const EBCellFAB& data    = (*a_data[lvl])[dit()];
      const EBISBox&   ebisbox = data.getEBISBox();
      const Box        box     = currentFab.box() & ebisbox.getDomain();

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

	  // set distance unless the length of the normal is zero
	  Real length = PolyGeom::dot(normal,normal);

	  if (length > 0){
	    Real dist = PolyGeom::computeAlpha(volFrac,normal) * dx;
	    currentFab(iv,indexDist) = -dist;
	  }
	}
      }
    }
  } 

  // Call the Chombo writer. 
  WriteAMRHierarchyHDF5(a_filename,
                        a_grids,
                        chomboData,
                        variableNamesHDF5,
                        a_coarsestDomain.domainBox(),
                        a_coarsestDx,
                        a_dt,
                        a_time,
                        a_refinementRatios,
                        a_numLevels);

  // Clean up memory.
  for (int lvl = 0; lvl < a_numLevels; lvl++) {
    delete chomboData[lvl];
  }  
}

#include <CD_NamespaceFooter.H>
