/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_LinearStencil.cpp
  @brief  Implementation of cd_LinearStencil.H
  @author Robert Marskar
*/

// Our includes
#include <CD_LinearStencil.H>
#include <CD_NamespaceHeader.H>

#define DEBUG_LINEARSTENCIL 1


  
Real LinearStencil::tolerance = 1.E-6;



bool LinearStencil::getLinearInterpStencil(VoFStencil&          a_stencil, 
					   const RealVect&      a_centroid,
					   const VolIndex&      a_vof,
					   const ProblemDomain& a_domain,
					   const EBISBox&       a_ebisbox){
  bool do_interp     = true;
  bool found_stencil = false;

  if(a_centroid.vectorLength() < LinearStencil::tolerance){
    do_interp = false;
  }

  if(do_interp){
#if CH_SPACEDIM == 2
    found_stencil = LinearStencil::computeInterpStencil2D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox, 99);
#elif CH_SPACEDIM == 3
    found_stencil = LinearStencil::computeInterpStencil3D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox);
#endif
  }
  else{ 
    found_stencil = true;
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);
  }
  
  return found_stencil;
}

bool LinearStencil::computeInterpStencil1D(VoFStencil&          a_stencil,
					   const RealVect&      a_centroid,
					   const VolIndex&      a_vof, 
					   const ProblemDomain& a_domain,
					   const EBISBox&       a_ebisbox,
					   const int            a_interp_dir){
  CH_TIME("LinearStencil::computeInterpStencil1D");

  bool found_stencil = false;
  
  // Switch coordinate direction in order to cover both Hi/Lo interpolation
  const int HiLoInterp = a_centroid[a_interp_dir] > 0. ? 1 : -1;
  const Real x         = a_centroid[a_interp_dir]*HiLoInterp;   // Interpolation distance always positive

  a_stencil.clear();
  const bool really1D = (x > LinearStencil::tolerance) ? true : false;

  if(really1D){
    // Check if vof1 is available
    Side::LoHiSide side;
    if(HiLoInterp == 1){
      side = Side::Hi;
    }
    else if (HiLoInterp == -1){
      side = Side::Lo;
    }

    // Get all Vofs to the correct side
    const Vector<VolIndex> otherVoFs = a_ebisbox.getVoFs(a_vof, a_interp_dir, side, 1);
    const int numOtherVoFs = otherVoFs.size();
    if(numOtherVoFs > 0){ // We can find a stencil
      found_stencil = true;
      a_stencil.clear();
      a_stencil.add(a_vof, 1-x);

      // Average over the other VoFs (other cell might be multivalued)
      for (int i = 0; i < numOtherVoFs; i++){
	a_stencil.add(otherVoFs[i], x/numOtherVoFs);
      }
    }
    else{ // Can't do deriv and can't find stencil
      a_stencil.clear();
      a_stencil.add(a_vof, 1.0);

      found_stencil = false;
    }
  }
  else{
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);
    found_stencil = true;
  }

#if DEBUG_LINEARSTENCIL // Debug hho. Check for negative weights and that the sum equals 1
  Real sumweights = 0.0;
  for (int i = 0; i < a_stencil.size(); i++){
    const Real w = a_stencil.weight(i);
    if(w < 0.0) MayDay::Warning("LinearStencil::computeInterpStencil1d - linear negative weight");
    sumweights += w;
  }
  if(Abs((sumweights - 1.0)) > 1.E-5) MayDay::Warning("LinearStencil::computeInterpStencil1d - sum of weights are not 1");
#endif

  return found_stencil;
}


bool LinearStencil::computeInterpStencil2D(VoFStencil&          a_stencil,
					   const RealVect&      a_centroid,
					   const VolIndex&      a_vof,
					   const ProblemDomain& a_domain,
					   const EBISBox&       a_ebisbox,
					   const int            a_nointerp_dir){
  CH_TIME("LinearStencil::computeInterpStencil2D");
  CH_assert(SpaceDim == 2 || SpaceDim == 3);

  // This code is a bit convoluted because of the many pathological cases that may happen with multi-valued cells. For example
  // we may have situations like this where we want to interpolate point (o) in cell A:
  //
  //    |----------|-----------|
  //    |          |           |
  //    | B multi  | C  reg    |
  //    |          |           |
  //    |----------| ----------|
  //    |        o |           |
  //    | A irreg  | D  reg    |
  //    |          |           |
  //    |----------| ----------|
  //
  // In order to interpolate to the centroid (o) in cell A, which is irregular, we must be able to account for the possibly
  // multiple degrees of freedom in cell B. So, we begin by checking that we can get 1D stencils from A-B and from A-D. If only
  // one of those exist then we must do 1D interpolation and we default to the one stencil that was computed. Since B is a
  // multi-valued cell, the 1D interpolation stencil from A-B uses all the available degrees of freedom. Next, we obtain 
  // all VoFs in cells B and D and compute interpolation stencils to C. If all of those stencils exist (they usually do), then
  // we are able to obtain a reasonable bilinear stencil that also account for multi-valued cells. There is the rare case that
  // we can find A-B and A-D, but not CD or BC. In those cases we _could_ take the interpolation stencil as the average of
  // AB and AD, but that can potentially give negative weights in 3D application so instead we default to the one of the stencil
  // AB or AD; we take the 1D interpolation along the coordinate axis where a_centroid is displaced furthest. If there are several
  // multivalued cells we should be able to handle those too. For example, if D is a multivalued cell then we obtain a list
  // of VoFs in D and compute the average interpolation stencil CD. 

  bool found_stencil = false;

  // Interpolation directions, in 2D it is always 0 (x) and 1 (y)
  int dir0 = 0;
  int dir1 = 1;

  // In 3D we may specify which direction we will interpolate in. No interpolation in z means bilinear interpolation in x,y
  // and no interpolation in x means bilinear interpolation in y,z and so on. 
#if CH_SPACEDIM == 3
  if(a_nointerp_dir == 2){
    dir0 = 0;
    dir1 = 1;
  }
  else if(a_nointerp_dir == 0){
    dir0 = 1;
    dir1 = 2;
  }
  else if(a_nointerp_dir == 1){
    dir0 = 0;
    dir1 = 2;
  }
#endif

  bool really2D = true;
  if(Abs(a_centroid[dir0]) <= LinearStencil::tolerance || Abs(a_centroid[dir1]) < LinearStencil::tolerance){
    really2D = false;
  }


  if(really2D){
    const int HiLoDir0 = a_centroid[dir0] >= 0. ? 1 : -1;
    const int HiLoDir1 = a_centroid[dir1] >= 0. ? 1 : -1;

    const Side::LoHiSide side0 = (HiLoDir0 == 1) ? Side::Hi : Side::Lo;
    const Side::LoHiSide side1 = (HiLoDir1 == 1) ? Side::Hi : Side::Lo;

    // Make sure we have VoFs to all sides
    const Vector<VolIndex> otherVoFsInDir0 = a_ebisbox.getVoFs(a_vof, dir0, side0, 1);
    const Vector<VolIndex> otherVoFsInDir1 = a_ebisbox.getVoFs(a_vof, dir1, side1, 1);

    // 0. Check if we can interpolate from a_vof in dir0 and dir1
    const bool canInterpDir0 = (otherVoFsInDir0.size() > 0);
    const bool canInterpDir1 = (otherVoFsInDir1.size() > 0);

    // If we can only interpolate in one of the direction, return the 1D stencil along that direction. 
    if(!canInterpDir0 && !canInterpDir1){ // This can happen if the EB cuts a domain "corner" and the normal points outwards. 
      a_stencil.clear();
      a_stencil.add(a_vof, 1.0);

      found_stencil = false;
    }
    else if(canInterpDir0 && !canInterpDir1){
      found_stencil = LinearStencil::computeInterpStencil1D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox, dir0);
    }
    else if(!canInterpDir0 &&  canInterpDir1){
      found_stencil = LinearStencil::computeInterpStencil1D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox, dir1);
    }
    else{ // OK, we know we can interpolate in both directions from a_vof. Let's check if we can find a bilinear stencil now
    
      // 1. First, starting on a_vof, compute the linear interpolation stencil in 1D along dir0 and dir1
      VoFStencil curStenDir0;
      VoFStencil curStenDir1;
      const bool foundFirstStencil0 = LinearStencil::computeInterpStencil1D(curStenDir0,a_centroid,a_vof,a_domain,a_ebisbox,dir0);
      const bool foundFirstStencil1 = LinearStencil::computeInterpStencil1D(curStenDir1,a_centroid,a_vof,a_domain,a_ebisbox,dir1);

#if LINEARSTENCIL_DEBUG
      if(!foundFirstStencil0 || !foundFirstStencil1){ // If we made it into this loop, this shouldn't happen
	MayDay::Warning("LinearStencil::computeInterpStencil2D - logic bust");
      }
#endif

      // 2. Check if we can make stencils from the neighboring cells and into the "corner cell". I.e. DC and BC in the figure
      //    above. We actually don't need to do both of them, but we do it anyways. 
      VoFStencil otherVoFsInterpStencilDir0;
      VoFStencil otherVoFsInterpStencilDir1;
    
      int numStencilsAdded0 = 0;
      int numStencilsAdded1 = 0;
    
      bool foundSecondStencil0 = false;
      bool foundSecondStencil1 = false;
    
      for (int i = 0; i<otherVoFsInDir1.size(); i++){
	VoFStencil sten;
	const bool found = LinearStencil::computeInterpStencil1D(sten,a_centroid,otherVoFsInDir1[i],a_domain,a_ebisbox,dir0);
      
	if(found){
	  foundSecondStencil0 = true;
	  otherVoFsInterpStencilDir0 += sten;
	  numStencilsAdded0++;
	}
      }
      for (int i = 0; i<otherVoFsInDir0.size(); i++){
	VoFStencil sten;
	const bool found = LinearStencil::computeInterpStencil1D(sten,a_centroid,otherVoFsInDir0[i],a_domain,a_ebisbox,dir1);
      
	if(found){
	  foundSecondStencil1 = true;
	  otherVoFsInterpStencilDir1 += sten;
	  numStencilsAdded1++;
	}
      }
    
      // 3. Check if we could get the interpolation stencil from the other VoFs. If we could, scale them with the number of VoFs
      //    we accessed. otherVoFsInterpStencilDir0 then contains the interpolation stencil along dir0 from the cell that
      //    was displaced along dir1 from a_vof. E.g. if dir0=x and dir1=y, then otherVoFsInterpStencilDir0 contains the
      //    interpolation stencil B-C in the figure above. 
      if(numStencilsAdded0 > 0) otherVoFsInterpStencilDir0 *= 1./numStencilsAdded0;
      if(numStencilsAdded1 > 0) otherVoFsInterpStencilDir1 *= 1./numStencilsAdded1;
    
#if DEBUG_LINEARSTENCIL // Debug, all stencils should have positive weights with weights summing to 1 at this point. 
      Real s0 = 0.0;
      Real s1 = 0.0;
      for (int i = 0; i < curStenDir0.size(); i++){
	const Real w = curStenDir0.weight(i);
	s0 += w;
	if(w < 0.0) MayDay::Warning("LinearStencil::computeInterpStencil2D - bilinear negative weight in curStenDir0");
      }
      for (int i = 0; i < curStenDir1.size(); i++){
	const Real w = curStenDir1.weight(i);
	s1 += w;
	if(w < 0.0) MayDay::Warning("LinearStencil::computeInterpStencil2D - bilinear negative weight in curStenDir1");
      }
      if((s0 - 1.0) > 1.E-5) MayDay::Warning("LinearStencil::computeInterpStencil2D - curStenDir0 weights do not sum to 1");
      if((s1 - 1.0) > 1.E-5) MayDay::Warning("LinearStencil::computeInterpStencil2D - curStenDir1 weights do not sum to 1");

      s0 = 0.0;
      s1 = 0.0;
      for (int i = 0; i < otherVoFsInterpStencilDir0.size(); i++){
	const Real w = otherVoFsInterpStencilDir0.weight(i);
	s0 += w;
	if(w < 0.0) MayDay::Warning("LinearStencil::computeInterpStencil2D - bilinear negative weight in otherVoFsDir0");
      }
      for (int i = 0; i < otherVoFsInterpStencilDir1.size(); i++){
	const Real w = otherVoFsInterpStencilDir1.weight(i);
	s1 += w;
	if(w < 0.0) MayDay::Warning("LinearStencil::computeInterpStencil2D - bilinear negative weight in otherVoFsDir1");
      }
      if((s0 - 1.0) > 1.E-5) MayDay::Warning("LinearStencil::computeInterpStencil2D - otherVoFsDir0 weights do not sum to 1");
      if((s1 - 1.0) > 1.E-5) MayDay::Warning("LinearStencil::computeInterpStencil2D - otherVoFsDir1 weights do not sum to 1");
#endif

      // 4. We might only be able to find one of the stencils
      a_stencil.clear();
      if(foundFirstStencil0 && !foundFirstStencil1){ // It shouldn't be possible for this to happen
#if DEBUG_LINEARSTENCIL
	MayDay::Warning("LinearStencil::computeInterpStencil2D - logic bust, foundFirstStencil0 && !foundFirstStencil1");
#endif
	a_stencil += curStenDir0;
	found_stencil = true;
      }
      else if(!foundFirstStencil0 && foundFirstStencil1){ // It shouldn't be possible for this to happen
#if DEBUG_LINEARSTENCIL
	MayDay::Warning("LinearStencil::computeInterpStencil2D - logic bust, !foundFirstStencil0 && foundFirstStencil1");
#endif
	a_stencil += curStenDir1;
	found_stencil = true;
      }
      else if(!foundFirstStencil0 && !foundFirstStencil1){ // Might already have covered this case but not sure
#if DEBUG_LINEARSTENCIL
	MayDay::Warning("LinearStencil::computeInterpStencil2D - logic bust, thought this case was covered");
#endif
	a_stencil.add(a_vof, 1.0);
	found_stencil = false;
      }
      else if(foundFirstStencil0 && foundFirstStencil1 && !(foundSecondStencil0 || foundSecondStencil1)){
#if DEBUG_LINEARSTENCIL
	MayDay::Warning("LinearStencil::computeInterpStencil2D - could not find stencil BC or CD. Defaulting to linear stencil");
#endif
	// Here, we could find stencil AB and AD but not BC or CD. This shouldn't happen but if it does, we default to 1D
	// interpolation along the axis with the largest displacement. 
	if(Abs(a_centroid[dir0]) > Abs(a_centroid[dir1])){ // 1D interpolation in the longest direction
	  a_stencil += curStenDir0;
	}
	else{
	  a_stencil += curStenDir1;
	}
	found_stencil = true;
      }
      else{ // Full bilinear interpolation. We should be able to do any combination we like (right?... RIGHT????)

	// Linear interpolation along dir0 followed by interpolation along dir1. The choice of directions SHOULDN'T matter,
	// but because the neighboring cells might be multi-valued and we do an averaging over possible 1D interpolation stencils
	// that may or may not be averages over multi-valued cells, the direction might actually matter. These pathological cases 
	// will be extremely rare, and I don't think they will ever be a problem. But if they do, please check this code. 
	const Real d = a_centroid[dir1]*HiLoDir1;

	const Real w0 = (1-d);
	const Real w1 = d;

	curStenDir0 *= w0;
	otherVoFsInterpStencilDir0 *= w1;
      
	a_stencil += curStenDir0;
	a_stencil += otherVoFsInterpStencilDir0;

	found_stencil = true;
      }
    }
  }
  else{ // Not really 2D code
    const bool interpDir0 = (Abs(a_centroid[dir0]) > LinearStencil::tolerance) ? true : false;
    const bool interpDir1 = (Abs(a_centroid[dir1]) > LinearStencil::tolerance) ? true : false;

    a_stencil.clear();
    if(!interpDir0){// Interp along dir1.
      found_stencil = LinearStencil::computeInterpStencil1D(a_stencil,a_centroid,a_vof,a_domain,a_ebisbox,dir1); 
    }
    else if(!interpDir1){
      found_stencil = LinearStencil::computeInterpStencil1D(a_stencil,a_centroid,a_vof,a_domain,a_ebisbox,dir0); 
    }
    else if(!interpDir0 && !interpDir1){
      a_stencil.add(a_vof, 1.0);
      found_stencil = true;
    }
    else{
      MayDay::Warning("LinearStencil::computeInterpStencil2d - logic bust");
    }
  }

#if DEBUG_LINEARSTENCIL
  Real sumweights = 0.0;
  for (int i = 0; i < a_stencil.size(); i++){
    const Real w = a_stencil.weight(i);
    if(w < 0.0) MayDay::Warning("LinearStencil::computeInterpStencil2D - bilinear negative weight");
    sumweights += w;
  }
  if(Abs((sumweights - 1.0)) > 1.E-5) MayDay::Warning("LinearStencil::computeInterpStencil2D - weights do not sum to 1");
#endif

  return found_stencil;
}

bool LinearStencil::computeInterpStencil3D(VoFStencil&          a_stencil,
					   const RealVect&      a_centroid,
					   const VolIndex&      a_vof,
					   const ProblemDomain& a_domain,
					   const EBISBox&       a_ebisbox){
  CH_TIME("computeInterpStencil3D");
  CH_assert(SpaceDim == 3);
  
  
  bool found_stencil = false;

  const int interpDir3D = 2;


  // Check if centroid lies in one of the coordinate planes connecting cell centers
  const Real tol = LinearStencil::tolerance;
  bool really3D = true;
  if(Abs(a_centroid[0]) <= tol || Abs(a_centroid[1]) <= tol || Abs(a_centroid[2]) <= tol){
    really3D = false;
  }

  // This code is crude as balls, and is only guaranteed to be correct if there aren't multivalued cells. When there are
  // multivalued cells in play, I'm not really sure how to do this because the 
  // code when we get the chance.
  if(really3D){

    //    MayDay::Warning("LinearStencil::computeInterpStencil3D - fix this code");
    const int HiLoInterpDir   = a_centroid[interpDir3D] > 0. ? 1 : -1;
    const Side::LoHiSide side = (a_centroid[interpDir3D] > 0.) ? Side::Hi : Side::Lo;

    // Compute bilinear stencils for two rows of VoFs. stenThisPlane is the bilinear stencil in the plane of this VoF,
    // and stenOtherPlane is the bilinear stencil in the plane of the neighboring vof
#if 0
    const IntVect iv0 = a_vof.gridIndex();
    const IntVect iv1 = iv0 + BASISV(interpDir3D)*HiLoDir2;
    const VolIndex vof0 = VolIndex(iv0,0);
    const VolIndex vof1 = VolIndex(iv1,0);
#endif


    // This is the bilinear stencil in the xy-plane of this VoF. This will always return true (right..?) because
    // of how we handle pathological cases for 2D
    VoFStencil stenThisPlane;
    const bool foundFirstStencil  = LinearStencil::computeInterpStencil2D(stenThisPlane,
									  a_centroid,
									  a_vof,
									  a_domain,
									  a_ebisbox,
									  interpDir3D);


    // This code tries to compute a bilinear stencil in the xy plane but for the cell (i,j,k+-1). There can be multiple if
    // the cell (i,j,k+1) is multivalued, hence the loop. We take the other bilinear stencil to be the average of those
    // stencils. 
    VoFStencil stenOtherPlane;       // This is the (possibly averaged) bilinear interpolation stencil in the other plane
    bool foundSecondStencil = false; // Stupid flag
    int numOtherStencils    = 0;     // Number of stencils that we got (possibly > 1 if multivalued)
    const Vector<VolIndex> otherVoFs = a_ebisbox.getVoFs(a_vof, interpDir3D, side, 1); 
    for (int ivof = 0; ivof < otherVoFs.size(); ivof++){
      VoFStencil otherStencil;
      const bool foundOtherStencil = LinearStencil::computeInterpStencil2D(otherStencil,
									   a_centroid,
									   otherVoFs[ivof],
									   a_domain,
									   a_ebisbox,
									   interpDir3D);

      if(foundOtherStencil){
	foundSecondStencil = true;
	stenOtherPlane += otherStencil;
	numOtherStencils++;
      }
    }
    if(foundSecondStencil){ // Do the average
      stenOtherPlane *= 1./numOtherStencils;
    }
    

    // First hook: If we got both stencils we can interpolate them safely. If this fails, we have to do bilinear interpolation
    if(foundFirstStencil && foundSecondStencil){
      found_stencil = true;
      a_stencil.clear();

      // If we made it here, we have bilinear interpolation stencils in both planes and we can
      // do linear interpolation of those two in order to get the 3D stencil
      const Real z  = a_centroid[2]*sign(side);//HiLoDir2;
      const Real w0 = (1 - z);
      const Real w1 = z;

      stenThisPlane  *= w0;
      stenOtherPlane *= w1;

      a_stencil += stenThisPlane;
      a_stencil += stenOtherPlane;

#if 0 // debugging hook
      std::cout << "vof = " << a_vof.gridIndex() << std::endl;
      std::cout << "stencil size = " << a_stencil.size() << std::endl;
      std::cout << "stencil0 size = " << stenThisPlane.size() << std::endl;
      std::cout << "stencil1 size = " << stenOtherPlane.size() << std::endl;
      for (int i = 0; i < a_stencil.size(); i++){
	//	std::cout << a_stencil.vof(i).gridIndex() << "\t" << a_stencil.weight(i) << std::endl;
      }

      for (int i = 0; i < stenThisPlane.size(); i++){
	std::cout << stenThisPlane.vof(i).gridIndex() << "\t" << stenThisPlane.weight(i) << std::endl;
      }
      for (int i = 0; i < stenOtherPlane.size(); i++){
	std::cout << stenOtherPlane.vof(i).gridIndex() << "\t" << stenOtherPlane.weight(i) << std::endl;
      }
      MayDay::Warning("stop");
#endif
    }
    else if(foundFirstStencil && !foundSecondStencil){ // This does bilinear interpolation instead of trilinear interpolation
      a_stencil.clear();
      a_stencil += stenThisPlane;
      found_stencil = true;
    }
    else if(!foundFirstStencil && !foundSecondStencil){ // Could not find any stencils
      found_stencil = false;
    }
    else{
      MayDay::Warning("LinearStencil::computeInterpStencil3d - logic bust");
    }
  }
  else if(!really3D){ // Find a direction in which we shouldn't interpolate and get a 2D stencil. This might be several
    // directions but that is taken care of in the 2D computation. 
    int noInterpDir = 0;
    for (int dir = 0; dir < SpaceDim; dir++){
      noInterpDir = Abs(a_centroid[dir]) <= tol ? dir : noInterpDir;
    }
    found_stencil = computeInterpStencil2D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox, noInterpDir);
  }

  // This can happen if the EB cuts a domain "corner" and the normal points outwards. This will cause the above code
  // to try to reach out of the domain, but those cells are covered and we don't want to reach into them anyways since they can
  // contain bogus data. 
  if(!found_stencil){
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);

    found_stencil = false;
  }

#if DEBUG_LINEARSTENCIL
  Real sumweights = 0.0;
  for (int i = 0; i < a_stencil.size(); i++){
    const Real w = a_stencil.weight(i);
    if(w < 0.0) MayDay::Warning("LinearStencil::computeInterpStencil3d - trilinear stencil weight < 0.0");
    sumweights += w;
  }
  if(Abs((sumweights - 1.0)) > 1.E-5) MayDay::Warning("LinearStencil::computeInterpStencil3d - sum of weights not equal to one");
#endif
  
  return found_stencil;
}

#include <CD_NamespaceFooter.H>
