/*!
  @file stencil_ops.cpp
  @brief Implementation of stencil_ops.H
  @author Robert Marskar
  @date Jan. 2018
*/

#include "stencil_ops.H"

Real stencil_ops::tolerance = 1.E-8;

bool stencil_ops::get_linear_interp_stencil(VoFStencil&          a_stencil, 
					    const RealVect&      a_centroid,
					    const VolIndex&      a_vof,
					    const ProblemDomain& a_domain,
					    const EBISBox&       a_ebisbox){
  bool do_interp     = true;
  bool found_stencil = false;

  if(a_centroid.vectorLength() < stencil_ops::tolerance){
    //    do_interp = false;
  }

  if(do_interp){
#if CH_SPACEDIM == 2
    found_stencil = stencil_ops::compute_bilinear_stencil(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox);
    //    found_stencil = compute_interp_stencil_2D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox);
#elif CH_SPACEDIM == 3
    //    found_stencil = stencil_ops::compute_interp_stencil_3D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox);
    found_stencil = stencil_ops::compute_interp_stencil_3D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox);
#endif
  }
  else{ 
    found_stencil = false;
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);
  }
  
  return found_stencil;
}

bool stencil_ops::compute_bilinear_stencil(VoFStencil&          a_stencil,
					   const RealVect&      a_centroid,
					   const VolIndex&      a_vof,
					   const ProblemDomain& a_domain,
					   const EBISBox&       a_ebisbox){

  const int comp = 0;
  
  bool found_stencil = false;

  // Get the quadrant and all the cells that are involved
  IntVect quadrant;
  for (int dir = 0; dir < SpaceDim; dir++){
    if(a_centroid[dir] > 0.){
      quadrant[dir] = 1;
    }
    else{
      quadrant[dir] = -1;
    }
  }

  RealVect pos;
  IntVect iv00, iv01, iv10, iv11;

  // Relative positions
  pos  = a_centroid;
  iv00 = a_vof.gridIndex();
  for (int dir = 0; dir < SpaceDim; dir++){
    if(quadrant[dir] == -1){
      pos[dir]  = 1 - pos[dir];
      iv00[dir] = iv00[dir] - 1;
    }
  }
  iv10 = iv00 + BASISV(0);
  iv01 = iv00 + BASISV(1);
  iv11 = iv00 + IntVect::Unit;

  const VolIndex vof00(iv00, comp);
  const VolIndex vof01(iv01, comp);
  const VolIndex vof10(iv10, comp);
  const VolIndex vof11(iv11, comp);

  // Add vofs to stencil. Exit if any is covered.
  const Real x = pos[0];
  const Real y = pos[1];
  
  a_stencil.clear();
  if(!a_domain.contains(iv00)  || !a_domain.contains(iv10)  || !a_domain.contains(iv01)  || !a_domain.contains(iv11)){
    found_stencil = false;
  }
  else if(a_domain.contains(iv00)){
    if(a_ebisbox.isCovered(iv00)){
      found_stencil = false;
    }
  }
  else if(a_domain.contains(iv10)){
    if(a_ebisbox.isCovered(iv10)){
      found_stencil = false;
    }
  }
  else if(a_domain.contains(iv01)){
    if(a_ebisbox.isCovered(iv01)){
      found_stencil = false;
    }
  }
  else if(a_domain.contains(iv11)){
    if(a_ebisbox.isCovered(iv11)){
      found_stencil = false;
    }
  }
  else{
    a_stencil.add(vof00, (1-x)*(1-y));
    a_stencil.add(vof10, x*(1-y));
    a_stencil.add(vof01, (1-x)*y);
    a_stencil.add(vof11, x*y);

    found_stencil = true;
  }

  return found_stencil;
}

bool stencil_ops::compute_interp_stencil_1D(VoFStencil&          a_stencil,
					    const RealVect&      a_centroid,
					    const VolIndex&      a_vof, 
					    const ProblemDomain& a_domain,
					    const EBISBox&       a_ebisbox,
					    const int            a_interp_dir){
  CH_TIME("compute_interp_stencil1D");

  bool found_stencil = false;
  
  // Switch coordinate direction in order to cover both Hi/Lo interpolation
  const int HiLoInterp = a_centroid[a_interp_dir] > 0. ? 1 : -1;
  const Real x         = a_centroid[a_interp_dir]*HiLoInterp;   // Interpolation distance always positive

  a_stencil.clear();

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
    a_stencil.clear();
    a_stencil.add(a_vof, 1-x);

    for (int i = 0; i < numOtherVoFs; i++){
      a_stencil.add(otherVoFs[i], x/numOtherVoFs);
    }
  }
  else{ // Can't do deriv in this direction.
#if 1
    std::cout << "could not find linear stencil for iv = " << a_vof.gridIndex() << std::endl;
#endif
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);
  }

  found_stencil = true;

#if 1 // Debug
  Real sumweights = 0.0;
  for (int i = 0; i < a_stencil.size(); i++){
    const Real w = a_stencil.weight(i);
    if(w < 0.0) MayDay::Abort("linear negative weight");
    sumweights += w;
  }
  if(Abs((sumweights - 1.0)) > 1.E-5) {
    std::cout << sumweights << std::endl;
    MayDay::Abort("linear sumweights is fucked up");
  }
#endif

  return found_stencil;
}

// 
bool stencil_ops::compute_interp_stencil_2D(VoFStencil&          a_stencil,
					    const RealVect&      a_centroid,
					    const VolIndex&      a_vof,
					    const ProblemDomain& a_domain,
					    const EBISBox&       a_ebisbox,
					    const int            a_nointerp_dir){
  CH_TIME("compute_interp_stencil2D");
  CH_assert(SpaceDim == 2 || SpaceDim == 3);

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

  // If we can only interpolate in one of the direction, return the 1D stnecil
  if(!canInterpDir0 && !canInterpDir1){
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);
  }
  else if( canInterpDir0 && !canInterpDir1){
    found_stencil = stencil_ops::compute_interp_stencil_1D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox, dir0);
  }
  else if(!canInterpDir0 &&  canInterpDir1){
    found_stencil = stencil_ops::compute_interp_stencil_1D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox, dir1);
  }
  else{ // OK, we know we can interpolate in both directions from a_vof. Let's check if we can find a bilinear stencil now

    // 1. First, starting on a_vof, compute the linear interpolation stencil in 1D along dir0 and dir1
    VoFStencil curStenDir0;
    VoFStencil curStenDir1;
    const bool foundFirstStencil0 = stencil_ops::compute_interp_stencil_1D(curStenDir0,a_centroid,a_vof,a_domain,a_ebisbox,dir0);
    const bool foundFirstStencil1 = stencil_ops::compute_interp_stencil_1D(curStenDir1,a_centroid,a_vof,a_domain,a_ebisbox,dir1);
    if(!foundFirstStencil0 || !foundFirstStencil1){ // If we made it into this loop, this shouldn't happen
      MayDay::Abort("stencil_ops - logic bust");
    }

    // 2. For each of the VoFs on appropriate side in dir0 and dir1 (usually just the one VoF), compute the average stencil
    VoFStencil otherVoFsInterpStencilDir0;
    VoFStencil otherVoFsInterpStencilDir1;
    
    int numStencilsAdded0 = 0;
    int numStencilsAdded1 = 0;
    
    bool foundSecondStencil0 = false;
    bool foundSecondStencil1 = false;
    
    for (int i = 0; i<otherVoFsInDir1.size(); i++){
      VoFStencil sten;
      const bool found = stencil_ops::compute_interp_stencil_1D(sten,a_centroid,otherVoFsInDir1[i],a_domain,a_ebisbox,dir0);
      
      if(found){
	foundSecondStencil0 = true;
	otherVoFsInterpStencilDir0 += sten;
	numStencilsAdded0++;
      }
    }
    for (int i = 0; i<otherVoFsInDir0.size(); i++){
      VoFStencil sten;
      const bool found = stencil_ops::compute_interp_stencil_1D(sten,a_centroid,otherVoFsInDir0[i],a_domain,a_ebisbox,dir1);
      
      if(found){
	foundSecondStencil1 = true;
	otherVoFsInterpStencilDir1 += sten;
	numStencilsAdded1++;
      }
    }
    
    // 3. Check if we could get the interpolation stencil from the other VoFs. If we could, scale them with the number of VoFs
    //    we accesses
    if(numStencilsAdded0 > 0) otherVoFsInterpStencilDir0 *= 1./numStencilsAdded0;
    if(numStencilsAdded1 > 0) otherVoFsInterpStencilDir1 *= 1./numStencilsAdded1;


  // 4. We now have interpolation stencils that should not have negative weights. Do some debugging
#if 1
  if(foundFirstStencil){
    Real sumweights = 0.0;
    for (int i = 0; i < curVoFInterpStencilDir0.size(); i++){
      const Real w = curVoFInterpStencilDir0.weight(i);
      if(w < 0.0) {
	MayDay::Warning("stencils_ops::bilinear - negative weight in sten0");
      }
      sumweights += w;
    }
    if(Abs(sumweights - 1.0) > 1.E-5){
      MayDay::Warning("stencil_ops::bilinear sten0weights do not sum to 1");
    }
  }
  
  if(foundSecondStencil){
    Real sumweights = 0.0;
    for (int i = 0; i < otherVoFsInterpStencilDir0.size(); i++){
      const Real w = otherVoFsInterpStencilDir0.weight(i);
      if(w < 0.0) {
	MayDay::Warning("stencils_ops::bilinear - negative weight in sten1");
      }
      sumweights += w;
    }
    if(Abs(sumweights - 1.0) > 1.E-5){
      MayDay::Warning("stencil_ops::bilinear stencil sten1 weights do not sum to 1");
    }
  }
#endif

  a_stencil.clear();
  if(foundFirstStencil){
    a_stencil += curVoFInterpStencilDir0; // OK, we can interpolate in dir0 from the current stencil

    if(foundSecondStencil){ // Can also interpolate in dir 1 because we found an interpolation stencil there
      const Real d = a_centroid[dir]*HiLoDir1;

      a_stencil.
    }
    else{
      a_stencil += curVoFInterStencilDir0; // Only 1D interpolation along dir0
    }
  }
  else{ // Sorry, no interpolation today
    a_stencil.add(a_vof, 1.0);
  }


  
  
#if 0
    // To compute the stencil, first do 1D linear interpolation along dir 0 for two rows of VoFs (one displaced along dir1)
    const IntVect iv0   = a_vof.gridIndex();
    const IntVect iv1   = iv0 + BASISV(dir1)*HiLoDir1;
    const VolIndex vof0 = VolIndex(iv0,0);
    const VolIndex vof1 = VolIndex(iv1,0);

    // Compute two linear stencils
    VoFStencil sten0;
    VoFStencil sten1;
    const bool foundFirstStencil  = stencil_ops::compute_interp_stencil_1D(sten0, a_centroid, vof0, a_domain, a_ebisbox, dir0);
    const bool foundSecondStencil = stencil_ops::compute_interp_stencil_1D(sten1, a_centroid, vof1, a_domain, a_ebisbox, dir0);

#if 1 // Debug
    for (int i = 0; i < sten0.size(); i++){
      const Real w = sten0.weight(i);
      if(w < 0.0) MayDay::Abort("stencils_ops STOP - negative weight in sten0");
    }
    for (int i = 0; i < sten1.size(); i++){
      const Real w = sten1.weight(i);
      if(w < 0.0) MayDay::Abort("stencils_ops STOP - negative weight in sten1");
    }
#endif

    // If both stencils were found, do linear interpolation again, this time along dir1 
    if(foundFirstStencil && foundSecondStencil){
      found_stencil = true;
      a_stencil.clear();

      const Real d = a_centroid[dir1]*HiLoDir1;

      const Real w0 = (1-d);
      const Real w1 = d;
      
      sten0 *= w0;
      sten1 *= w1;

      a_stencil += sten0;
      a_stencil += sten1;
    }
    else{
      found_stencil = false;
      a_stencil.clear();
      a_stencil.add(a_vof, 1.0);
    }
  }
  else if(!really2D){ // If this triggers, the centroid lies on a straight line connecting two VoFs
    // Get the interpolation direction and interpolate
    int interpDir;
    if(a_centroid[dir0] <= stencil_ops::tolerance){
      interpDir = dir1;
    }
    else if(a_centroid[dir1] <= stencil_ops::tolerance){
      interpDir = dir0;
    }
    found_stencil = stencil_ops::compute_interp_stencil_1D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox, interpDir);
  }
    
  //
  if(!found_stencil){
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);
  }

#if 1 // Debug
  Real sumweights = 0.0;
  for (int i = 0; i < a_stencil.size(); i++){
    const Real w = a_stencil.weight(i);
    if(w < 0.0) MayDay::Abort("bilinear negative weight");
    sumweights += w;
  }
  if(Abs((sumweights - 1.0)) > 1.E-5) {
    std::cout << sumweights << std::endl;
    MayDay::Warning("bilinear sumweights is fucked up");
  }

#endif
#endif

  return found_stencil;
}

bool stencil_ops::compute_interp_stencil_3D(VoFStencil&          a_stencil,
					    const RealVect&      a_centroid,
					    const VolIndex&      a_vof,
					    const ProblemDomain& a_domain,
					    const EBISBox&       a_ebisbox){
  CH_TIME("compute_interp_stencil_3D");
  CH_assert(SpaceDim == 3);


  a_stencil.clear();
  
  bool found_stencil = false;

  // Check if centroid lies in one of the coordinate planes connecting cell centers
  const Real tol = stencil_ops::tolerance;
  bool really3D = true;
  if(Abs(a_centroid[0]) <= tol || Abs(a_centroid[1]) <= tol || Abs(a_centroid[2]) <= tol){
    really3D = false;
  }

  if(really3D){

    const int HiLoDir2 = a_centroid[2] > 0. ? 1 : -1;

    // Compute bilinear stencils for two rows of VoFs 
    const IntVect iv0 = a_vof.gridIndex();
    const IntVect iv1 = iv0 + BASISV(2)*HiLoDir2;
    const VolIndex vof0 = VolIndex(iv0,0);
    const VolIndex vof1 = VolIndex(iv1,0);

    VoFStencil sten0;
    VoFStencil sten1;
    const bool foundFirstStencil  = stencil_ops::compute_interp_stencil_2D(sten0, a_centroid, vof0, a_domain, a_ebisbox, 2);
    const bool foundSecondStencil = stencil_ops::compute_interp_stencil_2D(sten1, a_centroid, vof1, a_domain, a_ebisbox, 2);

    if(foundFirstStencil && foundSecondStencil){
      found_stencil = true;
      a_stencil.clear();

      // Interpolation distances
      const Real z = a_centroid[2]*HiLoDir2;
      // Linearly interpolate the stencil
      const Real w0 = (1 - z);
      const Real w1 = z;

      sten0 *= w0;
      sten1 *= w1;

      a_stencil += sten0;
      a_stencil += sten1;
    }
  }
  else if(!really3D){ // Either bilinear or linear interpolation
    int noInterpDir = 0;
    for (int dir = 0; dir < SpaceDim; dir++){
      noInterpDir = Abs(a_centroid[dir]) <= tol ? dir : noInterpDir;
    }
    found_stencil = compute_interp_stencil_2D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox, noInterpDir);
  }

  if(!found_stencil){
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);
  }

#if 1 // Debug
  Real sumweights = 0.0;
  for (int i = 0; i < a_stencil.size(); i++){
    const Real w = a_stencil.weight(i);
    if(w < 0.0) MayDay::Abort("trilinear stencil weight < 0.0");
    sumweights += w;
  }
  if(Abs((sumweights - 1.0)) > 1.E-5) MayDay::Abort("trilinear sumweights is fucked up");
#endif

  return true;
}

