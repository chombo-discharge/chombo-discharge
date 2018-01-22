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
    do_interp = false;
  }

  if(do_interp){
#if CH_SPACEDIM == 2 
    found_stencil = compute_interp_stencil_2D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox);
#elif CH_SPACEDIM == 3
    found_stencil = compute_interp_stencil_3D(a_stencil, a_centroid, a_vof, a_domain, a_ebisbox);
#endif
  }
  else{ 
    found_stencil = true;
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);
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
  const IntVect iv0    = a_vof.gridIndex();
  const IntVect iv1    = iv0 + BASISV(a_interp_dir)*HiLoInterp;
  const Real x         = a_centroid[a_interp_dir]*HiLoInterp;   // Interpolation distance always positive

  if(a_domain.contains(iv0) && a_domain.contains(iv1)){
    if(!a_ebisbox.isCovered(iv0) && !a_ebisbox.isCovered(iv1)){
      found_stencil  = true;
      const VolIndex vof0 = VolIndex(iv0,0);
      const VolIndex vof1 = VolIndex(iv1,0);
      
      const Real w0 = (1-x);
      const Real w1 = x;
      
      a_stencil.clear();
      a_stencil.add(vof0, w0);
      a_stencil.add(vof1, w1);
    }
  }

  if(!found_stencil){
    a_stencil.clear();
  }

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


  bool really2D  = true;
  if(Abs(a_centroid[dir0]) <= stencil_ops::tolerance || Abs(a_centroid[dir1]) <= stencil_ops::tolerance){ // Must do 2D
    really2D = false; 
  }
  else if(Abs(a_centroid[dir0]) <= stencil_ops::tolerance && Abs(a_centroid[dir1]) <= stencil_ops::tolerance){ // Cell center
    found_stencil = true;
    a_stencil.clear();
    a_stencil.add(a_vof, 1.0);
    return found_stencil;
  }
  
  if(really2D){
    const int HiLoDir0 = a_centroid[dir0] > 0. ? 1 : -1;
    const int HiLoDir1 = a_centroid[dir1] > 0. ? 1 : -1;

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
  }

  return found_stencil;
}

bool stencil_ops::compute_interp_stencil_3D(VoFStencil&          a_stencil,
					    const RealVect&      a_centroid,
					    const VolIndex&      a_vof,
					    const ProblemDomain& a_domain,
					    const EBISBox&       a_ebisbox){
  CH_TIME("compute_interp_stencil_3D");
  CH_assert(SpaceDim == 3);
  
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
}
