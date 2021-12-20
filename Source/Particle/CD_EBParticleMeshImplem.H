/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBParticleMeshImplem.H
  @brief  Implementatoin of CD_EBParticleMesh.H
  @author Robert Marskar
*/

#ifndef CD_EBParticleMeshImplem_H
#define CD_EBParticleMeshImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBParticleMesh.H>
#include <CD_BoxLoops.H>
#include <CD_EBParticleMeshF_F.H>
#include <CD_NamespaceHeader.H>

template <class P, const Real&(P::*particleScalarField)() const>
void EBParticleMesh::deposit(const List<P>& a_particleList, EBCellFAB& a_rho, const DepositionType a_depositionType, const bool a_forceIrregNGP) const {
  CH_TIME("EBParticleMesh::deposit");

  // TLDR: This is a jack-of-all-trades deposition function. The user will use this function to supply a pointer to the field that will be
  //       deposited. As per API, this function must be of the type 'const Real& myParticleClass::myDepositionField() const'

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();    
    
    this->depositParticle(a_rho, 
			  m_probLo, 
			  m_dx, 
			  curPosition, 
			  curStrength, 
			  a_depositionType,
			  a_forceIrregNGP);
  }  
}

template <class P, Real(P::*particleScalarField)() const>
void EBParticleMesh::deposit(const List<P>& a_particleList, EBCellFAB& a_rho, const DepositionType a_depositionType, const bool a_forceIrregNGP) const {
  CH_TIME("EBParticleMesh::deposit");

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();    
    
    this->depositParticle(a_rho, 
			  m_probLo, 
			  m_dx, 
			  curPosition, 
			  curStrength, 
			  a_depositionType,
			  a_forceIrregNGP);
  }  
}

template <class P, const Real&(P::*particleScalarField)() const>
void EBParticleMesh::deposit2(const List<P>& a_particleList, EBCellFAB& a_rho, const DepositionType a_depositionType, const bool a_forceIrregNGP) const {
  CH_TIME("EBParticleMesh::deposit2");

  // TLDR: This is a jack-of-all-trades deposition function. The user will use this function to supply a pointer to the field that will be
  //       deposited. As per API, this function must be of the type 'const Real& myParticleClass::myDepositionField() const'  

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();    

    this->depositParticle2(a_rho, 
			   m_probLo, 
			   m_dx, 
			   curPosition, 
			   curStrength, 
			   a_depositionType,
			   a_forceIrregNGP);
  }  
}

template <class P,  Real(P::*particleScalarField)() const>
void EBParticleMesh::deposit2(const List<P>& a_particleList, EBCellFAB& a_rho, const DepositionType a_depositionType, const bool a_forceIrregNGP) const {
  CH_TIME("EBParticleMesh::deposit2");

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();    

    this->depositParticle2(a_rho, 
			   m_probLo, 
			   m_dx, 
			   curPosition, 
			   curStrength, 
			   a_depositionType,
			   a_forceIrregNGP);
  }  
}

template <class P, const Real&(P::*particleScalarField)() const>
void EBParticleMesh::deposit4(const List<P>& a_particleList, EBCellFAB& a_rho, const DepositionType a_depositionType, const bool a_forceIrregNGP) const {
  CH_TIME("EBParticleMesh::deposit4");

  // TLDR: This is a jack-of-all-trades deposition function. The user will use this function to supply a pointer to the field that will be
  //       deposited. As per API, this function must be of the type 'const Real& myParticleClass::myDepositionField() const'    

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();    
    
    this->depositParticle4(a_rho, 
			   m_probLo, 
			   m_dx, 
			   curPosition, 
			   curStrength, 
			   a_depositionType,
			   a_forceIrregNGP);
  }  
}

template <class P,  Real(P::*particleScalarField)() const>
void EBParticleMesh::deposit4(const List<P>& a_particleList, EBCellFAB& a_rho, const DepositionType a_depositionType, const bool a_forceIrregNGP) const {
  CH_TIME("EBParticleMesh::deposit4");

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    const P&        curParticle = lit();
    const RealVect& curPosition = curParticle.position();
    const Real&     curStrength = (curParticle.*particleScalarField)();    
    
    this->depositParticle4(a_rho, 
			   m_probLo, 
			   m_dx, 
			   curPosition, 
			   curStrength, 
			   a_depositionType,
			   a_forceIrregNGP);
  }  
}

template <class P, Real&(P::*particleScalarField)()>
void EBParticleMesh::interpolate(List<P>&             a_particleList,
				 const EBCellFAB&     a_meshScalarField,
				 const DepositionType a_interpType,
				 const bool           a_forceIrregNGP) const {
  CH_TIME("EBParticleMesh::interpolate(Real)");
  
  CH_assert(a_meshScalarField.nComp() == 1);

  // TLDR: This is a jack-of-all-trades interpolation function. The user will use this function to supply a pointer to the field that will be
  //       interpolated to. As per API, this function must be of the type 'Real& myParticleClass::myScalarVariable()'      

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    P&              curParticle  = lit();
    const RealVect& curPosition  = curParticle.position();

    // Note: We MUST initialize = 0 because the higher-order kernels do not reset the field. It doesn't matter if we do it here
    //       or in the Fortran kernel.     
    Real& curParticleField = (curParticle.*particleScalarField)();
    curParticleField = 0.0;
    this->interpolateParticle(curParticleField,
			      a_meshScalarField,
			      m_probLo,
			      m_dx,
			      curPosition,
			      a_interpType,
			      a_forceIrregNGP);
  }
}

template <class P, RealVect&(P::*particleVectorField)()>
void EBParticleMesh::interpolate(List<P>&             a_particleList,
				 const EBCellFAB&     a_meshVectorField,
				 const DepositionType a_interpType,
				 const bool           a_forceIrregNGP) const {
  CH_TIME("EBParticleMesh::interpolate(RealVect)");
  
  CH_assert(a_meshVectorField.nComp() == SpaceDim);

  // TLDR: This is a jack-of-all-trades interpolation function. The user will use this function to supply a pointer to the field that will be
  //       interpolated to. As per API, this function must be of the type 'RealVect& myParticleClass::myVectorVariable()'        

  for (ListIterator<P> lit(a_particleList); lit; ++lit){
    P&              curParticle  = lit();
    const RealVect& curPosition  = curParticle.position();

    // Note: We MUST initialize = 0 because the higher-order kernels do not reset the field. It doesn't matter if we do it here
    //       or in the Fortran kernel.     
    RealVect& curParticleField = (curParticle.*particleVectorField)();
    curParticleField = RealVect::Zero;
    this->interpolateParticle(curParticleField,
			      a_meshVectorField,
			      m_probLo,
			      m_dx,
			      curPosition,
			      a_interpType,
			      a_forceIrregNGP);
  }
}

inline
void EBParticleMesh::depositParticle(EBCellFAB&           a_rho,
				     const RealVect&      a_probLo,
				     const RealVect&      a_dx,
				     const RealVect&      a_position,
				     const Real&          a_mass,
				     const DepositionType a_depositionType,
				     const bool           a_forceIrregNGP) const {
  CH_TIME("EBParticleMesh::depositParticle");
  
  CH_assert(a_rho.nComp() == 1);

  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.  

  constexpr int comp  = 0;
  constexpr int nComp = 1;  

  // Nifty lambda for converting RealVect position to IntVect (lower-left corner). Note that the input vector
  // is the displacement from a_probLo and not the physical position. 
  auto getCellIndex = [&dx=a_dx] (const RealVect a_rv) -> IntVect {
    return IntVect(D_DECL(std::floor(a_rv[0]/dx[0]), std::floor(a_rv[1]/dx[1]), std::floor(a_rv[2]/dx[2])));
  };

  // Get grid cell corresponding to this particle. 
  const IntVect particleIndex = getCellIndex(a_position - a_probLo);

  // Assertion -- particle must live on this patch. 
  CH_assert(m_region.contains(particleIndex));

  // Factors needed for the kernels. density is the, well, density = mass/volume. The little lambda will give us the distance
  // between the cell we're looking at and the particle position (in units of the grid resolution)
  const Real density = a_mass/std::pow(a_dx[0], SpaceDim);

  auto particleDisplacement = [&] (const IntVect& iv) -> RealVect {
    return (a_probLo + (RealVect(iv) + 0.5*RealVect::Unit)*a_dx - a_position)/a_dx;
  };
  
  // Get regular data. 
  FArrayBox& rho = a_rho.getFArrayBox();  
  
  // We can force NGP deposition in cut-cells if we want. 
  if(m_ebisbox.isIrregular(particleIndex) && a_forceIrregNGP){
    rho(particleIndex, comp) += density;    
  }
  else{
    switch (a_depositionType){
    case DepositionType::NGP:
      {
	rho(particleIndex, comp) += density;
	break;
      }
    case DepositionType::CIC:
      {
	// Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
	// the box corresponding to the cloud. 
	const IntVect loIndex  = getCellIndex(a_position - a_probLo - 0.5*a_dx);
	const Box     cicBox = Box(loIndex, loIndex + IntVect::Unit);

	CH_assert(a_rho.box().contains(cicBox));		

	// This is the cloud-in-cell deposition kernel.
	auto cicKernel = [&](const IntVect& iv) -> void {
	  Real weight = density;

	  // L is the distance between the grid cell and the particle. 
	  const RealVect L = particleDisplacement(iv);
	  for (int dir = 0; dir < SpaceDim; dir++){
	    weight *= (1. - std::abs(L[dir]));
	  }

	  rho(iv,comp) += weight;
	};

	// Add mass to cells. 
	BoxLoops::loop(cicBox, cicKernel);

	break;
      }
    case DepositionType::TSC:
      {

	// Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
	// the box corresponding to the cloud. 
	const IntVect loIndex = getCellIndex(a_position - a_probLo - a_dx);
	const Box     tscBox  = Box(loIndex, loIndex + 2*IntVect::Unit);

	CH_assert(a_rho.box().contains(tscBox));

	// This is the triangle-shaped cloud kernel. 
	auto tscKernel = [&](const IntVect& iv) -> void {
	  Real weight = density;

	  const RealVect L = particleDisplacement(iv);	  
	  for (int dir = 0; dir < SpaceDim; dir++){
	    const Real l = std::abs(L[dir]);

	    if(l < 0.5){
	      weight *= 0.75 - l*l;
	    }
	    else{
	      weight *= 0.5 * (1.5 - l) * (1.5 - l);
	    }
	  }

	  rho(iv, comp) += weight;
	};

	// Run the kernel. 
	BoxLoops::loop(tscBox, tscKernel);

	break;
      }
    case DepositionType::W4:
      {

	// Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
	// the box corresponding to the cloud. 
	const IntVect loIndex = getCellIndex(a_position - a_probLo - 1.5*a_dx);
	const Box     w4Box   = Box(loIndex, loIndex + 3*IntVect::Unit);

	CH_assert(a_rho.box().contains(w4Box));

	// This is the fourth order kernel.
	auto w4Kernel = [&](const IntVect& iv) -> void {
	  Real weight = density;
	  
	  // L is the distance between the grid cell and the particle. 
	  const RealVect L = particleDisplacement(iv);	  	  
	  for (int dir = 0; dir < SpaceDim; dir++){
	    const Real l = std::abs(L[dir]);

	    if(l < 1.0) {
	      weight *= 1.0 - 2.5*l*l + 1.5*l*l*l;
	    }
	    else {
	      weight *= 0.5 * (2. - l) * (2. - l) * (1. - l);
	    }
	  }

	  rho(iv, comp) += weight;
	};

	// Run the kernel.
	BoxLoops::loop(w4Box, w4Kernel);

	break;
      }
    default:
      {
	MayDay::Error("EBParticleMesh::depositParticle - logic bust, unknown particle deposition.");
	break;
      }
    }
  }
}

inline
void EBParticleMesh::depositParticle2(EBCellFAB&           a_rho,
				      const RealVect&      a_probLo,
				      const RealVect&      a_dx,
				      const RealVect&      a_position,
				      const Real&          a_mass,
				      const DepositionType a_depositionType,
				      const bool           a_forceIrregNGP) const {
  CH_TIME("EBParticleMesh::depositParticle2");
  
  CH_assert(a_rho.nComp() == 1);

  constexpr int comp  = 0;
  constexpr int nComp = 1;    
  
  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.
  //
  //       Note that this is the version which deposits with particle widths that are 4 times the "usual" width. Currently, only NGP
  //       and CIC is supported.

  // Nifty lambda for converting RealVect position to IntVect (lower-left corner). Note that the input vector
  // is the displacement from a_probLo and not the physical position. 
  auto getCellIndex = [&dx=a_dx] (const RealVect a_rv) -> IntVect {
    return IntVect(D_DECL(std::floor(a_rv[0]/dx[0]), std::floor(a_rv[1]/dx[1]), std::floor(a_rv[2]/dx[2])));
  };

  // Get grid cell corresponding to this particle. 
  const IntVect particleIndex = getCellIndex(a_position - a_probLo);

  // Assertion -- particle must live on this patch. 
  CH_assert(m_region.contains(particleIndex));

  // Factors needed for the kernels. density is the, well, density = mass/volume. The little lambda will give us the distance
  // between the cell we're looking at and the particle position (in units of the grid resolution)
  const Real density = a_mass/std::pow(a_dx[0], SpaceDim);

  auto particleDisplacement = [&] (const IntVect& iv) -> RealVect {
    return (a_probLo + (RealVect(iv) + 0.5*RealVect::Unit)*a_dx - a_position)/a_dx;
  };
  
  // Get regular data. 
  FArrayBox& rho = a_rho.getFArrayBox();  

  // Force NGP in cut-cells if we want. 
  if(m_ebisbox.isIrregular(particleIndex) && a_forceIrregNGP){
    rho(particleIndex, comp) += density;        
  }
  else{
    switch (a_depositionType){
    case DepositionType::NGP:
      {
	rho(particleIndex, comp) += density;        	
	break;
      }
    case DepositionType::CIC:
      {
	// Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
	// the box corresponding to the cloud. 
	const IntVect loIndex = getCellIndex(a_position - a_probLo - 1.0*a_dx);
	const Box     cic2Box = Box(loIndex, loIndex + 2*IntVect::Unit);

	CH_assert(a_rho.box().contains(cic2Box));	

	// This is the kernel for a CIC particle with with 2*dx
	auto cic2Kernel = [&](const IntVect& iv) -> void {
	  Real weight = density;
	  
	  // L is the distance between the grid cell and the particle. 
	  const RealVect L = particleDisplacement(iv);
	  for (int dir = 0; dir < SpaceDim; dir++){
	    const Real l = std::abs(L[dir]);

	    if(l > 0.5) {
	      weight *= 0.5 * (1.5 - l);
	    }
	    else{
	      weight *= 0.5;
	    }
	  }
	  
	  rho(iv, comp) += weight;	    	  
	};

	BoxLoops::loop(cic2Box, cic2Kernel);
	
	break;
      }
    default:
      {
	MayDay::Error("EBParticleMesh::depositParticle2 - Invalid deposition type - only NGP and CIC supported for this deposition method. TSC/W4 have not been worked out.");
	break;
      }
    }
  }
}

inline
void EBParticleMesh::depositParticle4(EBCellFAB&           a_rho,
				      const RealVect&      a_probLo,
				      const RealVect&      a_dx,
				      const RealVect&      a_position,
				      const Real&          a_mass,
				      const DepositionType a_depositionType,
				      const bool           a_forceIrregNGP) const {
  CH_TIME("EBParticleMesh::depositParticle2");
  
  CH_assert(a_rho.nComp() == 1);

  constexpr int comp  = 0;
  constexpr int nComp = 1;    
  
  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.
  //
  //       Note that this is the version which deposits with particle widths that are 4 times the "usual" width. Currently, only NGP
  //       and CIC is supported.

  // Nifty lambda for converting RealVect position to IntVect (lower-left corner). Note that the input vector
  // is the displacement from a_probLo and not the physical position. 
  auto getCellIndex = [&dx=a_dx] (const RealVect a_rv) -> IntVect {
    return IntVect(D_DECL(std::floor(a_rv[0]/dx[0]), std::floor(a_rv[1]/dx[1]), std::floor(a_rv[2]/dx[2])));
  };

  // Get grid cell corresponding to this particle. 
  const IntVect particleIndex = getCellIndex(a_position - a_probLo);

  // Assertion -- particle must live on this patch. 
  CH_assert(m_region.contains(particleIndex));

  // Factors needed for the kernels. density is the, well, density = mass/volume. The little lambda will give us the distance
  // between the cell we're looking at and the particle position (in units of the grid resolution)
  const Real density = a_mass/std::pow(a_dx[0], SpaceDim);

  auto particleDisplacement = [&] (const IntVect& iv) -> RealVect {
    return (a_probLo + (RealVect(iv) + 0.5*RealVect::Unit)*a_dx - a_position)/a_dx;
  };
  
  // Get regular data. 
  FArrayBox& rho = a_rho.getFArrayBox();  

  // Force NGP in cut-cells if we want. 
  if(m_ebisbox.isIrregular(particleIndex) && a_forceIrregNGP){
    rho(particleIndex, comp) += density;        
  }
  else{
    switch (a_depositionType){
    case DepositionType::NGP:
      {
	rho(particleIndex, comp) += density;        	
	break;
      }
    case DepositionType::CIC:
      {
	// Compute the index of the cell that contains the lower-left corner of this particle cloud. Also compute
	// the box corresponding to the cloud. 
	const IntVect loIndex = getCellIndex(a_position - a_probLo - 2.0*a_dx);
	const Box     cic4Box = Box(loIndex, loIndex + 4*IntVect::Unit);

	// This is a warnign that will trigger if there are not enough ghost cells in the box. In principle, there needs to be at least four ghost cells
	// available in order to use this kernel. That being said -- this routine is called on refined coarse-grid data so that should not be a problem (normally). 
	CH_assert(a_rho.box().contains(cic4Box));

	// This is the kernel for a CIC particle with with 4*dx
	auto cic4Kernel = [&](const IntVect& iv) -> void {
	  Real weight = density;
	  
	  // L is the distance between the grid cell and the particle -- in units of dx. 
	  const RealVect L = particleDisplacement(iv);
	  for (int dir = 0; dir < SpaceDim; dir++){
	    const Real l = std::abs(L[dir]);

	    if(l >= 1.5) {
	      weight *= 0.25 * (2.5 - l);
	    }
	    else{
	      weight *= 0.25;
	    }
	  }
	  
	  rho(iv, comp) += weight;	    	  
	};

	BoxLoops::loop(cic4Box, cic4Kernel);

	break;
      }
    default:
      {
	MayDay::Error("EBParticleMesh::depositParticle4 - Invalid deposition type - only NGP and CIC supported for this deposition method. TSC/W4 have not been worked out.");
	break;
      }
    }
  }
}

inline
void EBParticleMesh::interpolateParticle(Real&                a_particleField,
					 const EBCellFAB&     a_meshField,
					 const RealVect&      a_probLo,
					 const RealVect&      a_dx,
					 const RealVect&      a_position,
					 const DepositionType a_interpType,
					 const bool           a_forceIrregNGP) const {
  CH_assert(a_meshField.nComp() == 1);

  // TLDR: This routine is just like usual particle-mesh interpolation. 

  constexpr int comp = 0;
  
  const RealVect rv  = (a_position - a_probLo)/a_dx;
  const IntVect  iv  = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  CH_assert(m_region.contains(iv));  

  const FArrayBox& meshField = a_meshField.getFArrayBox();      

  // Irregular cells always do an NGP deposit to prevent clouds leaking into the other side. 
  if(m_ebisbox.isIrregular(iv) && a_forceIrregNGP){
    FORT_NGP_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				CHF_CONST_FRA1(meshField,comp),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position));
  }
  else if(m_ebisbox.isCovered(iv)){ // Need to set to something. 
    a_particleField = 0.0;
  }  
  else{
    switch (a_interpType) {
    case DepositionType::NGP:
      FORT_NGP_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				  CHF_CONST_FRA1(meshField,comp),
				  CHF_CONST_REALVECT(a_probLo),
				  CHF_CONST_REALVECT(a_dx),
				  CHF_CONST_REALVECT(a_position));

      break;
    case DepositionType::CIC:
      FORT_CIC_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				  CHF_CONST_FRA1(meshField,comp),
				  CHF_CONST_REALVECT(a_probLo),
				  CHF_CONST_REALVECT(a_dx),
				  CHF_CONST_REALVECT(a_position));
      break;
    case DepositionType::TSC:
      FORT_TSC_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				  CHF_CONST_FRA1(meshField,comp),
				  CHF_CONST_REALVECT(a_probLo),
				  CHF_CONST_REALVECT(a_dx),
				  CHF_CONST_REALVECT(a_position));
      break;
    case DepositionType::W4:
      FORT_W4_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				 CHF_CONST_FRA1(meshField,comp),
				 CHF_CONST_REALVECT(a_probLo),
				 CHF_CONST_REALVECT(a_dx),
				 CHF_CONST_REALVECT(a_position));
      break;
    default:
      MayDay::Error("EBParticleMesh::interpolateParticle(Real) - Invalid interpolation type requested.");
    }
  }
}

inline
void EBParticleMesh::interpolateParticle(RealVect&            a_particleField,
					 const EBCellFAB&     a_meshField,
					 const RealVect&      a_probLo,
					 const RealVect&      a_dx,
					 const RealVect&      a_position,
					 const DepositionType a_interpType,
					 const bool           a_forceIrregNGP) const {
  CH_assert(a_meshField.nComp() == SpaceDim);

  // TLDR: This routine is just like usual particle-mesh interpolation. 

  constexpr int comp = 0;
  
  const RealVect rv  = (a_position - a_probLo)/a_dx;
  const IntVect  iv  = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  CH_assert(m_region.contains(iv));

  const FArrayBox& meshField = a_meshField.getFArrayBox();

  // Irregular cells always do an NGP deposit to prevent clouds leaking into the other side.
  if(m_ebisbox.isIrregular(iv) && a_forceIrregNGP){    
    FORT_NGP_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				CHF_CONST_FRA(meshField),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position));
  }
  else if(m_ebisbox.isCovered(iv)){ // Need to set to something. 
    a_particleField = RealVect::Zero;
  }
  else{
    switch (a_interpType) {
    case DepositionType::NGP:
      FORT_NGP_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				  CHF_CONST_FRA(meshField),
				  CHF_CONST_REALVECT(a_probLo),
				  CHF_CONST_REALVECT(a_dx),
				  CHF_CONST_REALVECT(a_position));

      break;
    case DepositionType::CIC:
      FORT_CIC_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				  CHF_CONST_FRA(meshField),
				  CHF_CONST_REALVECT(a_probLo),
				  CHF_CONST_REALVECT(a_dx),
				  CHF_CONST_REALVECT(a_position));
      break;
    case DepositionType::TSC:
      FORT_TSC_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				  CHF_CONST_FRA(meshField),
				  CHF_CONST_REALVECT(a_probLo),
				  CHF_CONST_REALVECT(a_dx),
				  CHF_CONST_REALVECT(a_position));
      break;
    case DepositionType::W4:
      FORT_W4_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				 CHF_CONST_FRA(meshField),
				 CHF_CONST_REALVECT(a_probLo),
				 CHF_CONST_REALVECT(a_dx),
				 CHF_CONST_REALVECT(a_position));
      break;
    default:
      MayDay::Error("EBParticleMesh::interpolateParticle(RealVect) - invalid interpolation type requested.");
    }
  }
}



#include <CD_NamespaceFooter.H>

#endif
