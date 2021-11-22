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

template <class P, Real&(P::*particleScalarField)()>
void EBParticleMesh::interpolate(List<P>&             a_particleList,
				 const EBCellFAB&     a_meshScalarField,
				 const DepositionType a_interpType,
				 const bool           a_forceIrregNGP) const {

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
  CH_assert(a_rho.nComp() == 0);
  
  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.

  constexpr int comp = 0;

  // Check if the cell is irregular. 
  const RealVect rv = (a_position - a_probLo)/a_dx;
  const IntVect  iv = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  FArrayBox& rho = a_rho.getFArrayBox();

  // We can force NGP deposition in cut-cells if we want. 
  if(m_ebisbox.isIrregular(iv) && a_forceIrregNGP){
    FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(rho, comp),
			    CHF_CONST_REALVECT(a_probLo),
			    CHF_CONST_REALVECT(a_dx),
			    CHF_CONST_REALVECT(a_position),
			    CHF_CONST_REAL(a_mass));
  }
  else{
    switch (a_depositionType){
    case DepositionType::NGP:
      {
	FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(rho, comp),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position),
				CHF_CONST_REAL(a_mass));
	break;
      }
    case DepositionType::CIC:
      {
	FORT_CIC_DEPOSIT_SCALAR(CHF_FRA1(rho, comp),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position),
				CHF_CONST_REAL(a_mass));
	break;
      }
    case DepositionType::TSC:
      {
	FORT_TSC_DEPOSIT_SCALAR(CHF_FRA1(rho, comp),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position),
				CHF_CONST_REAL(a_mass));
	break;
      }
    case DepositionType::W4:
      {
	FORT_W4_DEPOSIT_SCALAR(CHF_FRA1(rho, comp),
			       CHF_CONST_REALVECT(a_probLo),
			       CHF_CONST_REALVECT(a_dx),
			       CHF_CONST_REALVECT(a_position),
			       CHF_CONST_REAL(a_mass));
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
  CH_assert(a_rho.nComp() == 0);
  
  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.
  //
  //       Note that this is the version which deposits with particle widths that are 4 times the "usual" width. Currently, only NGP
  //       and CIC is supported.

  constexpr int comp = 0;
  
  const RealVect rv  = (a_position - a_probLo)/a_dx;
  const IntVect  iv  = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  FArrayBox& rho = a_rho.getFArrayBox();  

  // Force NGP in cut-cells if we want. 
  if(m_ebisbox.isIrregular(iv) && a_forceIrregNGP){
    FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(rho, comp),
			    CHF_CONST_REALVECT(a_probLo),
			    CHF_CONST_REALVECT(a_dx),
			    CHF_CONST_REALVECT(a_position),
			    CHF_CONST_REAL(a_mass));
  }
  else{
    switch (a_depositionType){
    case DepositionType::NGP:
      {
	FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(rho, comp),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position),
				CHF_CONST_REAL(a_mass));
	break;
      }
    case DepositionType::CIC:
      {
	FORT_CIC_DEPOSIT_SCALAR2(CHF_FRA1(rho, comp),
				 CHF_CONST_REALVECT(a_probLo),
				 CHF_CONST_REALVECT(a_dx),
				 CHF_CONST_REALVECT(a_position),
				 CHF_CONST_REAL(a_mass));
	break;
      }
    default:
      {
	MayDay::Error("EBParticleMesh::depositParticle2 - Invalid deposition type - only NGP and CIC supported for this deposition method.");
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
  
  CH_assert(a_rho.nComp() == 0);  

  // TLDR: This performs regular deposition as if the particle lives on regular mesh data. If the cell is irregular we can use a class option
  //       to enforce NGP deposition in those cells. If the particle lives in a multi-valued cell I have no idea how to handle deposition.
  //
  //       Note that this is the version which deposits with particle widths that are 4 times the "usual" width. Currently, only NGP
  //       and CIC is supported. 

  constexpr int comp = 0;
  
  const RealVect rv  = (a_position - a_probLo)/a_dx;
  const IntVect  iv  = IntVect(D_DECL(floor(rv[0]), floor(rv[1]), floor(rv[2])));

  FArrayBox& rho = a_rho.getFArrayBox();    

  // We can force NGP in cut-cells if we want. 
  if(m_ebisbox.isIrregular(iv) && a_forceIrregNGP){
    FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(rho, comp),
			    CHF_CONST_REALVECT(a_probLo),
			    CHF_CONST_REALVECT(a_dx),
			    CHF_CONST_REALVECT(a_position),
			    CHF_CONST_REAL(a_mass));
  }
  else{
    switch (a_depositionType){
    case DepositionType::NGP:
      {
	FORT_NGP_DEPOSIT_SCALAR(CHF_FRA1(rho, comp),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position),
				CHF_CONST_REAL(a_mass));
	break;
      }
    case DepositionType::CIC:
      {
	FORT_CIC_DEPOSIT_SCALAR4(CHF_FRA1(rho, comp),
				 CHF_CONST_REALVECT(a_probLo),
				 CHF_CONST_REALVECT(a_dx),
				 CHF_CONST_REALVECT(a_position),
				 CHF_CONST_REAL(a_mass));
	break;
      }
    default:
      {
	MayDay::Error("EBParticleMesh::depositParticle4 - Invalid deposition type - only NGP and CIC supported for this deposition method.");
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

  const FArrayBox& meshField = a_meshField.getFArrayBox();      

  // Irregular cells always do an NGP deposit to prevent clouds leaking into the other side. 
  if(m_ebisbox.isIrregular(iv) && a_forceIrregNGP){
    FORT_NGP_INTERPOLATE_SCALAR(CHF_REAL(a_particleField),
				CHF_CONST_FRA1(meshField,comp),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position));
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

  const FArrayBox& meshField = a_meshField.getFArrayBox();

  // Irregular cells always do an NGP deposit to prevent clouds leaking into the other side.
  if(m_ebisbox.isIrregular(iv) && a_forceIrregNGP){    
    FORT_NGP_INTERPOLATE_VECTOR(CHF_REALVECT(a_particleField),
				CHF_CONST_FRA(meshField),
				CHF_CONST_REALVECT(a_probLo),
				CHF_CONST_REALVECT(a_dx),
				CHF_CONST_REALVECT(a_position));
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
