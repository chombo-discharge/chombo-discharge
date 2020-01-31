#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "LoHiSide.H"
#include "BoxIterator.H"
#include "ConstantIBC.H"
#include "FourthOrderUtil.H"
#include "FourthOrderUtilF_F.H"
#include "FORT_PROTO.H"
//#include "BCUtilF_F.H"

#define OLDBC

#include "CoarseAverage.H"

// Null constructor
ConstantIBC::ConstantIBC()
{
   m_params_are_set = false;
}

// Sets parameters used by initial conditions
void ConstantIBC::setParams(Real a_mag)
{
   CH_assert(m_params_are_set == false);
   m_mag = a_mag;
   m_params_are_set = true;
}

/// set uniform velocity field
void
ConstantIBC::setUniformVel(const RealVect& a_vel)
{
   m_velType = UNIFORM;
   m_uniformVel = a_vel;
}

// set parameter for solid-body rotation
void
ConstantIBC::setSolidBodyRotation(const RealVect& a_rotationCenter,
                                  const Real a_omega)
{
  m_velType = SOLIDBODY;
  m_rotationCenter = a_rotationCenter;
  m_omega = a_omega;
}

// set parameters for translating-oscilating velocity
void
ConstantIBC::setTranslatingOscillation(const RealVect& a_translation_vel,
                                       const Real a_oscillation_amplitude)
{
  m_velType = TRANSOSC;
  m_uniformVel = a_translation_vel;
  m_oscAmp = a_oscillation_amplitude;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new BasicIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
BasicIBC* ConstantIBC::new_basicIBC()
{
  ConstantIBC* retval = new ConstantIBC();
  retval -> setParams( m_mag );
  //[NOTE: do this even though setParams() will set it true
  //       because it might not have been set in the source
  //       object (this).  Of course, that would be a bad idea
  //       because then any object created by the factor wont
  //       be usable.  Caveat usor.  -dbs]

  retval ->m_velType = m_velType;
  retval->m_uniformVel = m_uniformVel;
  retval->m_omega = m_omega;
  retval->m_oscAmp = m_oscAmp;
  retval->m_rotationCenter = m_rotationCenter;
  retval->m_params_are_set = m_params_are_set;

  return static_cast<BasicIBC*>(retval);
}

// Set up initial conditions
void ConstantIBC::initialize(LevelData<FArrayBox>& a_phiJ,
                             const ProblemDomain& a_domain,
                             const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                             Real a_dx,
                             Real a_time)
{
   CH_assert(m_params_are_set == true);
#if 1
   LevelData<FArrayBox> phiAvg( a_phiJ.getBoxes(),
                                a_phiJ.nComp(),
                                a_phiJ.ghostVect() + IntVect::Unit );
   DataIterator dit = phiAvg.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      pointVal( phiAvg[dit], a_domain, a_coordSys,
                phiAvg[dit].box(), false, a_dx, a_time );
   }

   // convert point values into 4th-order cell averages
   fourthOrderAverage(phiAvg);

   // construct cell average of phi*J
   const LevelData<FArrayBox>& cellAvgJ = a_coordSys.getJ();
   fourthOrderCellProd( a_phiJ, phiAvg, cellAvgJ );
#else
   bool temp = m_exactOnPhysicalDomain;
   setExactOnPhysicalDomain( false );
   exactSoln( a_phiJ, a_domain, a_coordSys, a_dx, a_time);
   setExactOnPhysicalDomain( temp );
#endif
}


void
ConstantIBC::advVel(LevelData<FArrayBox>& a_advVel,
                    const ProblemDomain& a_domain,
                    const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                    Real a_dXi,
                    Real a_time)
{
   DataIterator dit = a_advVel.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      if (m_velType == UNIFORM)
      {
         for (int dir=0; dir<SpaceDim; dir++)
         {
            a_advVel[dit].setVal(m_uniformVel[dir],dir);
         }
      }
      else if (m_velType == SOLIDBODY)
      {
         Real twoPiOmega = 8.0 * atan(1.0) * m_omega;
         FArrayBox& advVel = a_advVel[dit];

         BoxIterator bit(advVel.box());
         for (bit.begin();bit.ok();++bit)
         {
            // get cell center
            IntVect iv = bit();
            RealVect xi;
            for (int idir=0; idir<SpaceDim; idir++)
            {
               xi[idir] = ( 0.5 + iv[idir] ) * a_dXi;
            }
            RealVect x( a_coordSys.realCoord( xi ) - m_rotationCenter );
            D_TERM( advVel(iv,0) = -twoPiOmega * x[1];,
                    advVel(iv,1) =  twoPiOmega * x[0];,
                    advVel(iv,2) =  0; );
         }
      }
      else if (m_velType == TRANSOSC)
      {
         Real twoPi = 8.0 * atan(1.0);
         FArrayBox& advVel = a_advVel[dit];
         Real speed = sqrt( m_uniformVel.dotProduct( m_uniformVel ) );
         Real L = speed / m_uniformVel.dotProduct( RealVect::Unit );

         BoxIterator bit(advVel.box());
         for (bit.begin();bit.ok();++bit)
         {
            // get cell center
            IntVect iv = bit();
            RealVect xi;
            for (int idir=0; idir<SpaceDim; idir++)
            {
               xi[idir] = ( 0.5 + iv[idir] ) * a_dXi;
            }
            RealVect x( a_coordSys.realCoord( xi ) );
            Real zeta = x.dotProduct( m_uniformVel ) / ( L * speed );
            Real osc = m_oscAmp * cos( twoPi * zeta );
            D_TERM( advVel(iv,0) = L*(m_uniformVel[0]-m_uniformVel[1]*osc);,
                    advVel(iv,1) = L*(m_uniformVel[1]+m_uniformVel[0]*osc);,
                    advVel(iv,2) = 0; );
         }

      }
      else
      {
         MayDay::Error("ConstantBC::advVel -- bad velType");
      }
   }
   fourthOrderAverage(a_advVel);
}

/// fill ghost cell values at domain boundaries
void
ConstantIBC::ghostCellBC(LevelData<FArrayBox>& a_phi,
                         const ProblemDomain& a_domain,
                         const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                         Real a_dx,
                         Real a_time)
{
   IntVect ghostVect = a_phi.ghostVect();
   const DisjointBoxLayout& grids = a_phi.getBoxes();
   const Box& domainBox = a_domain.domainBox();

   LevelData<FArrayBox> advVel( grids, SpaceDim, ghostVect );
   this->advVel( advVel, a_domain, a_coordSys, a_dx, a_time);
   LevelData<FluxBox> faceVel( grids, SpaceDim, ghostVect );
   fourthOrderCellToFace( faceVel, advVel );
   faceVel.exchange();

   DataIterator dit = a_phi.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FArrayBox& thisPhi = a_phi[dit];
      //FluxBox& thisFaceVel = faceVel[dit];
      if (!a_domain.contains(thisPhi.box()) )
      {
         for (int dir=0; dir<SpaceDim; dir++)
         {
            // this is designed to ensure that
            // corner cells get filled
            IntVect tanGrow(ghostVect);
            tanGrow[dir] = 0;
            if (!a_domain.isPeriodic(dir))
            {
               SideIterator sit;
               for (sit.begin(); sit.ok(); ++sit)
               {
                  Box ghostBox;
                  if (sit() == Side::Lo)
                  {
                     ghostBox = adjCellLo(domainBox, dir, ghostVect[dir]);
                  }
                  else
                  {
                     ghostBox = adjCellHi(domainBox, dir, ghostVect[dir]);
                  }
                  // ensure that corner cells get filled
                  ghostBox.grow(tanGrow);
                  ghostBox &= thisPhi.box();

                  if (!ghostBox.isEmpty())
                  {
                     // need to grow the ghost box by
                     // one for 4th-order averaging
                     ghostBox.grow(1);
                     FArrayBox ghostData(ghostBox,
                                         thisPhi.nComp());
                     pointVal(ghostData,
                              a_domain,
                              a_coordSys,
                              ghostBox,
                              true, a_dx, a_time);

                     fourthOrderAverageCell(ghostData);
                     ghostBox.grow(-1);
#ifndef OLDBC
                     const int side = sit();
                     FORT_BILINEXTRAP(CHF_CONST_INT(dir),
                                      CHF_CONST_INT(side),
                                      CHF_BOX(ghostBox),
                                      CHF_CONST_INTVECT(ghostVect),
                                      CHF_CONST_FRA(thisFaceVel[dir]),
                                      CHF_CONST_FRA(thisPhi),
                                      CHF_FRA(ghostData));
#endif
                     // now copy into ghost cells
                     thisPhi.copy(ghostData, ghostBox);
                  }  // end if there are domain ghost cells here
               } // end loop over hi-lo
            } // end if not periodic in this direction
         } // end loop over directions
      } // end if phi sticks out of domain
   } // end loop over grids
}

static void computeFineCellAveragedJacobian(
   LevelData<FArrayBox>& a_fineJ,
   const CoordSys<FArrayBox,FluxBox>& a_coordSys,
   Real a_dx)
{
   LevelData<FArrayBox> Jpw( a_fineJ.getBoxes(),
                             a_fineJ.nComp(),
                             a_fineJ.ghostVect() + IntVect::Unit );
   DisjointBoxLayout finerGrids = a_fineJ.getBoxes();
   DataIterator dit = finerGrids.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      FArrayBox& thisJpw = Jpw[dit];
      BoxIterator bit(thisJpw.box());
      for (bit.begin();bit.ok();++bit)
      {
         IntVect iv = bit();
         RealVect xi( iv );
         xi += 0.5;
         xi *= a_dx;
         RealVect x( a_coordSys.realCoord( xi ) );
         thisJpw(iv) = a_coordSys.pointwiseJ( x );
      }
   }

   // convert point values to 4th-order cell averages
   fourthOrderAverage( Jpw );

   // copy into m_fineJ
   for (dit.begin(); dit.ok(); ++dit)
   {
      a_fineJ[dit].copy( Jpw[dit], a_fineJ[dit].box() );
   }
}

/// compute exact solution
/**
 */
void
ConstantIBC::exactSoln(
   LevelData<FArrayBox>& a_phiJ,
   const ProblemDomain& a_domain,
   const CoordSys<FArrayBox,FluxBox>& a_coordSys,
   Real a_dx,
   Real a_time)
{
   // fine grid domain
   const int nRef = 8;
   const DisjointBoxLayout& grids = a_phiJ.getBoxes();
   DisjointBoxLayout finerGrids;
   refine( finerGrids, grids, nRef );
   const Real fineDx = a_dx / nRef;
   ProblemDomain finerDomain( a_domain );
   finerDomain.refine( nRef );

   // fine grid data objects
   const int ncomp( a_phiJ.nComp() );
   const IntVect ghostVect( a_phiJ.ghostVect() );
   LevelData<FArrayBox> finerPhi( finerGrids, ncomp, ghostVect );
   LevelData<FArrayBox> finerPhiJ( finerGrids, ncomp, ghostVect );

   // compute 4th-order cell-averaged phi solution on finer computatonal mesh
   DataIterator dit = finerPhi.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      pointVal(finerPhi[dit], finerDomain, a_coordSys,
               finerPhi[dit].box(), false, fineDx, a_time );
   }
   fourthOrderAverage( finerPhi );

   // get the cell-averaged Jacobian on the finer computational mesh
   if ( !m_fineJ.isDefined() )
   {
      m_fineJ.define( finerGrids, 1, ghostVect );
      computeFineCellAveragedJacobian( m_fineJ, a_coordSys, fineDx );
      m_crseJ.define( grids, 1, ghostVect );
      CoarseAverage averager( finerGrids, grids, 1, nRef );
      averager.averageToCoarse( m_crseJ, m_fineJ );
   }

   // get the cell-averaged phi*J on the finer computational mesh
   fourthOrderCellProd( finerPhiJ, finerPhi, m_fineJ );

   // average down
   CoarseAverage averager( finerGrids, grids, ncomp, nRef );
   averager.averageToCoarse( a_phiJ, finerPhiJ );

   // get the cell-averaged phi on the finer physical mesh, if necessary
   if (m_exactOnPhysicalDomain)
   {
      DataIterator dit = a_phiJ.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
      {
         a_phiJ[dit].divide( m_crseJ[dit], grids[dit], 0, 0, ncomp );
      }
   }

   return;
}


// Set up initial conditions
void ConstantIBC::pointVal(FArrayBox& a_phi,
                           const ProblemDomain& a_domain,
                           const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                           const Box& a_box,
                           bool a_includeJ,
                           Real a_dXi,
                           Real a_time)
{
   CH_assert(m_params_are_set == true);

   // compute initial values
   a_phi.setVal( m_mag );

   if (a_includeJ)
   {
     // 0.5 for cell-centered, 0 for node-centering
     RealVect meshOffset(0.5*RealVect::Unit);
     for (int dir=0; dir<SpaceDim; dir++)
       {
         if (a_box.type(dir) == IndexType::NODE) meshOffset[dir] = 0.0;
       }

      BoxIterator bit( a_box );
      for (bit.begin();bit.ok();++bit)
      {
         IntVect iv = bit();
         RealVect Xi( iv );
         Xi += meshOffset;
         Xi *= a_dXi;
         RealVect X( a_coordSys.realCoord( Xi ) );
         Real J = a_coordSys.pointwiseJ( X );

         for (int icomp=0; icomp<a_phi.nComp(); icomp++)
         {
            a_phi(iv,icomp) *= J;
         }
      }
   }
}

