#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include "PetscCompGridVTO.H"
#include "FluxBox.H"
#include "IntVectSet.H"
#include "NamespaceHeader.H"

// derived ViscousTensor Operator class

#ifdef CH_USE_PETSC
void
PetscCompGridVTO::clean()
{
  PetscCompGrid::clean();
}

void 
PetscCompGridVTO::createOpStencil(IntVect a_iv, int a_ilev,const DataIndex &a_di, StencilTensor &a_sten)
{
  CH_TIME("PetscCompGridVTO::createOpStencil");
  Real dx=m_dxs[a_ilev][0],idx2=1./(dx*dx);
  Real eta_x0,eta_x1,lam_x0,lam_x1,eta_y0,eta_y1,lam_y0,lam_y1;
  IntVect ivwst=IntVect::Zero,ivest=IntVect::Zero,ivsth=IntVect::Zero,ivnth=IntVect::Zero;
  //IntVect ivSW(-1,-1),ivSE(1,-1),ivNE(1,1),ivNW(-1,1);
  IntVect ivSW(IntVect::Unit),ivSE(IntVect::Unit),ivNE(IntVect::Unit),ivNW(IntVect::Unit);
  ivSW[1] = ivSW[0] = -1; // >= 2D
  ivSE[1] = -1; // >= 2D
  ivNW[0] = -1;

  ivwst[0] = -1; ivest[0] = 1; ivsth[1] = -1; ivnth[1] = 1;
  {    // get stencil coeficients
    Real beta_hinv2=m_beta*idx2;
    const FluxBox &etaFab = (*m_eta[a_ilev])[a_di];
    const FArrayBox &eta_x = etaFab[0];
    const FArrayBox &eta_y = etaFab[1];
    const FluxBox &lamFab = (*m_lamb[a_ilev])[a_di];
    const FArrayBox &lam_x = lamFab[0];
    const FArrayBox &lam_y = lamFab[1];
    eta_x0 = beta_hinv2*eta_x(a_iv,0);           eta_y0 = beta_hinv2*eta_y(a_iv,0);
    eta_x1 = beta_hinv2*eta_x(ivest+a_iv,0);     eta_y1 = beta_hinv2*eta_y(ivnth+a_iv,0);
    lam_x0 = beta_hinv2*lam_x(a_iv,0);           lam_y0 = beta_hinv2*lam_y(a_iv,0);
    lam_x1 = beta_hinv2*lam_x(ivest+a_iv,0);     lam_y1 = beta_hinv2*lam_y(ivnth+a_iv,0);
  }

  // loop over two equations, should probably just hard wire this
  for (int kk=0,vv=1;kk<2;kk++,vv--)
    {
      Real vd;
      // add one eta everywhere -- same for u and v
      vd = eta_x0;
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivwst+a_iv,a_ilev)];  v1.define(CH_SPACEDIM);
        v1.addValue(kk,kk,vd);
      }
      vd = eta_x1;
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivest+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,kk,vd);
      }
      vd = eta_y0;
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivsth+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,kk,vd);
      }
      vd = eta_y1;
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivnth+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,kk,vd);
      }
      vd = -(eta_x0 + eta_y0 + eta_x1 + eta_y1);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,kk,vd);
      }

      // add extra eta and the lambda for each direction
      if (kk==0)
        {
          vd = eta_x0 + lam_x0;
          { 
            StencilTensorValue &v1 = a_sten[IndexML(ivwst+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
            v1.addValue(kk,kk,vd);
          }
          vd = eta_x1 + lam_x1;
          { 
            StencilTensorValue &v1 = a_sten[IndexML(ivest+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
            v1.addValue(kk,kk,vd);
          }
          vd = -(eta_x0 + eta_x1 + lam_x0 + lam_x1);
          { 
            StencilTensorValue &v1 = a_sten[IndexML(a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
            v1.addValue(kk,kk,vd);
          }
        }
      else {
        vd = eta_y0 + lam_y0;
        { 
          StencilTensorValue &v1 = a_sten[IndexML(ivsth+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
          v1.addValue(kk,kk,vd);
        }
        vd = eta_y1 + lam_y1;
        { 
          StencilTensorValue &v1 = a_sten[IndexML(ivnth+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
          v1.addValue(kk,kk,vd);
        }
        vd = -(eta_y0 + eta_y1 + lam_y0 + lam_y1);
        { 
          StencilTensorValue &v1 = a_sten[IndexML(a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
          v1.addValue(kk,kk,vd);
        }
      }

      // u -- v coupling terms -- 8 point stencil
      // S
      if (kk==0) vd = 0.25*(-lam_x1 + lam_x0);
      else vd = 0.25*(-eta_x1 + eta_x0);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivsth+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // S.W.
      if (kk==0) vd = 0.25*(eta_y0 + lam_x0);
      else vd = 0.25*(eta_x0 + lam_y0);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivSW+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // S.E.
      if (kk==0) vd = 0.25*(-lam_x1 - eta_y0);
      else vd = 0.25*(-lam_y0 - eta_x1);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivSE+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // W
      if (kk==0) vd = 0.25*(-eta_y1 + eta_y0);
      else vd = 0.25*(-lam_y1 + lam_y0);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivwst+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // E
      if (kk==0) vd = 0.25*(eta_y1 - eta_y0);
      else vd = 0.25*(lam_y1 - lam_y0);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivest+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // N.W.
      if (kk==0) vd = 0.25*(-lam_x0 - eta_y1);
      else vd = 0.25*(-lam_y1 - eta_x0);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivNW+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // N
      if (kk==0) vd = 0.25*(lam_x1 - lam_x0);
      else vd = 0.25*(eta_x1 - eta_x0);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivnth+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // N.E.
      if (kk==0) vd = 0.25*(eta_y1 + lam_x1);
      else vd = 0.25*(eta_x1 + lam_y1);
      { 
        StencilTensorValue &v1 = a_sten[IndexML(ivNE+a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,vv,vd);
      }
      // diagonal alpha term
      vd = m_alpha * (!(m_a[a_ilev]) ? 1. : (*m_a[a_ilev])[a_di](a_iv,0));
      { 
        StencilTensorValue &v1 = a_sten[IndexML(a_iv,a_ilev)]; v1.define(CH_SPACEDIM);
        v1.addValue(kk,kk,vd);
      }
    } // kk=1:2 loop
  
  if (0)
    {
      Real summ=0.;
      StencilTensor::const_iterator end = a_sten.end(); 
      for ( StencilTensor::const_iterator it = a_sten.begin(); it != end; ++it) {
        for (int i=0; i<CH_SPACEDIM; ++i) 
          {
            for (int j=0; j<CH_SPACEDIM; ++j) 
              {
                summ += it->second.value(i,j);
              }
          }
      }
      
      if (abs(summ)>1.e-10)
        {
          for (StencilTensor::const_iterator it = a_sten.begin(); it != end; ++it) {
            pout() << it->first << "  ERROR: \n";
            for (int i=0; i<CH_SPACEDIM; ++i) 
              {
                pout() << "\t\t";
                for (int j=0; j<CH_SPACEDIM; ++j) 
                  {
                    pout() << it->second.value(i,j) << ", ";
                  }
                pout() << endl;
              }
          }
          pout() << "\t ERROR summ: " << summ << endl; 
        }
    }
}

// clean out BC from stencil of cell a_iv
void 
PetscCompGridVTO::applyBCs( IntVect a_iv, int a_ilev, const DataIndex &di_dummy, Box a_dombox, 
                            StencilTensor &a_sten )
{
  CH_TIME("PetscCompGridVTO::applyBCs");
  Vector<Vector<StencilNode> > new_vals; // StencilTensorValue
  
  RefCountedPtr<BCFunction> bc = m_bc.getBCFunction();
  CompGridVTOBC *mybc = dynamic_cast<CompGridVTOBC*>(&(*bc));
  
  // count degrees, add to 'corners'
  Vector<IntVectSet> corners(CH_SPACEDIM);
  StencilTensor::const_iterator end2 = a_sten.end(); 
  for (StencilTensor::const_iterator it = a_sten.begin(); it != end2; ++it) 
    {
      const IntVect &jiv = it->first.iv();
      if (!a_dombox.contains(jiv) && it->first.level()==a_ilev) // a BC
        {
          int degree = CH_SPACEDIM-1;
          for (int dir=0;dir<CH_SPACEDIM;++dir)
            {
              if (jiv[dir] >= a_dombox.smallEnd(dir) && jiv[dir] <= a_dombox.bigEnd(dir)) degree--;
            }
          CH_assert(degree>=0);
          corners[degree] |= jiv;
        }
    }
  // move ghosts starting at with high degree corners and cascade to lower degree
  for (int ideg=CH_SPACEDIM-1;ideg>=0;ideg--)
    {
      // iterate through list
      for (IVSIterator ivit(corners[ideg]); ivit.ok(); ++ivit)
        {
          const IntVect &jiv = ivit(); // get rid of this bc ghost
          // get layer of ghost
          Box gbox(IntVect::Zero,IntVect::Zero); gbox.shift(jiv);
          Box inter = a_dombox & gbox;    CH_assert(inter.numPts()==0);
          int igid = -1; // find which layer of ghost I am
          do{
            igid++;
            gbox.grow(1);
            inter = a_dombox & gbox;
          }while (inter.numPts()==0);
	  if (igid!=0) MayDay::Error("PetscCompGridVTO::applyBCs layer???");
          for (int dir=0;dir<CH_SPACEDIM;++dir)
            {
              if (jiv[dir] < a_dombox.smallEnd(dir) || jiv[dir] > a_dombox.bigEnd(dir))
		{ // have a BC, get coefs on demand
		  int iside = 1; // hi
		  int isign = 1;
                  if (jiv[dir] < a_dombox.smallEnd(dir)) isign = -1;
		  if (jiv[dir] < a_dombox.smallEnd(dir)) iside = 0; // lo
		  if (!mybc) MayDay::Error("PetscCompGridVTO::applyBCs: wrong BC!!!!");
		  new_vals.resize(1);
		  new_vals[0].resize(1); 
		  //new_vals[i][j].second.setValue(mybc->getCoef(j,i));
                  for (int comp=0; comp<SpaceDim; comp++)
                    {
		      new_vals[0][0].second.define(1);
                      if (mybc->m_bcDiri[dir][iside][comp]) new_vals[0][0].second.setValue(comp,comp,-1.); // simple diagonal 
                      else new_vals[0][0].second.setValue(comp,comp,1.);
                    } // end loop over components
                  new_vals[0][0].first.setLevel(a_ilev);
		  
		  IndexML kill(jiv,a_ilev);
		  IntVect biv = jiv;
		  biv.shift(dir,-isign*(igid+1));                 
		  new_vals[igid][0].first.setIV(biv);
		  if (ideg>0 && !a_dombox.contains(biv)) // this could be a new stencil value for high order BCs
		    {
		      corners[ideg-1] |= biv; // send down to lower list for later removal
		    }
		  else CH_assert(a_dombox.contains(biv));
		
		  StencilProject(kill,new_vals[igid],a_sten);
		  int nrm = a_sten.erase(kill); CH_assert(nrm==1);
		  break;
		} // BC
	    } // spacedim
	} // ghosts
    } // degree
}

#endif // ifdef petsc

// CompGridVTOBC
void 
CompGridVTOBC::createCoefs()
{  
  m_isDefined=true;
  // is this needed ?
  m_Rcoefs[0] = -1.;
}

// 
void 
CompGridVTOBC::operator()( FArrayBox&           a_state,
			   const Box&           a_valid,
			   const ProblemDomain& a_domain,
			   Real                 a_dx,
			   bool                 a_homogeneous)
{
  const Box& domainBox = a_domain.domainBox();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (!a_domain.isPeriodic(idir))
        {
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Side::LoHiSide side = sit();              
              if (a_valid.sideEnd(side)[idir] == domainBox.sideEnd(side)[idir])
                {
                  // Dirichlet BC
                  int isign = sign(side);
                  Box toRegion = adjCellBox(a_valid, idir, side, 1);
                  // include corner cells if possible by growing toRegion in transverse direction
                  toRegion.grow(1);
                  toRegion.grow(idir, -1);
                  toRegion &= a_state.box();
                  for (BoxIterator bit(toRegion); bit.ok(); ++bit)
                    {
                      IntVect ivTo = bit();                     
                      IntVect ivClose = ivTo - isign*BASISV(idir);
                      for (int ighost=0;ighost<m_nGhosts[0];ighost++,ivTo += isign*BASISV(idir))
                        {
                          //for (int icomp = 0; icomp < a_state.nComp(); icomp++) a_state(ivTo, icomp) = 0.0;
                          IntVect ivFrom = ivClose;
                          
                          // hardwire to linear BCs for now
                          for (int icomp = 0; icomp < a_state.nComp() ; icomp++)
                            {
                              if (m_bcDiri[idir][side][icomp]) 
                                {
                                  a_state(ivTo, icomp) = (-1.0)*a_state(ivFrom, icomp);
                                }
                              else 
                                {
                                  a_state(ivTo, icomp) = (1.0)*a_state(ivFrom, icomp);
                                }
                            }   
                        }
                    } // end loop over cells
                } // if ends match
            } // end loop over sides
        } // if not periodic in this direction
    } // end loop over directions    
}

#include "NamespaceFooter.H"

#ifdef CH_USE_PETSC

#endif
