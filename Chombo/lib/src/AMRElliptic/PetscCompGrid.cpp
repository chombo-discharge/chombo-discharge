#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PetscCompGrid.H"
#include "IntVectSet.H"
#include "CoarseAverage.H"
#include "SPMD.H"

#include "NamespaceHeader.H"


std::ostream& operator<< (std::ostream& os, GID_type a_type)
{
  switch (a_type)
    {
    case GHOST: os << "ghost";                       break;
    case FINE_COVERED: os << "covered (with fine)" ; break;
    case DIRI_BC: os << "BC (Diri)";                 break;
    case NEUM_BC: os << "BC (Neum)";                 break;
    case ANY_BC: os << "BC";                         break;
    case UNKNOWN: os << "unknown";                   break;
    default: os << "index " << (int)a_type;          break;
    }
  return os;
}

//
// My BC function
//
CompBC::~CompBC()
{
#ifdef CH_USE_PETSC
  PetscFree(m_Rcoefs);
#endif
}

CompBC::CompBC(int a_order, IntVect a_nGhosts) 
{
#ifdef CH_USE_PETSC
  m_Rcoefs = 0;
#endif

  m_isDefined = false;
  define(a_order,a_nGhosts);
}

void 
CompBC::define(int a_order, IntVect a_nGhosts)
{  
#ifdef CH_USE_PETSC
  if (m_Rcoefs) PetscFree(m_Rcoefs);
#endif
  if (a_nGhosts[0]!=1 && a_nGhosts[0]!=2)
    MayDay::Error("Unsupported number of ghosts in CompBC");
  
  m_nGhosts = a_nGhosts; m_nSources = a_order;

#ifdef CH_USE_PETSC
  PetscMalloc1(m_nSources*m_nGhosts[0],&m_Rcoefs);
#else
  m_Rcoefs.resize(m_nSources*m_nGhosts[0]);
#endif

  m_isDefined = false; // needs to be set
}

#ifdef CH_USE_PETSC
PetscReal CompBC::getCoef(int a_iSrc, int a_iGhost)
#else
Real CompBC::getCoef(int a_iSrc, int a_iGhost)
#endif
{
  if (!m_isDefined) createCoefs();
  return m_Rcoefs[a_iGhost*m_nSources + a_iSrc];
}

//
void 
ConstDiriBC::createCoefs()
{  
  m_isDefined=true;
  if (m_nSources==1) 
    {
      m_Rcoefs[0] = -1.;
      if (m_nGhosts[0]==2) m_Rcoefs[1] = -3.;   
    }
  else if (m_nSources==2) 
    { // s =  (8.0/3.0)*a_inhomogVal + (1.0/3.0)*(a_farVal) -2*(a_nearVal);
      m_Rcoefs[0] = -2.; m_Rcoefs[1] = 1./3.;
      if (m_nGhosts[0]==2) 
        {
          m_Rcoefs[2] = -21./2; m_Rcoefs[3] = 5./2;
        }
    }
  else if (m_nSources==3) 
    { // s = 12, 48
      m_Rcoefs[0] = -13./3; m_Rcoefs[1] = 5./3; m_Rcoefs[2] = -1./3;
      if (m_nGhosts[0]==2) 
        {
          m_Rcoefs[3] = -70./3; m_Rcoefs[4] = 32./3;  m_Rcoefs[5] = -7./3; 
        }
    }
  else if (m_nSources==4) 
    { // s = 60, 300
      m_Rcoefs[0] = -77./12; m_Rcoefs[1] = 43./12; m_Rcoefs[2] = -17./12; m_Rcoefs[3] = 3./12;
      if (m_nGhosts[0]==2) 
        {
          m_Rcoefs[4] = -505./12; m_Rcoefs[5] = 335./12;  m_Rcoefs[6] = -145./12; m_Rcoefs[7] = 27./12; 
        }
    }
  else
    MayDay::Error("Unsupported degree in ConstDiriBC");
}

void 
ConstDiriBC::operator()( FArrayBox&           a_state,
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
                  toRegion &= a_state.box();
                  for (BoxIterator bit(toRegion); bit.ok(); ++bit)
                    {
                      IntVect ivTo = bit();                     
                      IntVect ivClose = ivTo - isign*BASISV(idir);
                      for (int ighost=0;ighost<m_nGhosts[0];ighost++,ivTo += isign*BASISV(idir))
                        {
                          for (int icomp = 0; icomp < a_state.nComp(); icomp++) a_state(ivTo, icomp) = 0.0;
                          IntVect ivFrom = ivClose;
                          for (int i=0;i<m_nSources;i++,ivFrom -= isign*BASISV(idir))
                            {
                              for (int icomp = 0; icomp < a_state.nComp(); icomp++)
                                {
                                  a_state(ivTo, icomp) += m_Rcoefs[ighost*m_nSources + i]*a_state(ivFrom, icomp);
                                }
                            }
                        }
                    }
                } // if ends match
            } // end loop over sides
        } // if not periodic in this direction
    } // end loop over directions    
}

//
// PetscCompGrid
//

#ifdef CH_USE_PETSC
void
PetscCompGrid::clean()
{
  if (m_mat) 
    {
      MatDestroy(&m_mat);
      m_mat = 0; // this can get called multiple times
      if (m_domains.size() > 1)
        {
          const Box& stencilBox = m_FCStencils.box();
          for (BoxIterator bit(stencilBox); bit.ok(); ++bit)
            {
              IntVect offset = bit();
              delete m_FCStencils(offset, 0);
            }   
        }
       if (m_Pmat) 
         {
           MatDestroy(&m_Pmat);
           m_Pmat = 0;
         }
       if (m_from_petscscat) 
         {
           VecScatterDestroy(&m_from_petscscat);
           VecDestroy(&m_origvec);
           CH_assert(!m_from_petscscat);
         }
    }
  m_patch_size = 0;
  m_gid0 = 0;
  m_num_extra_nnz = 0;
}

//
PetscCompGrid::~PetscCompGrid()
{
  clean();
}

// helper function
void
PetscCompGrid::setCFCoverMaps( int a_nlev )
{
  // construct covering gids maps
  for (int ilev=1;ilev<a_nlev;ilev++)
    {
      Interval inter = m_GIDs[ilev]->interval();
      DisjointBoxLayout crsenedFineDBL,refinedCrsDBL;
      LevelData<BaseFab<PetscInt> > *pl;
      int refRatio = m_refRatios[ilev-1];

      // form m_crsSupportGIDs[ilev-1]:
      const DisjointBoxLayout& fdbl = m_grids[ilev];
      coarsen(crsenedFineDBL,fdbl,refRatio);
      // rad+1 needed to get stencils coarse cell coarse cover
      const IntVect nCrsSuppCells = (m_CFStencilRad+1)*IntVect::Unit;
      pl = new LevelData<BaseFab<PetscInt> >(crsenedFineDBL,1,nCrsSuppCells);
      m_crsSupportGIDs[ilev-1] = RefCountedPtr<LevelData<BaseFab<PetscInt> > >(pl);
      // set coarse gids == UNKNOWN
      for (DataIterator dit = fdbl.dataIterator(); dit.ok(); ++dit)
        {
          BaseFab<PetscInt>& gidfab = (*pl)[dit];
          gidfab.setVal(UNKNOWN); // enum this.  BCs and covered by this 
        }
      Copier rcopier(m_grids[ilev-1],crsenedFineDBL,nCrsSuppCells);
      m_GIDs[ilev-1]->copyTo(inter,*m_crsSupportGIDs[ilev-1],inter,rcopier);

      // from m_fineCoverGIDs[ilev]: fine support of stencils for coarse grid ilev-1
      const DisjointBoxLayout& cdbl = m_grids[ilev-1];
      refine(refinedCrsDBL,cdbl,refRatio);
      //
      if (getGhostVect()[0]>refRatio)
        {
          MayDay::Error("PetscCompGrid::define getGhostVect>refRatio - not setup for this");
        }
      // ghost cells for fine support governed by CF radius (12=nr*nr(rad+1), could use 1 with nesting radius==4 maybe
      const IntVect nFineProcGhosts = refRatio*refRatio*(m_CFStencilRad+1)*IntVect::Unit;
      pl = new LevelData<BaseFab<PetscInt> >(refinedCrsDBL,1,nFineProcGhosts); 
      m_fineCoverGIDs[ilev] = RefCountedPtr<LevelData<BaseFab<PetscInt> > >(pl);
      // set gids == UNKNOWN
      for (DataIterator dit = cdbl.dataIterator(); dit.ok(); ++dit)
        {
          BaseFab<PetscInt>& gidfab = (*pl)[dit];
          gidfab.setVal(UNKNOWN); // enum this. fine grid that is not there (true data here)
        }
      Copier pcopier(fdbl,refinedCrsDBL,nFineProcGhosts);
      m_GIDs[ilev]->copyTo(inter,*m_fineCoverGIDs[ilev],inter,pcopier);
      if (m_verbose>3) 
        {
          pout() << "PetscCompGrid::define: m_fineCoverGIDs["<< ilev-1 <<"] boxes: " << endl; 
          for (DataIterator dit = cdbl.dataIterator(); dit.ok(); ++dit)
            {
              BaseFab<PetscInt>& gidfab = (*pl)[dit];
              if (m_verbose>5)
                {
                  pout() << "PetscCompGrid::define: refined box="
                         << gidfab.box() 
                         << ", coarse box = " << m_grids[ilev-1][dit] << endl;
                }
            }
        }
    }
}

// Petsc composite grid solver - a matrix with solve methods that builds itself with a Chombo operator and a hierarchy of grids.
void
PetscCompGrid::define( const ProblemDomain &a_cdomain,
                       Vector<DisjointBoxLayout> &a_grids, 
                       Vector<int> &a_refRatios, 
                       BCHolder a_bc,
                       const RealVect &a_cdx,
                       int a_numLevels/* =-1 */, int a_ibase/* =0 */)
{
  CH_TIME("PetscCompGrid::define");
  PetscMPIInt nprocs,nblockpts=0;
  if (a_numLevels<=0) nprocs = a_grids.size() - 1;
  else nprocs = a_numLevels - 1;
  const int numLevs = nprocs - a_ibase + 1;

  m_bc = a_bc;

  PetscCompGrid::clean(); // this is virtual so lets not step on derived classes data

#ifndef CH_MPI
  m_repartition = PETSC_FALSE;
  nprocs = 1;
#else
  MPI_Comm_size(Chombo_MPI::comm,&nprocs);
#endif
    
  if (m_verbose>5) 
    {
      pout() << "PetscCompGrid::define: a_numLevels=" << a_numLevels << ", numLevs=" << 
        numLevs << ", a_ibase=" << a_ibase << endl; 
    }
  if (numLevs>1)
    {
      int degree = 3, nref = a_refRatios[0]; // assume same coarsening for all directions and levels
      IntVect interpUnit = IntVect::Unit;
      Box stencilBox( -m_CFStencilRad*interpUnit,
                      m_CFStencilRad*interpUnit );
      m_FCStencils.define(stencilBox,1);
      for (BoxIterator bit(stencilBox); bit.ok(); ++bit)
        {
          IntVect offset = bit();
          m_FCStencils(offset, 0) =
            new FourthOrderInterpStencil(offset, nref, degree);
        }
    }
  
  // copy in data
  m_domains.resize(numLevs);
  m_grids.resize(numLevs);
  m_refRatios.resize(numLevs);
  m_GIDs.resize(numLevs);
  m_crsSupportGIDs.resize(numLevs);
  m_fineCoverGIDs.resize(numLevs);
  m_dxs.resize(numLevs);

  RealVect dx = a_cdx;
  ProblemDomain dom = a_cdomain;
  for (int ilev=0; ilev<a_ibase ; ilev++) 
    {
      dom.refine(a_refRatios[ilev]);
      dx /= a_refRatios[ilev];
    }
  for (int ilev=a_ibase, ii=0; ilev < numLevs ; ilev++, ii++)
    {
      m_domains[ii] = dom;
      m_grids[ii] = a_grids[ilev];
      m_dxs[ii] = dx;
      if (ilev != numLevs-1) 
        {
          m_refRatios[ii] = a_refRatios[ilev];
          dom.refine(a_refRatios[ilev]);
          dx /= a_refRatios[ilev];
        }
      if (m_verbose>5) 
        {
          pout() << "PetscCompGrid::define: level=" << ilev << ", " << 
            m_grids[ii].size() << " patches" << endl; 
        }
    }

  // iterate over uncovered cells, count real cells and patches
  PetscInt my0s[2] = {0,0}; // count of gross cells, and non-empty patches
  for (int ilev=0;ilev<numLevs;ilev++)
    { 
      const DisjointBoxLayout& dbl = m_grids[ilev];
      // 3 is a hack for treb (?); refRat*(rad+1) happens at corners (not clear def of "radius")
      IntVect nProcessGhost = (ilev==0) ? /*getGhostVect()*/ 3*IntVect::Unit : 
        m_refRatios[ilev-1]*(m_CFStencilRad+1)*IntVect::Unit;
      m_GIDs[ilev] = RefCountedPtr<LevelData<BaseFab<PetscInt> > >
        (new LevelData<BaseFab<PetscInt> >(dbl,1,nProcessGhost));
      // set gids == GHOST
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          gidfab.setVal(GHOST);
        }

      // zero covered if not finest
      if (ilev!=numLevs-1) 
        {
          // zero out fine cover
          DisjointBoxLayout dblCoarsenedFine;
          Copier copier;
          coarsen(dblCoarsenedFine, m_grids[ilev+1], m_refRatios[ilev]); // coarsens entire grid
          copier.define(dblCoarsenedFine, dbl, IntVect::Zero);
          for (CopyIterator it(copier, CopyIterator::LOCAL); it.ok(); ++it)
            {
              const MotionItem& item = it();
              (*m_GIDs[ilev])[item.toIndex].setVal(FINE_COVERED, item.toRegion, 0);           
            }
          for (CopyIterator it(copier, CopyIterator::TO); it.ok(); ++it)
            {
              const MotionItem& item = it();
              (*m_GIDs[ilev])[item.toIndex].setVal(FINE_COVERED, item.toRegion, 0);
            }
        }
      
      // add extra covered
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit]; // make explicit that we modify this
          addExtraCovered(FINE_COVERED,ilev,dit(),gidfab);
        }
      // count
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          Box region = dbl[dit]; // no ghosts
          BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          int nreal = 0; // count real
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit(); 
              if (gidfab(iv,0) == GHOST){ nreal++; my0s[0]++; }// not covered, so real
            }
          if (m_repartition && nreal) 
            {
              my0s[1]++; 
              if (!nblockpts) nblockpts = nreal;
              else if (nblockpts != nreal) MayDay::Error("PetscCompGrid::define block size changed -- must use a 'tree' grid");
            }
        }
    } // levels
  if (nblockpts) 
    {
      m_patch_size = (PetscInt)pow((double)nblockpts,(double)1.0/(double)SpaceDim);
      if (m_verbose>0) pout() << "\tPetscCompGrid::define m_patch_size : " <<  m_patch_size << ", nblockpts : " << nblockpts << endl;
    }

  // get global indices - scan
  m_nlocrealpatches = my0s[1]; // output for partitioning
#ifdef CH_MPI
  MPI_Comm wcomm = Chombo_MPI::comm;
  PetscInt result[2];
  MPI_Datatype mtype;
  PetscDataTypeToMPIDataType(PETSC_INT,&mtype);
  MPI_Scan(my0s, result, 2, mtype, MPI_SUM, Chombo_MPI::comm);
  m_gid0 = result[0] - my0s[0];
  m_patchid0 = result[1] - my0s[1]; // output for partitioning
#else
  MPI_Comm wcomm = PETSC_COMM_SELF; // not used really
  m_gid0 = 0; m_patchid0 = 0;
#endif

  // set GIDs & exchange to get process ghost ids
  for (PetscInt ilev=0,gid=m_gid0;ilev<numLevs;ilev++)
    {
      const DisjointBoxLayout& dbl = m_grids[ilev];
      // count and set ids
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          Box region = dbl[dit]; // no ghosts
          BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit(); 
              if (gidfab(iv,0) == GHOST) gidfab(iv,0) = gid++; // not covered, so real
              else CH_assert(gidfab(iv,0) == FINE_COVERED);   // fine covered
            }
        }
      m_GIDs[ilev]->exchange(); // this has to be virtualized for COGENT (or does it ignore block patch ghosts?)
    }

  // set helper GID maps
  setCFCoverMaps(numLevs);

  // create glabal matrix here (we have local size)
  MatCreate(wcomm,&m_mat);
  MatSetSizes(m_mat,my0s[0]*m_dof,my0s[0]*m_dof,PETSC_DECIDE,PETSC_DECIDE);
  MatSetBlockSize(m_mat,m_dof);
  MatSetType(m_mat,MATAIJ);
  MatSetFromOptions(m_mat);
}

// main method of this class -- make a matrix
#undef __FUNCT__
#define __FUNCT__ "createMatrix"
PetscErrorCode
PetscCompGrid::createMatrix(int a_makePmat/*=0*/)
{
  CH_TIME("PetscCompGrid::createMatrix");
  PetscErrorCode ierr;
  PetscInt nloc,nglob,mygidnext,ilev,idx,patchidx,max_size;
  int nGrids=m_domains.size(),nblockpts=pow((Real)m_patch_size,SpaceDim);
  IntVect nghost = getGhostVect();
  Vector<StencilTensor> stenVect;
  Vector<StencilTensor> patchStencil;
  PetscInt *d_nnz, *o_nnz;
  PetscFunctionBeginUser;

  // set preallocation
#if defined(PETSC_USE_LOG)
  PetscLogEventBegin(m_event0,0,0,0,0);
#endif
  ierr = MatGetLocalSize(m_mat,&nloc,0);CHKERRQ(ierr);
  nloc /= m_dof;
  //ierr = MatGetSize(m_mat,&nglob,0);CHKERRQ(ierr); -- PETSc has not done the reduction yet!
  //nglob /= m_dof;
#ifdef CH_MPI
  MPI_Datatype mtype;
  PetscDataTypeToMPIDataType(PETSC_INT,&mtype);
  ierr = MPI_Allreduce(&nloc, &nglob, 1, mtype, MPI_SUM, Chombo_MPI::comm);CHKERRQ(ierr);
#else
  nglob = nloc;
#endif

  /* count nnz */
  PetscMalloc1((nloc+1)*m_dof, &d_nnz);
  PetscMalloc1((nloc+1)*m_dof, &o_nnz);
  stenVect.resize(nloc);
  if (m_repartition)
    {
      mygidnext = m_gid0 + nloc; 
      patchStencil.resize(m_nlocrealpatches);
    }
  // 
  if (a_makePmat)
    {
      ierr = MatDuplicate(m_mat,MAT_DO_NOT_COPY_VALUES,&m_Pmat);CHKERRQ(ierr);
    }
  // add data
  for (ilev=0,idx=0,patchidx=0;ilev<nGrids;ilev++)
    {
      const DisjointBoxLayout& dbl = m_grids[ilev];
      Box testbox,dombox(m_domains[ilev].domainBox());
      for (int dir=0; dir<CH_SPACEDIM; ++dir) 
        {
          if (m_domains[ilev].isPeriodic(dir)) {
            dombox.grow(dir,nghost[dir]); // nullify
          }
        }
      testbox = dombox;
      testbox.grow(-nghost); // non-boundary cells
      max_size=0;
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        { 
          BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          Box region = dbl[dit];
          int nreal = 0;
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if (gidfab(iv,0) >= 0) // not covered
                { 
                  nreal++;
                  // doit 
                  StencilTensor &sten2 = stenVect[idx];
                  createOpStencil(iv,ilev,dit(),sten2); // get raw stencil from app
                  if (!testbox.contains(iv)) // test is just for speed
                    {
                      applyBCs(iv,ilev,dit(),dombox,sten2);
                    }
                  if (ilev != 0) // test is just for speed
                    {
                      InterpToCoarse(iv,ilev,dit(),sten2);
                    }
                  if (ilev != nGrids-1) // test is just for speed
                    {
                      InterpToFine(iv,ilev,dit(),sten2);
                    }
                  // COGENT -- add function to transform block patch ghosts

                  // collect prealloc sizes
                  o_nnz[idx*m_dof] = d_nnz[idx*m_dof] = sten2.size();
                  if (o_nnz[idx*m_dof]>max_size) { max_size = o_nnz[idx*m_dof];}
                  d_nnz[idx*m_dof] += m_num_extra_nnz;
                  o_nnz[idx*m_dof] += m_num_extra_nnz;
                  if (d_nnz[idx*m_dof]>nloc) d_nnz[idx*m_dof] = nloc;
                  if (o_nnz[idx*m_dof]>(nglob-nloc)) o_nnz[idx*m_dof]=(nglob-nloc);
                  d_nnz[idx*m_dof] *= m_dof;
                  o_nnz[idx*m_dof] *= m_dof;
                  for (int i=1; i<m_dof; ++i) 
                    {
                      d_nnz[idx*m_dof+i] = d_nnz[idx*m_dof];
                      o_nnz[idx*m_dof+i] = o_nnz[idx*m_dof];
                    }
                  idx++;
                  // make graph for repartitioning
                  if (m_repartition)
                    {
                      BaseFab<PetscInt>& this_gidfabJ = (*m_GIDs[ilev])[dit()];
                      StencilTensor::const_iterator end = sten2.end(); 
                      for (StencilTensor::const_iterator it = sten2.begin(); it != end; ++it) 
                        {
                          const IndexML &ivJ = it->first;
                          int jLevel = ivJ.level();
                          BaseFab<PetscInt>& gidfabJ = (jLevel==ilev) ? this_gidfabJ : 
                            (jLevel==ilev+1) ? (*m_fineCoverGIDs[ilev+1])[dit()] : 
                            (*m_crsSupportGIDs[ilev-1])[dit()];
                          
                          if (!gidfabJ.box().contains(ivJ.iv())){
                            pout() << "ERROR row  " << IndexML(iv,ilev) << ", this box:" << 
                              m_grids[ilev][dit()] << ", fine box (failure): " << gidfabJ.box() << 
                              " failed to find "<< ivJ << endl;
                            MayDay::Error("PetscCompGrid::createMatrix failed to find cell");
                          }

                          PetscInt gidj = gidfabJ(ivJ.iv(),0)*m_dof;
                          if (gidj<0) 
                            {
                              pout() << "\tFAILED TO FIND iv=" << ivJ << "\t gidj type:" << (GID_type)gidj << "\tiLevel=" << ilev << "\tiv=" << iv << endl;
                              MayDay::Error("PetscCompGrid::createMatrix failed to find gid");
                            }
                          else if (gidj < m_gid0 || gidj >= mygidnext) // off proc value
                            {
                              StencilTensor &sten = patchStencil[patchidx]; // add to this row
                              // find patch that has this gidj
                              IntVect iv = ivJ.iv();
                              iv.coarsen(m_patch_size);          // get into block space
                              StencilTensorValue &v0 = sten[IndexML(iv,ivJ.level())]; // will create if not there
                              if (v0.value(0,0)==0.0) v0.setValue((Real)(gidj/nblockpts+1)); // put global index into stencil -- one based
                              else CH_assert((PetscInt)v0.value(0,0)==(gidj/nblockpts)+1);
                            }
                        }                     
                    } // repart
                } // real cell 
            } // box
          if (nreal && m_repartition)
            {
              if (m_verbose>2) pout() << patchidx << "\tPetscCompGrid::createMatrix num graph edges for repartitioning : " <<  patchStencil[patchidx].size() << endl;
              patchidx++;
            }
        } // level

      if (m_verbose>1) 
        {
          if (m_verbose>=5) 
            {
#ifdef CH_MPI
              MPI_Comm wcomm = Chombo_MPI::comm;
              PetscInt n = max_size;
              MPI_Datatype mtype;
              PetscDataTypeToMPIDataType(PETSC_INT,&mtype);
              ierr = MPI_Allreduce(&n, &max_size, 1, mtype, MPI_MAX, wcomm);CHKERRQ(ierr);
#endif
            }
          pout() << "\t PetscCompGrid::createMatrix level " << ilev+1 << 
            "/" << nGrids << ". domain " << m_domains[ilev] 
                 << ". max. stencil size: " 
                 << max_size << endl;
        }
    } // grids
  CH_assert(idx==nloc);  
  if (m_verbose>0 && m_repartition)
    {
      pout() << "\tPetscCompGrid::createMatrix: repartitioning with " << patchidx << " local patches ("<< m_nlocrealpatches << ")" << endl;
      CH_assert(patchidx==m_nlocrealpatches);
    }

  ierr = MatSeqAIJSetPreallocation(m_mat,0,d_nnz);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(m_mat,0,d_nnz,0,o_nnz);CHKERRQ(ierr);
  if (a_makePmat)
    {
      ierr = MatSeqAIJSetPreallocation(m_Pmat,0,d_nnz);CHKERRQ(ierr);
      ierr = MatMPIAIJSetPreallocation(m_Pmat,0,d_nnz,0,o_nnz);CHKERRQ(ierr);
    }
  if (m_num_extra_nnz) MatSetOption(m_mat,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE);
  
  // assemble
  for (ilev=0,idx=0;ilev<nGrids;ilev++)
    {
      const DisjointBoxLayout& dbl = m_grids[ilev];
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
          BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          Box region = dbl[dit];
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if ( gidfab(iv,0) >= 0)
                {
                  ierr = AddStencilToMat(iv,ilev,dit(),stenVect[idx],m_mat);CHKERRQ(ierr);
                  idx++;
                }
            }
        }
    }
  CH_assert(idx==nloc);
  PetscFree(d_nnz);  PetscFree(o_nnz);

  // total assesmbly
  ierr = MatAssemblyBegin(m_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);  
  if (a_makePmat)
    {
      ierr = MatAssemblyBegin(m_Pmat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(m_Pmat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);  
    }

// finish off repartitioning
if (m_repartition) 
    {
#if PETSC_VERSION_LT(3,5,4)
      ierr = MatGetVecs(m_mat,0,&m_origvec);CHKERRQ(ierr); // keep old vec around
#else
      ierr = MatCreateVecs(m_mat,0,&m_origvec);CHKERRQ(ierr); // keep old vec around
#endif

#if defined(PETSC_USE_LOG)
      PetscLogEventBegin(m_event2,0,0,0,0);
#endif
      // partition, move matrix, change meta data
      ierr = permuteDataAndMaps(patchStencil);CHKERRQ(ierr);

#if defined(PETSC_USE_LOG)
      PetscLogEventEnd(m_event2,0,0,0,0);
#endif
    }
  
  if (m_writeMatlab)
    {
      PetscViewer viewer;
      char suffix[30];
#ifdef CH_MPI
      MPI_Comm comm = Chombo_MPI::comm;
#else
      MPI_Comm comm = PETSC_COMM_SELF;
#endif
      sprintf(suffix, "A%d.m",nGrids);
      ierr = PetscViewerASCIIOpen(comm, suffix, &viewer);CHKERRQ(ierr);
      ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
      ierr = MatView(m_mat,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);
    }
#if defined(PETSC_USE_LOG)
  PetscLogEventEnd(m_event0,0,0,0,0);
#endif
  PetscFunctionReturn(0);
}
//
//  Helper method: partition graph (input), covert matrix, create map objects
//
#undef __FUNCT__
#define __FUNCT__ "permuteDataAndMaps"
PetscErrorCode 
PetscCompGrid::permuteDataAndMaps(Vector<StencilTensor> &patchStencil)
{
  PetscErrorCode ierr;
  Mat graph, graph_new, new_mat;
  IS old_new_patchids,new_old_patchids;
  PetscInt npatch_new,nloc_new,nnewlocpatches,*tidx,nblockpts=pow((Real)m_patch_size,SpaceDim),old_nloceq;
  PetscMPIInt nprocs;
  const PetscInt *idx;
  PetscFunctionBeginUser;
  CH_assert(patchStencil.size()==m_nlocrealpatches);
#ifdef CH_MPI
  ierr = MPI_Comm_size(Chombo_MPI::comm,&nprocs);CHKERRQ(ierr);
  CH_assert(nprocs>1);
  ierr = VecGetLocalSize(m_origvec,&old_nloceq);CHKERRQ(ierr);

  // make symmetric graph
  { 
    PetscScalar v=1.0;
    ierr = MatCreate(Chombo_MPI::comm,&graph);CHKERRQ(ierr);
    ierr = MatSetSizes(graph,m_nlocrealpatches,m_nlocrealpatches,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
    MatSetType(graph,MATAIJ);
    ierr = MatMPIAIJSetPreallocation(graph,250,0,150,0);CHKERRQ(ierr);
    for (PetscInt ii=0,gidi=m_patchid0;ii<m_nlocrealpatches;ii++,gidi++)
      {
        StencilTensor &sten = patchStencil[ii]; // add to this row
        StencilTensor::const_iterator end = sten.end(); 
        for (StencilTensor::const_iterator it = sten.begin(); it != end; ++it) 
          {
            //const IndexML &ivJ = it->first;
            PetscInt gidj = (PetscInt)it->second.value(0,0) - 1; // zero based
            ierr = MatSetValues(graph,1,&gidi,1,&gidj,&v,INSERT_VALUES);CHKERRQ(ierr);
          }
      }
    ierr = MatAssemblyBegin(graph,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(graph,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 

    // symetrize
    ierr = MatTranspose(graph, MAT_INITIAL_MATRIX, &graph_new);CHKERRQ(ierr);
    ierr = MatAXPY(graph,1.0,graph_new,DIFFERENT_NONZERO_PATTERN);
    MatDestroy(&graph_new);
  }
  if (false) {    
    PetscViewer viewer;
    char suffix[30] = "G0.m";
    ierr = PetscViewerASCIIOpen(Chombo_MPI::comm,suffix,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = MatView(graph,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  // partition -- old_new_patchids, new_old_patchids, new sizes
  {
    MatPartitioning mpart;
    IS              proc_is;
    PetscMPIInt     rank;
    PetscInt       *counts;

    ierr = MPI_Comm_rank(Chombo_MPI::comm, &rank);CHKERRQ(ierr);        
    ierr = MatPartitioningCreate(Chombo_MPI::comm, &mpart);CHKERRQ(ierr);
    ierr = MatPartitioningSetAdjacency(mpart,graph);CHKERRQ(ierr);      
    ierr = MatPartitioningSetFromOptions(mpart);CHKERRQ(ierr);
    ierr = MatPartitioningSetNParts(mpart,nprocs);CHKERRQ(ierr);
    ierr = MatPartitioningApply(mpart, &proc_is);CHKERRQ(ierr);
    ierr = MatPartitioningDestroy(&mpart);CHKERRQ(ierr);
    ierr = ISPartitioningToNumbering(proc_is, &new_old_patchids);CHKERRQ(ierr);
    // Determine how many equations/vertices are assigned to each processor -- new sizes
    PetscMalloc1(nprocs, &counts);
    ierr = ISPartitioningCount(proc_is, nprocs, counts);CHKERRQ(ierr);
    ierr = ISDestroy(&proc_is);CHKERRQ(ierr);
    npatch_new = counts[rank];
    nnewlocpatches = counts[rank];
    nloc_new = nnewlocpatches*nblockpts;
    PetscFree(counts);
    // Invert for MatGetSubMatrix -- expensive
    ierr = ISInvertPermutation(new_old_patchids, npatch_new, &old_new_patchids);CHKERRQ(ierr);
    
    ierr = ISSort(old_new_patchids);CHKERRQ(ierr); /* is this needed? */
    // MatGetSubMatrix -- not needed
    ierr = MatGetSubMatrix(graph, old_new_patchids, old_new_patchids, MAT_INITIAL_MATRIX, &graph_new);CHKERRQ(ierr);
    if (false) {          
      PetscViewer viewer;
      char suffix[30] = "G2.m";
      ierr = PetscViewerASCIIOpen(Chombo_MPI::comm,suffix,&viewer);CHKERRQ(ierr);
      ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
      ierr = MatView(graph_new,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }
  }
  // get new matrix, get maps to CH
  {
    IS old_new_geqs;
    PetscInt ii,jj;

    // old_new_geqs: the old GIDS of my new locals
    ierr = PetscMalloc1(nloc_new*m_dof, &tidx);CHKERRQ(ierr);
    ierr = ISGetIndices(old_new_patchids, &idx);CHKERRQ(ierr);
    for (ii=0,jj=0; ii<npatch_new; ii++) 
      {
        PetscInt geq = idx[ii]*nblockpts*m_dof; // equation
        for (int kk=0; kk<nblockpts*m_dof; kk++, jj++, geq++) tidx[jj] = geq;
      }
    CH_assert(jj==nloc_new*m_dof);
    ierr = ISRestoreIndices(old_new_patchids, &idx);CHKERRQ(ierr);
    ierr = ISDestroy(&old_new_patchids);CHKERRQ(ierr);
    ierr = ISCreateGeneral(Chombo_MPI::comm,nloc_new*m_dof,tidx,PETSC_COPY_VALUES,&old_new_geqs);CHKERRQ(ierr);
    ierr = PetscFree(tidx);CHKERRQ(ierr);
    // get new matrix
    ierr = MatGetSubMatrix(m_mat, old_new_geqs, old_new_geqs, MAT_INITIAL_MATRIX, &new_mat);CHKERRQ(ierr);

    Vec newvec;
    ierr = VecCreate(Chombo_MPI::comm,&newvec);CHKERRQ(ierr);
    ierr = VecSetType(newvec,VECSTANDARD);CHKERRQ(ierr);
    ierr = VecSetSizes(newvec,nloc_new*m_dof,PETSC_DECIDE);CHKERRQ(ierr);

    // make map from PETSc to Chombo
    ierr = VecScatterCreate(newvec,NULL,m_origvec,old_new_geqs,&m_from_petscscat);CHKERRQ(ierr);
    //
    ierr = VecDestroy(&newvec);CHKERRQ(ierr);
    ierr = ISDestroy(&old_new_geqs);CHKERRQ(ierr);
    ierr = ISDestroy(&new_old_patchids);
  }
  // set new maps, cleanup
  m_nlocrealpatches = 0; // could use args for this
  m_patchid0 = 0;        // could use args for this
  ierr = MatDestroy(&graph_new);CHKERRQ(ierr);
  ierr = MatDestroy(&graph);CHKERRQ(ierr);    
  ierr = MatDestroy(&m_mat);CHKERRQ(ierr);
  m_mat = new_mat;
  if (false) {
    PetscViewer viewer;
    char suffix[30] = "A2.m";
    ierr = PetscViewerASCIIOpen(Chombo_MPI::comm,suffix,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = MatView(m_mat,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
#else
  CH_assert(0);
#endif
  PetscFunctionReturn(0);
}

static bool s_check_row_sum = false;

void 
PetscCompGrid::applyBCs( IntVect a_iv, int a_ilev, const DataIndex &a_dummy, Box a_dombox, 
                         StencilTensor &a_sten )
{
  CH_TIME("PetscCompGrid::applyBCs");
  int nSources = -1;
  Vector<Vector<StencilNode> > new_vals;

  // count degrees, add to 'corners'
  Vector<IntVectSet> corners(CH_SPACEDIM);
  StencilTensor::const_iterator end2 = a_sten.end(); 
  for (StencilTensor::const_iterator it = a_sten.begin(); it != end2; ++it) 
    {
      const IntVect &jiv = it->first.iv();
      if (!a_dombox.contains(jiv) && it->first.level()==a_ilev) 
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
      for (IVSIterator ivit(corners[ideg]); ivit.ok(); ++ivit)
        {
          const IntVect &jiv = ivit(); // get rid of this bc ghost
          // get layer of ghost
          Box gbox(IntVect::Zero,IntVect::Zero); gbox.shift(jiv);
          Box inter = a_dombox & gbox;    CH_assert(inter.numPts()==0);
          PetscInt igid = -1; // find which layer of ghost I am
          do{
            igid++;
            gbox.grow(1);
            inter = a_dombox & gbox;
          }while (inter.numPts()==0);
          for (int dir=0;dir<CH_SPACEDIM;++dir)
            {
              if (jiv[dir] < a_dombox.smallEnd(dir) || jiv[dir] > a_dombox.bigEnd(dir)) 
                {
                  // have a BC, get coefs on demand
                  if (nSources == -1) { 
                    IntVect ghostVect = getGhostVect();
                    RefCountedPtr<BCFunction> bc = m_bc.getBCFunction();
                    ConstDiriBC *mybc = dynamic_cast<ConstDiriBC*>(&(*bc));
                    if (mybc)
                      {
                        CH_assert(mybc->nGhosts()==ghostVect);
                        nSources = mybc->nSources();
                        new_vals.resize(ghostVect[0]);
                        for (int i = 0 ; i < ghostVect[0]; i++)
                          {
                            new_vals[i].resize(nSources); 
                            for (int j = 0 ; j < nSources; j++ )
                              {
                                new_vals[i][j].second.define(1); // have to define because it can be a tensor now, yuck
                                new_vals[i][j].second.setValue(mybc->getCoef(j,i));
                                new_vals[i][j].first.setLevel(a_ilev);
                              }
                          }
                      }
                    else
                      {
                        CH_assert(ghostVect[0]==1); 
                        Real denom=2.;   new_vals.resize(1); new_vals[0].resize(2);
                        new_vals[0][0].second.define(1);           new_vals[0][1].second.define(1); // have to define because it can be a tensor now, yuck
                        new_vals[0][0].second.setValue(-5./denom); new_vals[0][1].second.setValue(1./denom);
                        new_vals[0][0].first.setLevel(a_ilev);     new_vals[0][1].first.setLevel(a_ilev);
                      }             
                    // for debugging -- neumann
                    if (s_check_row_sum)
                      {
                        new_vals.resize(ghostVect[0]);
                        for (int idx=0;idx<ghostVect[0];idx++)
                          {
                            new_vals[idx].resize(1);
                            new_vals[idx][0].second.setValue(1.);
                          }
                        nSources = 1;
                      }
                  } // init coefs
                  
                  IndexML kill(jiv,a_ilev);
                  int isign = 1;
                  if (jiv[dir] < a_dombox.smallEnd(dir)) isign = -1;
                  IntVect biv = jiv;
                  biv.shift(dir,-isign*(igid+1));                 
                  for (int j=0; j<nSources; j++, biv.shift(dir,-isign))
                    {
                      new_vals[igid][j].first.setIV(biv);
                      if (ideg>0 && !a_dombox.contains(biv)) // this could be a new stencil value for high order BCs
                        {
                          corners[ideg-1] |= biv;
                        }
                      else CH_assert(a_dombox.contains(biv));
                    }
                  StencilProject(kill,new_vals[igid], a_sten);
                  int nrm = a_sten.erase(kill);CH_assert(nrm==1);
                  break;
                }
            }
        } // ghosts
    } // degree
  
  // debug
  if (0)
    {
      double summ=0.;
      StencilTensor::const_iterator end = a_sten.end(); 
      for (StencilTensor::const_iterator it = a_sten.begin(); it != end; ++it) {
        for (int i=0; i<m_dof; ++i) 
          {
            for (int j=0; j<m_dof; ++j) 
              {
                summ += it->second.value(i,j);
              }
          }
      }
      if (Abs(summ)>1.e-10)
        {
          for (StencilTensor::const_iterator it = a_sten.begin(); it != end; ++it) {
            pout() << it->first << ": \n";
            for (int i=0; i<m_dof; ++i) 
              {
                pout() << "\t\t";
                for (int j=0; j<m_dof; ++j) 
                  {
                    pout() << it->second.value(i,j) << ", ";
                  }
                pout() << endl;
              }
          }
          pout() << "\t ERROR summ (if Neumann): " << summ << endl; 
        }
    }
}
//
// BCs have been removed (should not care).
//
void 
PetscCompGrid::InterpToCoarse(IntVect a_iv,int a_ilev,const DataIndex &a_di, StencilTensor &a_sten)
{
  CH_TIME("PetscCompGrid::InterpToCoarse");
  int nref=m_refRatios[a_ilev-1];
  int nfine = (int)pow((Real)nref,SpaceDim);
  BaseFab<PetscInt>& supgidfab = (*m_crsSupportGIDs[a_ilev-1])[a_di];
  BaseFab<PetscInt>& gidfab = (*m_GIDs[a_ilev])[a_di];
  ProblemDomain cdom = m_domains[a_ilev-1];
  std::list<StencilTensor::iterator> killList;

  // remove ghosts on domain
  StencilTensor::iterator end = a_sten.end(); 
  for (StencilTensor::iterator it = a_sten.begin(); it != end; ++it) 
    {
      const IndexML &mlivJ = it->first;
      if (mlivJ.level()==a_ilev && gidfab(mlivJ.iv(),0) == GHOST)
        {
          const IntVect &jiv = mlivJ.iv(); // kill
          IntVect cjiv = jiv; cjiv.coarsen(nref);
          const IntVect fivBase = cjiv*nref;
          const IntVect fividx = jiv - fivBase;
          IntVect bndryIndex = getCFStencil(cdom,cjiv);
          const FourthOrderInterpStencil *interp = m_FCStencils(bndryIndex,0);
          const FArrayBox &coarseToFineFab = interp->m_coarseToFineFab;
          Vector<StencilNode> new_vals(interp->m_stencilSize);
          // OLD WAY: FourthOrderInterpStencil::m_coarseBaseIndices
          // WAS Vector<int>, is now Vector<IntVect>.
          //          const Vector<int> &coarseBaseIndices = interp->m_coarseBaseIndices;
          //          Vector<StencilNode> new_vals(interp->m_stencilSize);
          //          for (PetscInt cidx=0,kk=0;cidx < interp->m_stencilSize;cidx++)
          //            {
          //              IntVect civ = cjiv;
          //              for (int dir=0;dir<CH_SPACEDIM;++dir,kk++)
          //                {
          //                  civ[dir] += coarseBaseIndices[kk];
          //                }
          //              Real val = coarseToFineFab(fividx,cidx);
          //              NodeDefine(new_vals[cidx],civ,a_ilev-1,val);
          //            }
          const Vector<IntVect> &coarseBaseIndices = interp->m_coarseBaseIndices;
          for (PetscInt cidx=0;cidx < interp->m_stencilSize;cidx++)
            {
              IntVect civ = cjiv + coarseBaseIndices[cidx];
              Real val = coarseToFineFab(fividx,cidx);
              NodeDefine(new_vals[cidx],civ,a_ilev-1,val);
            }
          StencilProject(mlivJ,new_vals,a_sten);
          killList.push_back (it);
        }
    }

  // get rid of coarse covered with simple averaging
  end = a_sten.end(); 
  for (StencilTensor::iterator it = a_sten.begin(); it != end; ++it) 
    {
      const IndexML &mlivJ = it->first;
      IntVect civ = mlivJ.iv();
      if (mlivJ.level()==a_ilev-1 && supgidfab(civ,0)<0)
        {
          if (supgidfab(civ,0)!=FINE_COVERED) {
            pout() << "    ERROR row " << IndexML(a_iv,a_ilev) << ", column " << mlivJ << ", type = " << supgidfab(civ,0) << endl;
            pout() << "   Use nesting radius >= 3!!!" << endl;
            CH_assert(supgidfab(civ,0) == FINE_COVERED);
          }
          Real fineInterp = 1./(Real)nfine;
          const IntVect fivBase = civ*nref; // fine box low end to interp to
          Box fbox(IntVect::Zero,(nref-1)*IntVect::Unit); fbox.shift(fivBase);
          Vector<StencilNode> new_vals(nfine);
          BoxIterator bit(fbox);
          for (int i=0; bit.ok(); ++bit,i++)
            {
              NodeDefine(new_vals[i],bit(),a_ilev,fineInterp);
            }
          StencilProject(mlivJ,new_vals,a_sten);
          killList.push_back (it);
        }
    }
  // erase killed
  for (std::list<StencilTensor::iterator>::iterator lit=killList.begin(); 
       lit != killList.end(); 
       ++lit)
    {
      a_sten.erase(*lit);
    }
}

// F: talk to fine grid covers -- refluxing
void 
PetscCompGrid::InterpToFine(IntVect a_iv, int a_ilev,const DataIndex &a_di, StencilTensor &a_sten )
{
  CH_TIME("PetscCompGrid::InterpToFine");
  int nref=m_refRatios[a_ilev];
  BaseFab<PetscInt>& gidfab = (*m_GIDs[a_ilev])[a_di];
  BaseFab<PetscInt>& covergidfab = (*m_fineCoverGIDs[a_ilev+1])[a_di];
  Box coverfabbox = covergidfab.box();
  std::list<StencilTensor::iterator> killList;
  
  // get rid of coarse covered
  StencilTensor::iterator end = a_sten.end(); 
  for (StencilTensor::iterator it = a_sten.begin(); it != end; ++it) 
    {
      const IndexML &mlivJ = it->first;
      const IntVect &jiv = mlivJ.iv();
      if (mlivJ.level()==a_ilev && gidfab(jiv,0)<0)     
        {
          // simple averaging
          const IntVect fivBase = jiv*nref; // fine box low end to interp to
          Box fbox(IntVect::Zero,(nref-1)*IntVect::Unit); fbox.shift(fivBase);
          int nfine = (int)pow((Real)nref,SpaceDim);
          Vector<StencilNode> new_vals(nfine);
          Real fineInterp = 1./(Real)nfine;
          BoxIterator bit(fbox);
          for (int i=0; bit.ok(); ++bit,i++)
            {
              NodeDefine(new_vals[i],bit(),a_ilev+1,fineInterp);
            }
          StencilProject(mlivJ,new_vals,a_sten);
          killList.push_back (it);
        }
    }
  // erase killed - would like to have StencilProject do this
  for (std::list<StencilTensor::iterator>::iterator lit=killList.begin(); 
       lit != killList.end(); 
       ++lit)
    {
      a_sten.erase(*lit);
    }
}

// add row of matrix (stencil)
#undef __FUNCT__
#define __FUNCT__ "AddStencilToMat"
PetscErrorCode 
PetscCompGrid::AddStencilToMat(IntVect a_iv, int a_ilev,const DataIndex &a_di, StencilTensor &a_sten, Mat a_mat)
{
  PetscErrorCode ierr;
  PetscScalar vals[4096],zerovals[4096];
  PetscInt cols[4096],iLevel=a_ilev,vidx[STENCIL_MAX_DOF],ncols=m_dof*a_sten.size();
  BaseFab<PetscInt>& this_gidfabJ = (*m_GIDs[iLevel])[a_di];
  double summ=0.,abssum=0.;
  PetscFunctionBeginUser;

  if (a_sten.size()*m_dof > 4096) MayDay::Error("PetscCompGrid::AddStencilToMat buffer (4096) too small");
  // get cols & vals
  for (PetscInt ni=0,geq = this_gidfabJ(a_iv,0)*m_dof;ni<m_dof;ni++,geq++) vidx[ni] = geq;

  StencilTensor::const_iterator end = a_sten.end(); 
  int ci=0,jj=0;
  for (StencilTensor::const_iterator it = a_sten.begin(); it != end; ++it,jj++) 
    {
      const IndexML &ivJ = it->first;
      int jLevel = ivJ.level();
      BaseFab<PetscInt>& gidfabJ = (jLevel==iLevel) ? this_gidfabJ : 
        (jLevel==iLevel+1) ? (*m_fineCoverGIDs[iLevel+1])[a_di] : 
        (*m_crsSupportGIDs[iLevel-1])[a_di];
                  
      if (!gidfabJ.box().contains(ivJ.iv())){
        pout() << "ERROR row  " << IndexML(a_iv,a_ilev) << ", this box:" << 
          m_grids[a_ilev][a_di] << ", fine box (failure): " << gidfabJ.box() << 
          " failed to find "<< ivJ << endl;
        jj=0;
        for (StencilTensor::const_iterator it2 = a_sten.begin(); it2 != end; ++it2) pout()<<++jj<<") j="<<it2->first<<endl;
        MayDay::Error("PetscCompGrid::AddStencilToMat failed to find cell");
      }

      PetscInt gidj = gidfabJ(ivJ.iv(),0)*m_dof;
      if (gidj<0) 
        {
          pout() << "\tFAILED TO FIND iv=" << ivJ << "\t gidj type:" << (GID_type)gidj << "\tiLevel=" << iLevel << "\tiv=" << a_iv << endl;
          MayDay::Error("PetscCompGrid::AddStencilToMat failed to find gid");
        }
      const Real *vv = it->second.getVals();
      if (*vv!=*vv) MayDay::Error("PetscCompGrid::AddStencilToMatit->second.getVals(0,0) is a NaN");
      for (int nj=0;nj<m_dof;nj++,gidj++,ci++) cols[ci] = gidj;   // columns
      for (int ni=0;ni<m_dof;ni++) 
        {
          for (int nj=0;nj<m_dof;nj++) 
            {
              double tt = vv[ni*m_dof + nj];
              if (tt!=tt) MayDay::Error("PetscCompGrid::AddStencilToMat found a NaN");
              if (abs(tt)>1.e300)  MayDay::Error("PetscCompGrid::AddStencilToMat found a BigNum");
              vals[ni*ncols + jj*m_dof + nj] = tt;
              summ += tt;
              abssum += Abs(tt);
              if (m_num_extra_nnz) zerovals[ni*ncols + jj*m_dof + nj] = 0.;
            }
        }
    }

  ierr = MatSetValues(a_mat,m_dof,vidx,ncols,cols,vals,ADD_VALUES);CHKERRQ(ierr);
  if (m_num_extra_nnz) 
    {
      ierr = MatSetValues(a_mat,ncols,cols,m_dof,vidx,zerovals,ADD_VALUES);CHKERRQ(ierr);
    }

  // debug
  if (s_check_row_sum && Abs(summ)/abssum>1.e-5)
    {
      pout() << " error gid=" << vidx[0]+1 << " row iv =" << a_iv << endl;
      for (int idx=0;idx<ncols;idx++)
        {
          if (vals[idx]!=0.0)
            {
              PetscInt gidj = cols[idx];
              pout() << "\t" << gidj+1 << ", v = " << vals[idx] << endl;
            }
        }
      SETERRQ4(PETSC_COMM_WORLD,1,"sum=%e on level %d, iv=[%d,%d]",summ,iLevel,a_iv[0],a_iv[1]);
    }

  PetscFunctionReturn(0);
}

//
IntVect PetscCompGrid::getCFStencil(const ProblemDomain &a_cdom, const IntVect a_ivc)
{
  const Box& coarseDomainBox = a_cdom.domainBox();
  const IntVect& coarseDomainLo = coarseDomainBox.smallEnd();
  const IntVect& coarseDomainHi = coarseDomainBox.bigEnd();
  IntVect dist = IntVect::Zero;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (!a_cdom.isPeriodic(idir))
        {
          int offLo = coarseDomainLo[idir] - a_ivc[idir] - 1;
          int offHi = coarseDomainHi[idir] - a_ivc[idir] + 1;
          if (offLo < 0 && offHi >0) // condition means ivc in coarseDomain
            {
              if ((offLo >= -m_CFStencilRad) &&
                  (offHi <= m_CFStencilRad))
                { // both of these:  very narrow domain, you are in trouble
                  MayDay::Error("getCFStencil::define bad boxes");
                }
              if (offLo >= -m_CFStencilRad) // -1 or -2
                dist[idir] = offLo;
              if (offHi <= m_CFStencilRad) // 1 or 2
                dist[idir] = offHi;
              // Otherwise, dist[idir] = 0.
            }
        }
    }
  return dist;
}
#undef __FUNCT__
#define __FUNCT__ "putChomboInPetsc"
PetscErrorCode
PetscCompGrid::putChomboInPetsc(const Vector<LevelData<FArrayBox> * > &a_rhs, Vec a_b)const
{
  CH_TIME("PetscCompGrid::putChomboInPetsc");
  PetscInt gid,nGrids=a_rhs.size(),*idxj,idx;
  PetscErrorCode ierr;
  PetscScalar *vals;
  Vec bb;
  PetscFunctionBeginUser;

  if (m_repartition) bb = m_origvec;
  else bb = a_b;

  // iterate over real cells, get gids, make ops
  for (int ilev=nGrids-1;ilev>=0;ilev--)
    {
      const int nc=a_rhs[ilev]->nComp();
      const DisjointBoxLayout &dbl = a_rhs[ilev]->getBoxes();
      for (DataIterator dit = a_rhs[ilev]->dataIterator(); dit.ok(); ++dit)
        {
          Box region = dbl[dit];
          const BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          const FArrayBox& phiFAB = (*a_rhs[ilev])[dit];
          idx = nc*region.numPts();
          PetscMalloc1(idx, &idxj);
          PetscMalloc1(idx, &vals);
          idx = 0;
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if ( (gid=gidfab(iv,0)) >= 0)
                {
                  for (PetscInt i=nc*gid,n=0;n<nc;n++,i++)
                    {
                      idxj[idx] = i;
                      vals[idx] = phiFAB(iv,n);
                      idx++;
                    }
                }
            }
          ierr = VecSetValues(bb,idx,idxj,vals,INSERT_VALUES);CHKERRQ(ierr);
          PetscFree(vals);
          PetscFree(idxj);
        }
    }
  ierr = VecAssemblyBegin(bb);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(bb);CHKERRQ(ierr);
  
  if (m_repartition) 
    {
      ierr = VecScatterBegin(m_from_petscscat,m_origvec,a_b,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
      ierr = VecScatterEnd  (m_from_petscscat,m_origvec,a_b,INSERT_VALUES,SCATTER_REVERSE);CHKERRQ(ierr);
    }

  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "putPetscInChombo"
PetscErrorCode 
PetscCompGrid::putPetscInChombo(Vec a_x, Vector<LevelData<FArrayBox> * > &a_phi)const
{
  CH_TIME("PetscCompGrid::putPetscInChombo");
  PetscInt gid,nGrids=a_phi.size();
  PetscErrorCode ierr;
  const PetscInt nc=a_phi[0]->nComp(), my0eq=nc*m_gid0; CH_assert(nc==m_dof);
  const PetscScalar *avec;
  Vec xx;
  PetscFunctionBeginUser;

  if (m_repartition) 
    {
      ierr = VecScatterBegin(m_from_petscscat,a_x,m_origvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
      ierr = VecScatterEnd  (m_from_petscscat,a_x,m_origvec,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
      xx = m_origvec;
    }
  else xx = a_x;

  ierr = VecGetArrayRead(xx,&avec);CHKERRQ(ierr);
  // iterate over real cells, get gids, make ops
  for (int ilev=nGrids-1;ilev>=0;ilev--)
    {
      const DisjointBoxLayout &dbl = a_phi[ilev]->getBoxes();
      for (DataIterator dit = a_phi[ilev]->dataIterator(); dit.ok(); ++dit)
        {
          Box region = dbl[dit];
          const BaseFab<PetscInt>& gidfab = (*m_GIDs[ilev])[dit];
          FArrayBox& phiFAB = (*a_phi[ilev])[dit];
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if ( (gid=gidfab(iv,0)) >= 0)
                {
                  for (PetscInt i=nc*gid-my0eq,n=0;n<nc;n++,i++) phiFAB(iv,n) = avec[i];
                }
              else
                {
                  for (int n=0;n<nc;n++) phiFAB(iv,n) = 0.;
                }
            }
        }
      //a_phi[ilev]->exchange();
    }
  // average fine level solution onto coarse levels
  if (m_averageFineSolutionToCoarse)
    {
      for (int lev = a_phi.size() -1; lev > 0; lev--)
        {
          CoarseAverage avg(a_phi[lev]->disjointBoxLayout(),a_phi[lev-1]->disjointBoxLayout(),
                            a_phi[lev]->nComp(),m_refRatios[lev-1]);
          avg.averageToCoarse(*a_phi[lev-1], *a_phi[lev]);
        }
    }
  
  ierr = VecRestoreArrayRead(xx,&avec);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#endif

#include "NamespaceFooter.H"


