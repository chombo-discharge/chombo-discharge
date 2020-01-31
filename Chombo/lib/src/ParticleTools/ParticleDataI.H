#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

 //  ANAG, LBNL

#ifndef _PARTICLEDATAI_H_
#define _PARTICLEDATAI_H_

#include "MPI_util.H"
#include "NamespaceHeader.H"

using std::vector;
using std::map;
using std::binary_function;

//
// Implementation
//

template <class P>
void collectValidParticles(List<P>&                         a_partValid,
                           ParticleData<P>&                 a_PD,
                           const LevelData<BaseFab<bool> >* a_mask,
                           const RealVect                   a_meshSpacing,
                           const int                        a_refRatio,
                           const bool                       a_flip,
                           const RealVect                   a_origin)
{
  CH_TIMERS("collectValidParticles(List)");

  BoxLayout fromBL;
  if   (a_refRatio>0) refine (fromBL,a_PD.getBoxes(), a_refRatio);
  else                coarsen(fromBL,a_PD.getBoxes(),-a_refRatio);

  BoxLayoutData<BaseFab<bool> > valid(fromBL,1);
  for (DataIterator di(fromBL); di.ok(); ++di)
    {
      valid[di].setVal(false);
    }
  a_mask->copyTo(valid);

  for (DataIterator di(fromBL); di.ok(); ++di)
    {
      for (ListIterator<P> li (a_PD[di].listItems()); li.ok(); )
        {
          const RealVect x=li().position();
          const IntVect iv(D_DECL6((int)floor((x[0] - a_origin[0])/a_meshSpacing[0]),
				   (int)floor((x[1] - a_origin[1])/a_meshSpacing[1]),
				   (int)floor((x[2] - a_origin[2])/a_meshSpacing[2]),
				   (int)floor((x[3] - a_origin[3])/a_meshSpacing[3]),
				   (int)floor((x[4] - a_origin[4])/a_meshSpacing[4]),
				   (int)floor((x[5] - a_origin[5])/a_meshSpacing[5])));

          if (valid[di](iv) ^ a_flip) a_partValid.transfer(li);
          else                      ++li;
        }
    }
}

template <class P>
void collectValidParticles(ParticleData<P>&                 a_PDValid,
                           ParticleData<P>&                 a_PD,
                           const LevelData<BaseFab<bool> >* a_mask,
                           const RealVect                   a_meshSpacing,
                           const int                        a_refRatio,
                           const bool                       a_flip,
                           const RealVect                   a_origin)
{
  CH_TIMERS("collectValidParticles(ParticleData)");

  BoxLayout fromBL;
  if   (a_refRatio>0) refine (fromBL,a_PD.getBoxes(), a_refRatio);
  else                coarsen(fromBL,a_PD.getBoxes(),-a_refRatio);

  BoxLayoutData<BaseFab<bool> > valid(fromBL,1);
  for (DataIterator di(fromBL); di.ok(); ++di)
    {
      valid[di].setVal(false);
    }
  a_mask->copyTo(valid);

  for (DataIterator di(fromBL); di.ok(); ++di)
    {
      for (ListIterator<P> li (a_PD[di].listItems()); li.ok(); )
        {
          const RealVect x=li().position();
          const IntVect iv(D_DECL6((int)floor((x[0] - a_origin[0])/a_meshSpacing[0]),
				   (int)floor((x[1] - a_origin[1])/a_meshSpacing[1]),
				   (int)floor((x[2] - a_origin[2])/a_meshSpacing[2]),
				   (int)floor((x[3] - a_origin[3])/a_meshSpacing[3]),
				   (int)floor((x[4] - a_origin[4])/a_meshSpacing[4]),
				   (int)floor((x[5] - a_origin[5])/a_meshSpacing[5])));

          if (valid[di](iv) ^ a_flip)
            a_PDValid[di].listItems().transfer(li);
          else
            ++li;
        }
    }
}

inline void myswap(pair<const IntVect,boxids>& aa, pair<const IntVect,boxids>& bb)
{
  pair<IntVect,boxids> a=(pair<IntVect,boxids>)aa;
  pair<IntVect,boxids> b=(pair<IntVect,boxids>)bb;
  pair<IntVect,boxids> t=make_pair(b.first,b.second);
  a.first=b.first;
  a.second=b.second;

  b.first=t.first;
  b.second=t.second;
}

template <class ForwardIterator>
  void myrotate ( ForwardIterator first, ForwardIterator middle,
		  ForwardIterator last )
{
  ForwardIterator next = middle;
  while (first!=next)
  {
    myswap (*first++,*next++);
    if (next==last) next=middle;
    else if (first == middle) middle=next;
  }
}

typedef multimap<IntVect, boxids>::iterator RmIt;

template <class P>
class CompCounts : public binary_function<pair<IntVect,boxids>,
                                          pair<IntVect,boxids>,bool>
{
public:
  // constructor
  explicit CompCounts(vector<map<unsigned,List<P> > >& a_pp)
    : pp(a_pp) {}

  bool operator() (const pair<const IntVect,boxids>& a,
                   const pair<const IntVect,boxids>& b) const
  {
    return (pp[a.second.pid][a.second.idx].length() <
            pp[b.second.pid][b.second.idx].length());
  }

protected:
  vector<map<unsigned, List<P> > >& pp;
};

template <class P>
ParticleData<P>::ParticleData() 
  : 
  LayoutData<ListBox<P> >() 
{
  m_isDefined = false;
}

template <class P>
ParticleData<P>::ParticleData( const BoxLayout&         a_dp,
			       const ProblemDomain&     a_domain,
			       const int&               a_fixedBoxSize,
                               const RealVect&          a_meshSpacing,
			       const RealVect&          a_origin)
{
  define(a_dp, a_domain, a_fixedBoxSize, a_meshSpacing, a_origin);
}

template <class P>
void ParticleData<P>::define( const BoxLayout&         a_dp,
			      const ProblemDomain&     a_domain,
			      const int&               a_fixedBoxSize,
			      const RealVect&          a_meshSpacing,
			      const RealVect&          a_origin)
{
  LayoutData<ListBox<P> >::define(a_dp);
  m_physDomain = a_domain;
  m_fixedBoxSize = a_fixedBoxSize;
  m_meshSpacing = a_meshSpacing;
  m_origin = a_origin;
  m_factory.define(m_meshSpacing, m_origin);
  allocateVector();
  m_isDefined = true;
}


template <class P>
const BoxLayout& ParticleData<P>::getBoxes() const
{
  return LayoutData<ListBox<P> >::boxLayout();
}

template <class P>
void ParticleData<P>::clear()
{
  if (m_isDefined)
    {
      for (DataIterator di=this->dataIterator(); di.ok(); ++di)
        {
          (*this)[di].listItems().clear();
        }
      m_outcast.clear();
    }
}

template <class P>
size_t ParticleData<P>::numParticlesLocal() const
{
  if (!m_isDefined) return 0;
  else
    {
      size_t numParticles = 0;
      for (DataIterator di=this->dataIterator(); di.ok(); ++di)
        {
          numParticles += (*this)[di].numItems();
        }
      numParticles +=m_outcast.length();

      return numParticles;
    }
}

template <class P>
size_t ParticleData<P>::numParticles() const
{
  if (!m_isDefined) return 0;
  else
    {
      size_t countLocal;
      countLocal = numParticlesLocal();
      
      // gather all countLocals onto process 0
      int srcProc = 0;
      Vector<size_t> allCounts(numProc());
      gather(allCounts, countLocal, srcProc);
      
      size_t count = 0;
      if(procID() == srcProc)
	{
	  for(int ivec = 0; ivec < numProc(); ivec++)
	    {
	      count += allCounts[ivec];
	    }
	}
      
      // broadcast right answer to all procs
      broadcast(count, srcProc);
      return count;
    }
}

template <class P>
size_t ParticleData<P>::numValidLocal() const
{
  if (!m_isDefined) return 0;
  else
    {
      size_t numValid = 0;
      for (DataIterator di=this->dataIterator(); di.ok(); ++di)
        {
          numValid += (*this)[di].numItems();
        }
      return numValid;
    }
}

template <class P>
size_t ParticleData<P>::numValid() const
{
  if (!m_isDefined) return 0;
  else
    {
      size_t countLocal;
      countLocal = numValidLocal();
      
      // gather all countLocals onto process 0
      int srcProc = 0;
      Vector<size_t> allCounts(numProc());
      gather(allCounts, countLocal, srcProc);
      
      size_t count = 0;
      if(procID() == srcProc)
	{
	  for(int ivec = 0; ivec < numProc(); ivec++)
	    {
	      count += allCounts[ivec];
	    }
	}
      
      // broadcast right answer to all procs
      broadcast(count, srcProc);
      return count;
    }
}

template <class P>
size_t ParticleData<P>::numOutcastLocal() const
{
  if (!m_isDefined) return 0;
  else
    {
      size_t numOutcast = m_outcast.length();
      return numOutcast;
    }
}

template <class P>
size_t ParticleData<P>::numOutcast() const
{
  if (!m_isDefined) return 0;
  else
    {
      size_t countLocal;
      countLocal = numOutcastLocal();
      
      // gather all countLocals onto process 0
      int srcProc = 0;
      Vector<size_t> allCounts(numProc());
      gather(allCounts, countLocal, srcProc);
      
      size_t count = 0;
      if(procID() == srcProc)
	{
	  for(int ivec = 0; ivec < numProc(); ivec++)
	    {
	      count += allCounts[ivec];
	    }
	}

      // broadcast right answer to all procs
      broadcast(count, srcProc);
      return count;
    }
}

template <class P>
List<P>& ParticleData<P>::outcast()
{
  return m_outcast;
}

template <class P>
bool ParticleData<P>::isClosed() const
{ 
  return m_outcast.isEmpty(); 
}

template <class P>
const ProblemDomain& ParticleData<P>::physDomain() const
{
  return m_physDomain;
}

template <class P>
const RealVect& ParticleData<P>::meshSpacing() const
{
  return m_meshSpacing;
}

template <class P>
const RealVect& ParticleData<P>::origin() const
{
  return m_origin;
}

template <class P>
const int& ParticleData<P>::fixedBoxSize() const
{
  return m_fixedBoxSize;
}

template <class P>
inline void ParticleData<P>::allocateVector()
{
  if (this->m_callDelete == true)
    {
      for (unsigned int i=0; i<this->m_vector.size(); ++i)
	{
	  delete this->m_vector[i];
	  this->m_vector[i] = NULL;
	}
    }

  this->m_callDelete = m_factory.callDelete();

  DataIterator it(this->dataIterator()); int nbox=it.size();
  this->m_vector.resize(it.size(), NULL);
  for(int i=0; i<nbox; i++)
    {
      unsigned int index = it[i].datInd();
      Box abox = this->box(it[i]);
      this->m_vector[index] = m_factory.create(abox, 1, it[i]);
      if (this->m_vector[index] == NULL)
        {
	  MayDay::Error("OutOfMemory in ParticleData::Vector");
        }
    }
}

template <class P>
void ParticleData<P>::gatherOutcast()
{
  BoxLayout grids = getBoxes();
  for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      (*this)[dit].getInvalidDestructive(m_outcast,
					 grids[dit]);
    }
}

template <class P>
void ParticleData<P>::remapOutcast()
{
  CH_TIMERS("ParticleData::remapOutcast"); // barrier before this

  // this proc id
  const unsigned myPID = procID();

  // storage for particles to be delivered to each box/proc
  vector<map<unsigned, List<P> > > pp(numProc());
  CompCounts<P> compCounts(pp);

  BoxLayout grids = getBoxes();
  List<P>& p = m_outcast;

  if (!p.isEmpty())
    {
      // pseudo grid spacing
      const RealVect dx = m_meshSpacing*(Real)m_fixedBoxSize;
      set<IntVect,CompIntVect> bivs;
      map<IntVect,unsigned,CompIntVect> mil;
      multimap<IntVect,boxids,CompIntVect>   mip;
      for (LayoutIterator li = grids.layoutIterator(); li.ok(); ++li)
        {
          // shrink box to get valid region
          const Box sbox = coarsen(grids[li],
                                   m_fixedBoxSize);
          const IntVect  biv=sbox.smallEnd();
#ifdef CH_MPI
          const unsigned pid=grids.procID(li());
          const unsigned idx=grids.index(li());
          if (pid==myPID)
            {
              mil[biv] = idx;
            }
          else
            {
              mip.insert(pair<IntVect,boxids>(biv,boxids(idx,pid)));
            }
          bivs.insert(biv);
#else
          mil[biv] = grids.index(li());
#endif
        }
#ifdef CH_MPI
      // rotate
      {
        set<IntVect,CompIntVect>::iterator it;
        for ( it=bivs.begin() ; it != bivs.end(); ++it )
          {
            const int count=mip.count(*it);
            if (count>1)
              {
                unsigned midn= myPID%count;
                pair<RmIt,RmIt> rit=mip.equal_range(*it);
                RmIt mid=rit.first;
                while (midn-->0) mid++;
                myrotate(rit.first,mid,rit.second);
              }
          }
      }
#endif
      // map particles to local boxes
      map<IntVect,unsigned,CompIntVect>::iterator lit=mil.begin();
      // map particles to pid-boxes
      map<IntVect,boxids,CompIntVect>::iterator pit;

      pair<RmIt,RmIt> rit=mip.equal_range(mip.begin()->first);

      for (ListIterator<P> pi(p); pi.ok(); )
        {
          const IntVect ip=locateBin(pi().position(), dx, m_origin);
          if (ip==lit->first)
            {
              pp[myPID][lit->second].transfer(pi);
            }
          else
            {
              lit=mil.find(ip);
              if (lit!=mil.end())
                {
                  pp[myPID][lit->second].transfer(pi);
                }
#ifdef CH_MPI
              else if (ip==rit.first->first)
                {
                  pit = min_element(rit.first,rit.second,compCounts);
                  pp[pit->second.pid][pit->second.idx].transfer(pi);
		  lit=mil.begin();
                }
              else
                {
                  if (mip.count(ip)>0)
                    {
                      rit=mip.equal_range(ip);
                      pit= min_element(rit.first,rit.second,compCounts);
                      pp[pit->second.pid][pit->second.idx].transfer(pi);
		      lit=mil.begin();
                    }
                  else ++pi;
                }
#else
              else ++pi;
#endif
            }
        }
    }

  // collect local mapping
  map<unsigned,List<P> > lp;
  typename map<unsigned,List<P> >::iterator mi;
  for (mi=pp[myPID].begin(); mi!=pp[myPID].end(); ++mi)
    {
      lp[mi->first].catenate(mi->second);
    }
  pp[myPID].clear();

#ifdef CH_MPI
  // distribute particles
  mpi_scatter_part(lp,pp);
#endif

  // finally assign particles to boxes
  for (DataIterator di(grids); di.ok(); ++di)
    {
      unsigned idx = grids.index(di());
      if (lp.find(idx) != lp.end())
        {
          (*this)[di].listItems().catenate(lp[idx]);
        }
    }
}

template <class P>
IntVect ParticleData<P>::enforcePeriodic(const IntVect& a_index) const
{
  Box domainBox = m_physDomain.domainBox();
  Box cbox = coarsen(domainBox, m_fixedBoxSize);
  IntVect lo = cbox.smallEnd();
  IntVect hi = cbox.bigEnd();
  IntVect dim = hi - lo + 1;
  IntVect new_index = a_index;

  for (int d = 0; d < CH_SPACEDIM; d++)
    {
      if (m_physDomain.isPeriodic())
	{

	  if (a_index[d] < lo[d])
	    {
	      new_index[d] += dim[d];
	    }
      
	  if (a_index[d]  > hi[d])
	    {
	      new_index[d] -= dim[d];
	    }
	}
    }
  return new_index;
}

template <class P>
void ParticleData<P>::fillGhosts(ParticleData<P>& a_particlesWithGhosts,
				 const int a_numGhost) const
{
    CH_TIMERS("ParticleData::fillGhosts");
    BoxLayout grids = getBoxes();
    BoxLayout ghostGrids = a_particlesWithGhosts.getBoxes();

    const Box& domainBox = m_physDomain.domainBox();
    IntVect domLo = domainBox.smallEnd();
    IntVect domHi = domainBox.bigEnd();
    RealVect domLE = ((RealVect) domLo) * m_meshSpacing + m_origin;
    RealVect domRE = ((RealVect) domHi + IntVect::Unit) * m_meshSpacing + m_origin;
    RealVect domainWidth = domRE - domLE;
    
    // Because we are dealing with fixed size BoxLayouts, each box has a unique IntVect index. 
    // Here, we construct a mapping between those IntVect indices, the integer indices, and the
    // process numbers.
    map<IntVect, boxids, CompIntVect> IntVectToBoxIDMap;
    for (LayoutIterator lit = grids.layoutIterator(); lit.ok(); ++lit)
      {
	const Box sbox = coarsen(grids[lit],
				 m_fixedBoxSize);
	
	IntVectToBoxIDMap[sbox.smallEnd()].idx = grids.index(lit());
	IntVectToBoxIDMap[sbox.smallEnd()].pid = grids.procID(lit());
      }

    // define and initialize an iterator over that map
    map<IntVect, boxids, CompIntVect>::iterator IntVectToBoxIDIterator;
    IntVectToBoxIDIterator = IntVectToBoxIDMap.begin();

    // each proc will store the ghost particle lists that it needs to send to each other proc.
    vector<map<unsigned, List<P> > > ghostsToSend(numProc());

    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
      {
	const Box& thisBox = grids[dit];
	const ListBox<P>& thisListBox = (*this)[dit];
	unsigned thisIndex = grids.index(dit());
	
	// theseGhosts stores the particles that will need to be copied and sent elsewhere 
	List<P> theseGhosts;
	thisListBox.getInvalid(theseGhosts, grow(thisBox, -a_numGhost)); 
	const RealVect boxSpacing = m_meshSpacing*(Real)m_fixedBoxSize;

	// keep track of where we are sending each ghost to avoid duplication
	vector<unsigned> alreadySentTo;

	// For each ghost particle, we shift its position in each direction
	// to figure out where it should go.
	for (ListIterator<P> pi(theseGhosts); pi.ok(); ++pi)
	  {
	    alreadySentTo.clear();

	    D_TERM6( for (int i = -1; i < 2; i++) {,
		       for (int j = -1; j < 2; j++) {, 
		         for (int k = -1; k < 2; k++) {,
			   for (int l = -1; l < 2; l++) {,
			     for (int m = -1; m < 2; m++) {,
			       for (int n = -1; n < 2; n++) {)

		    IntVect offset(D_DECL6(i, j, k, l, m, n));
		    RealVect shift = m_meshSpacing * a_numGhost * offset;
		    RealVect shiftedPosition = pi().position() + shift;

		    IntVect srcIntVectIndex;
		    srcIntVectIndex = locateBin(pi().position(),
						boxSpacing,
						m_origin);

		    IntVect destIntVectIndex;
		    destIntVectIndex = locateBin(shiftedPosition, 
						 boxSpacing,
						 m_origin);
	
		    bool doSend = false;
		    if (destIntVectIndex != srcIntVectIndex)
		      {
			doSend = true;
		      }

		    IntVect periodicIntVectIndex;
		    periodicIntVectIndex = enforcePeriodic(destIntVectIndex);
		    
		    IntVectToBoxIDIterator = IntVectToBoxIDMap.find(periodicIntVectIndex);
		    unsigned destPID = IntVectToBoxIDIterator->second.pid;
		    unsigned destBoxID = IntVectToBoxIDIterator->second.idx;

		    // do not send the same particle to given box more than once
		    for (unsigned id=0; id < alreadySentTo.size(); id++)
		      {
		    	if (destBoxID == alreadySentTo[id])
		    	  {
		    	    doSend = false;
		    	  }
		      }

		    if (doSend)
		      {
			alreadySentTo.push_back(destBoxID);
			P p = pi();

			// if the particle is being shipped across the domain due to 
			// periodic boundary conditions, we shift its position here
			// so that a naive distance calculation between the particle 
			// and a point on its box will return the *shortest* distance
			// between them.
			if (periodicIntVectIndex != destIntVectIndex)
			  {
			    RealVect x = pi().position();
			    for (int idim = 0; idim < SpaceDim; idim++)
			      {
				if (periodicIntVectIndex[idim] < destIntVectIndex[idim]) {
				  x[idim] -= domainWidth[idim];
				}
				if (periodicIntVectIndex[idim] > destIntVectIndex[idim]) {
				  x[idim] += domainWidth[idim];
				}
			      }
			    p.setPosition(x);
			  }
			
			ghostsToSend[destPID][destBoxID].add(p);
		      }
		    D_TERM6( }, }, }, }, }, })
	      }
      }

    // after communication, localGhosts will have the ghost particles that belong on each proc
    // for particles going to local boxes, we perform this move here
    unsigned myPID = procID();
    map<unsigned, List<P> > localGhosts;
    typename map<unsigned, List<P> >::iterator ghostMapIterator;
    for (ghostMapIterator = ghostsToSend[myPID].begin(); ghostMapIterator != ghostsToSend[myPID].end(); ++ghostMapIterator)
      {
	unsigned myBoxID = ghostMapIterator->first;
	List<P>& myGhosts = ghostMapIterator->second; 
	localGhosts[myBoxID].catenate(myGhosts);
      }
    ghostsToSend[myPID].clear();

#ifdef CH_MPI
    // distribute particles
    mpi_scatter_part(localGhosts, ghostsToSend);
#endif

    // put the particles, together with Ghosts, in a_particlesWithGhosts
    for (DataIterator dit(grids); dit.ok(); ++dit)
    {
      unsigned idx = grids.index(dit());
      if (localGhosts.find(idx) != localGhosts.end())
	{
	  a_particlesWithGhosts[dit].listItems().join((*this)[dit].listItems());
          a_particlesWithGhosts[dit].listItems().catenate(localGhosts[idx]);
        }
    }
}

template <class P>
bool ParticleData<P>::isDefined() const
{
  return m_isDefined;
}

#include "NamespaceFooter.H"

#endif
