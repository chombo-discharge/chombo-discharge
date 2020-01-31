#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "RefCellTagger.H"
#include "Box.H"
#include "BoxIterator.H"
#include "RealVect.H"
#include "Ancillae.H"
#include "LoHiCenter.H"
#include "LoHiSide.H"
#include "LGintegrator.H"
#include "SelfGravF_F.H"


RefCellTagger::RefCellTagger()
{
  define();
}

RefCellTagger::RefCellTagger(const ProblemDomain& a_problemDomain,
                             const Real&          a_dx,
                             const int&           a_refRatio)
{
  define(a_problemDomain, a_dx, a_refRatio);
}

//
void RefCellTagger::define()
{
  m_isDefined = false;
  m_useRefineGradient = false;
  m_isRefineGradientSet = false;
  m_useRefineShocks = false;
  m_isRefineShocksSet = false;
  m_useRefineVorticity = false;
  m_isRefineVorticitySet = false;
  m_useRefineOverDensity = false;
  m_isRefineOverDensitySet= false;
  m_useRefineJeansLength  = false;
  m_isRefineJeansLengthSet= false;
  m_useRefineRegion = false;
  m_isRefineRegionSet =  false;
}


void RefCellTagger::define(const ProblemDomain& a_problemDomain,
                           const Real&          a_dx,
                           const int&           a_refRatio)
{
  m_dx = a_dx;
  m_refRatio = a_refRatio;
  m_problemDomain = a_problemDomain;

  m_isDefined = true;
}

void RefCellTagger::setParameters(const bool           a_useRefineGradient,
                                  const Real&          a_gradientThreshold,
                                  const Interval&      a_gradVarInterval,
                                  const bool           a_useRefineOverDensity,
                                  const Real&          a_cellMassThreshold,
                                  const bool           a_useRefineShocks,
                                  const Real&          a_presJumpThreshold,
                                  const bool           a_useRefineVorticity,
                                  const Real&          a_vorticityThreshold,
                                  const bool           a_useRefineJeansLength,
                                  const Real&          a_jeansResolThreshold,
                                  const bool           a_useRefineRegion,
                                  const Vector<Box>&   a_refineRegion,
                                  const Vector<RefineMode>& a_refineMode)
{
  m_gradientThreshold  = a_gradientThreshold;
  m_gradVarInterval    = a_gradVarInterval;
  m_useRefineGradient  = a_useRefineGradient;
  m_isRefineGradientSet= true;

  m_presJumpThreshold= a_presJumpThreshold;
  m_useRefineShocks  = a_useRefineShocks;
  m_isRefineShocksSet= true;

  m_vorticityThreshold  = a_vorticityThreshold;
  m_useRefineVorticity  = a_useRefineVorticity;
  m_isRefineVorticitySet= true;

  m_cellMassThreshold     = a_cellMassThreshold;
  m_useRefineOverDensity  = a_useRefineOverDensity;
  m_isRefineOverDensitySet= true;

  m_jeansResolThreshold   = a_jeansResolThreshold;
  m_useRefineJeansLength  = a_useRefineJeansLength;
  m_isRefineJeansLengthSet= true;

  m_useRefineRegion   = a_useRefineRegion;
  m_refineRegion      = a_refineRegion;
  m_refineMode        = a_refineMode;
  m_isRefineRegionSet = true;
}

// Factory method - this object is its own factory
RefCellTagger* RefCellTagger::new_refCellTagger() const
{
  // Make the new object
  RefCellTagger* retval = static_cast<RefCellTagger*>(new RefCellTagger());

  retval->setParameters(m_useRefineGradient,
                        m_gradientThreshold,
                        m_gradVarInterval,
                        m_useRefineOverDensity,
                        m_cellMassThreshold,
                        m_useRefineShocks,
                        m_presJumpThreshold,
                        m_useRefineVorticity,
                        m_vorticityThreshold,
                        m_useRefineJeansLength,
                        m_jeansResolThreshold,
                        m_useRefineRegion,
                        m_refineRegion,
                        m_refineMode);

  // Return the new object
  return retval;
}

//
void RefCellTagger::setRefineShocks(const bool a_useRefineShocks,
                                    const Real& a_presJumpThreshold)
{
  m_presJumpThreshold = a_presJumpThreshold;
  m_useRefineShocks   = a_useRefineShocks;
  m_isRefineShocksSet = true;
}

//
void RefCellTagger::setRefineVorticity(const bool a_useRefineVorticity,
                                       const Real& a_vorticityThreshold)
{
  m_vorticityThreshold  = a_vorticityThreshold;
  m_useRefineVorticity  = a_useRefineVorticity;
  m_isRefineVorticitySet= true;
}

// set the bool and cell mass threshold for refinement
void RefCellTagger::setRefineOverdense(const bool  a_useRefineOverDensity,
                                       const Real& a_cellMassThreshold)
{
  m_cellMassThreshold      = a_cellMassThreshold;
  m_useRefineOverDensity   = a_useRefineOverDensity;
  m_isRefineOverDensitySet = true;
}

// set the bool and the jeansResol threshold for refinement
void RefCellTagger::setRefineJeans(const bool  a_useRefineJeansLength,
                                   const Real& a_jeansResolThreshold)
{
  m_jeansResolThreshold    = a_jeansResolThreshold;
  m_useRefineJeansLength   = a_useRefineJeansLength;
  m_isRefineJeansLengthSet = true;
}

// set the bool, the vars interval and the gradient threshold for refinement
void RefCellTagger::setRefineGradient(const bool      a_useRefineGradient,
                                      const Real&     a_gradientThreshold,
                                      const Interval& a_gradVarInterval)
{
  m_gradientThreshold   = a_gradientThreshold;
  m_gradVarInterval     = a_gradVarInterval;
  m_useRefineGradient   = a_useRefineGradient;
  m_isRefineGradientSet = true;
}


// setup for refinement of a region: a mode can be union or intersection for
void RefCellTagger::setRefineRegion(const bool            a_useRefineRegion,
                                    const Vector<Box>&     a_refineRegion,
                                    const Vector<RefineMode>& a_refineMode)
{
  m_useRefineRegion = a_useRefineRegion;
  m_refineRegion = a_refineRegion;
  m_refineMode  = a_refineMode;
  m_isRefineRegionSet = true;
}

//
bool RefCellTagger::refineGradient(const std::string a_var)
{
  if (!m_useRefineGradient) return m_useRefineGradient;
  CH_assert(m_gradVarInterval.size()>0);

  if (a_var == "gas")
  {
    return (m_gradVarInterval.begin() < UNUM);
  }
  if (a_var == "phi")
  {
    return m_gradVarInterval.contains(UNUM);
  }
  else
  {
    pout() << " input var " << a_var << '\n';
    MayDay::Warning("RefCellTagger::refineGradient: invalid input string");
  }
  return false;
}


//
IntVectSet RefCellTagger::tagShocks(const LevelData<FArrayBox>& a_U)
{
  // it assumes that ghost cells abutting a coarser lev have been filled in
  CH_assert(m_isRefineShocksSet);

  IntVectSet shockTags;

  // variable used for refinement criterion
  //  Interval presInt(1+SpaceDim,1+SpaceDim);
  //  a_U.exchange(presInt);

  const DisjointBoxLayout& grids = a_U.disjointBoxLayout();

  for (DataIterator di = grids.dataIterator(); di.ok(); ++di)
  {
    const Box& box = grids[di()];
    const Box& ubox = a_U[di].box();
    const Box gbox = box & grow(ubox,-1);

    FArrayBox pjump(box,1);
    FArrayBox wrk(ubox,SpaceDim+1);

    pjump.setVal(zero);
    FORT_SHOCKPRESJUMPF(CHF_FRA1(pjump,0),
                        CHF_FRA(wrk),
                        CHF_CONST_FRA(a_U[di]),
                        CHF_BOX(gbox),
                        CHF_BOX(ubox));

    for (BoxIterator bi(box); bi.ok(); ++bi)
    {
      if (pjump(bi()) > m_presJumpThreshold)
      {
        shockTags |= bi();
      }
    }
  }
  return shockTags;
}

//
IntVectSet RefCellTagger::tagVorticity(const LevelData<FArrayBox>& a_U)
{
  // it assumes that ghost cells abutting a coarser lev have been filled in
  CH_assert(m_isRefineVorticitySet);

  IntVectSet vorticityTags;

  const DisjointBoxLayout& grids = a_U.disjointBoxLayout();
  for (DataIterator di = grids.dataIterator(); di.ok(); ++di)
  {
    const Box& box = grids[di()];
    const Box& ubox = a_U[di].box();
    const Box gbox = box & grow(ubox,-1);

    FArrayBox vorticity(box,1);
    FArrayBox wrk(ubox,SpaceDim);

    vorticity.setVal(zero);
    FORT_VORTICITYF(CHF_FRA1(vorticity,0),
                    CHF_FRA(wrk),
                    CHF_CONST_FRA(a_U[di]),
                    CHF_BOX(gbox),
                    CHF_BOX(ubox));

    for (BoxIterator bi(box); bi.ok(); ++bi)
    {
      if (vorticity(bi()) > m_vorticityThreshold)
      {
        vorticityTags |= bi();
      }
    }
  }
  return vorticityTags;
}

//
IntVectSet RefCellTagger::tagOverdense(const LevelData<FArrayBox>& a_rho)
{
  CH_assert(m_isRefineOverDensitySet);
  const Real dv = pow(m_dx,SpaceDim);

  IntVectSet overDensTags;

  const DisjointBoxLayout& grids = a_rho.disjointBoxLayout();
  for (DataIterator di=grids.dataIterator(); di.ok(); ++di)
  {
    for (BoxIterator bi(grids[di()]); bi.ok(); ++bi)
    {
      const Real m = a_rho[di](bi()) * dv;
      if (m > m_cellMassThreshold)
      {
        overDensTags |= bi();
      }
    }
  }
  return overDensTags;
}


//
IntVectSet RefCellTagger::tagGradient(const LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isRefineGradientSet);

  IntVectSet gradTags;

  // variable used for refinement criterion
  int firstComp = 0;
  int numComp   = 1;

  if (m_gradVarInterval.begin() < UNUM)
  {
    firstComp = m_gradVarInterval.begin();
    numComp   = m_gradVarInterval.size();
  }

  const DisjointBoxLayout& grids = a_U.disjointBoxLayout();
  for (DataIterator di=grids.dataIterator(); di.ok(); ++di)
  {
    const Box& box = grids[di()];
    FArrayBox gradMagFab(box,1);

    for (int var=firstComp; var<numComp; ++var)
    {
      relativeGradient(gradMagFab,a_U[di],var,box);

      for (BoxIterator bi(box); bi.ok(); ++bi)
      {
        if (gradMagFab(bi()) >= m_gradientThreshold)
        {
          gradTags |= bi();
        }
      }
    }
  }
  return gradTags;
}

//
IntVectSet RefCellTagger::tagJeans(const LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isRefineJeansLengthSet && m_isDefined);
    MayDay::Error("RefCellTagger::tagJeans: not defined yet ");

  IntVectSet jeansTags;

  const DisjointBoxLayout& grids = a_U.getBoxes();
  for (DataIterator di=grids.dataIterator(); di.ok(); ++di)
  {
    for (BoxIterator bi(grids[di()]); bi.ok(); ++bi)
    {
      Real v2 = zero;
      for (int dir=0; dir<SpaceDim; ++dir)
      {
        v2 += a_U[di](bi(),UMOMX+dir)*a_U[di](bi(),UMOMX+dir);
      }
      Real ljeans= sqrt(v2/a_U[di](bi(),URHO));

      if (ljeans < m_jeansResolThreshold*m_dx)
      {
        jeansTags |= bi();
      }
    }
  }
  return jeansTags;
}


//
IntVectSet RefCellTagger::tagRegion(const RealVect&          a_low,
                                    const RealVect&          a_high,
                                    const DisjointBoxLayout& a_grids)
{
  CH_assert(m_isDefined);

  IntVectSet regionTags;

  IntVect low, high;

  for (int dir=0; dir<SpaceDim; ++dir)
  {
    int lval = int(floor(a_low[dir]/m_dx));
    low.setVal(dir,lval);

    int hval = int(floor(a_high[dir]/m_dx));
    high.setVal(dir,hval);
  }

  const Box region(low, high);

  for (DataIterator di=a_grids.dataIterator(); di.ok(); ++di)
  {
    const Box box = a_grids[di()];
    regionTags |= (box & region);
  }

  return regionTags;
}

// apply tag cells within a box according to the refinement mode
IntVectSet RefCellTagger::tagRegion(const IntVectSet& a_tags,
                                    const int         a_level,
                                    const DisjointBoxLayout& a_grids)
{
  CH_assert(m_isRefineRegionSet==true);

  // do nothing in this case
  if (a_level>m_refineMode.size()-1) return a_tags;

  //
  IntVectSet regionTags;
  for (DataIterator di=a_grids.dataIterator(); di.ok(); ++di)
  {
    regionTags |= (m_refineRegion[a_level] & a_grids[di()]);
  }

  switch(m_refineMode[a_level])
  {
  case FIX:
    break;
  case AND:
    regionTags &= a_tags;
    break;
  case OR:
    regionTags |= a_tags;
    break;
  default:
    MayDay::Error("RefCellTagger::tagRegion: wrong mode");
  }
  return regionTags;
}


//
void RefCellTagger::relativeGradient(FArrayBox&           a_gradMagFab,
                                     const FArrayBox&     a_varFab,
                                     const int            a_comp,
                                     const Box&           a_box)
{
  CH_assert(a_varFab.interval().contains(a_comp));

  const Box& b = a_box;
  FArrayBox gradFab(b,SpaceDim);

  for (int dir = 0; dir < SpaceDim; ++dir)
  {
    const Box bCenter = b & grow(m_problemDomain,-BASISV(dir));
    const Box bLo     = b & adjCellLo(bCenter,dir);
    const int hasLo = ! bLo.isEmpty();
    const Box bHi     = b & adjCellHi(bCenter,dir);
    const int hasHi = ! bHi.isEmpty();
    FORT_GETRELGRADF(CHF_FRA1(gradFab,dir),
                     CHF_CONST_FRA1(a_varFab,a_comp),
                     CHF_CONST_INT(dir),
                     CHF_BOX(bLo),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(bHi),
                     CHF_CONST_INT(hasHi),
                     CHF_BOX(bCenter));
  }

  FORT_MAGNITUDEF(CHF_FRA1(a_gradMagFab,0),
                  CHF_CONST_FRA(gradFab),
                  CHF_BOX(b));
}

