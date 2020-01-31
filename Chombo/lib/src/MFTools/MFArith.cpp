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

#include "MFArith.H"
#include "MFAliasFactory.H"
#include <iomanip>
#include "NamespaceHeader.H"

/******/
int
MFArith::
orderScript(int icomp, int inorm, int ncomp)
{
  return icomp + inorm*ncomp;
}
/******/
void
MFArith::
compareError(Vector<Real>&                            a_orders,
             const Vector< LevelData<MFCellFAB>* >&   a_errorFine,
             const Vector< LevelData<MFCellFAB>* >&   a_errorCoar,
             const Vector< DisjointBoxLayout >&       a_gridsFine,
             const Vector< DisjointBoxLayout >&       a_gridsCoar,
             const Vector< Vector<EBISLayout> >&      a_ebislvFine,
             const Vector< Vector<EBISLayout> >&      a_ebislvCoar,
             const Vector<int>&                       a_refRat,
             const Box&                               a_coarseDom,
             int a_testverbosity)
{
  CH_TIME("MFArith::compareError");
  const int ncomp = a_errorFine[0]->nComp();
  const int nnorm = 3;
  Real* normsCoar = new Real[ncomp*nnorm];
  Real* normsFine = new Real[ncomp*nnorm];
  Real* regNormsCoar = new Real[ncomp*nnorm];
  Real* regNormsFine = new Real[ncomp*nnorm];
  Real* irrNormsCoar = new Real[ncomp*nnorm];
  Real* irrNormsFine = new Real[ncomp*nnorm];
  a_orders.resize(ncomp*nnorm, 0.0);
  Real* orders    = &(a_orders[0]);
  for (int icomp = 0; icomp < ncomp; icomp++)
    {
      orders[icomp] = 0.0;
      for (int inorm = 0; inorm < nnorm; inorm++)
        {
          normsCoar[MFArith::orderScript(icomp, inorm, ncomp)] = 0;
          normsFine[MFArith::orderScript(icomp, inorm, ncomp)] = 0;
          regNormsCoar[MFArith::orderScript(icomp, inorm, ncomp)] = 0;
          regNormsFine[MFArith::orderScript(icomp, inorm, ncomp)] = 0;
          irrNormsCoar[MFArith::orderScript(icomp, inorm, ncomp)] = 0;
          irrNormsFine[MFArith::orderScript(icomp, inorm, ncomp)] = 0;
        }
    }

  if (a_testverbosity > 1)
    pout() << "==============================================" << endl;
  for (int icomp = 0; icomp < ncomp; icomp++)
    {
      if (a_testverbosity > 1)
        pout() << "Comparing error in variable  " << icomp << endl;
      if (a_testverbosity > 1)
        pout() << "==============================================" << endl;
      for (int itype = 0; itype < 3; itype++)
        {
          EBNormType::NormMode normtype;
          if (itype == 0)
            {
              normtype = EBNormType::OverBoth;
              if (a_testverbosity > 1)
                pout() << endl << "Using all uncovered cells." << endl  ;
            }
          else if (itype == 1)
            {
              normtype = EBNormType::OverOnlyRegular;
              if (a_testverbosity > 1)
                pout() << endl << "Using only regular cells." << endl ;
            }
          else
            {
              normtype = EBNormType::OverOnlyIrregular;
              if (a_testverbosity > 1)
                pout() << endl << "Using only irregular cells." << endl;
            }

          for (int inorm = 0; inorm <= 2; inorm++)
            {
              if (inorm == 0)
                {
                  if (a_testverbosity > 1)
                    pout() << endl << "Using max norm." << endl;
                }
              else
                {
                  if (a_testverbosity > 1)
                    pout() << endl << "Using L-" << inorm << "norm." << endl;
                }
              Real coarnorm = MFArith::norm(a_errorCoar,
                                            a_gridsCoar,
                                            a_ebislvCoar,
                                            a_refRat,
                                            icomp, inorm, normtype);

              Real finenorm = MFArith::norm(a_errorFine,
                                            a_gridsFine,
                                            a_ebislvFine,
                                            a_refRat,
                                            icomp, inorm, normtype);
              if (a_testverbosity > 1)
                pout() << "Coarse Error Norm = " << coarnorm << endl;
              if (a_testverbosity > 1)
                pout() << "Fine   Error Norm = " << finenorm << endl;
              if (itype == 0)
                {
                  normsCoar[MFArith::orderScript(icomp, inorm, ncomp)] = coarnorm;
                  normsFine[MFArith::orderScript(icomp, inorm, ncomp)] = finenorm;
                }
              else if (itype == 1)
                {
                  regNormsCoar[MFArith::orderScript(icomp,inorm,ncomp)] = coarnorm;
                  regNormsFine[MFArith::orderScript(icomp,inorm,ncomp)] = finenorm;
                }
              else
                {
                  irrNormsCoar[MFArith::orderScript(icomp,inorm,ncomp)] = coarnorm;
                  irrNormsFine[MFArith::orderScript(icomp,inorm,ncomp)] = finenorm;
                }
              if ((Abs(finenorm) > 1.0e-10) && (Abs(coarnorm) > 1.0e-10))
                {
                  Real order = log(Abs(coarnorm/finenorm))/log(2.0);
                  //pout() << "Order of scheme = " << order << endl;
                  if (itype == 0)
                    {
                      orders[MFArith::orderScript(icomp,inorm,ncomp)] = order;
                    }
                }
            }
        }
      if (a_testverbosity > 1)
        pout() << "==============================================" << endl ;;
    }

  //output in latex format to be safe
  int nfine = 2*a_coarseDom.size(0);
  if (a_testverbosity > 0)
    {
      pout() << setw(12)
             << setprecision(6)
             << setiosflags(ios::showpoint)
             << setiosflags(ios::scientific) ;
      for (int inorm = 0; inorm <= 2; inorm++)
        {
          pout() << "\\begin{table}[h]" << endl;
          pout() << "\\begin{center}" << endl;
          pout() << "\\begin{tabular}{|c|c|c|c|c|c|} \\hline" << endl;
          pout() << "Variable & Resolution & Error-All & Order & Error-Regular & Error-Irregular\\\\" << endl;;
          pout() << "\\hline \\hline " << endl;
          for (int icomp = 0; icomp < ncomp; icomp++)
            {
              int iindex = MFArith::orderScript(icomp,inorm,ncomp);
              pout() << "var" << icomp << " &    \t "
                     << "$R_{" << nfine/2 << "}$ & "
                     << setw(12)
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << normsCoar[iindex] << " & "
                     << "            " << " & "
                     << setw(12)
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << regNormsCoar[iindex] << " & "
                     << setw(12)
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << irrNormsCoar[iindex];
              pout() << " \\\\ " << endl;
              pout() << "     &    \t "
                     << "$R_{" << nfine << "}$ & "
                     << setw(12)
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << normsFine[iindex] << " & "
                     << setw(12)
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << orders[iindex] << " & "
                     << setw(12)
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << regNormsFine[iindex] << " & "
                     << setw(12)
                     << setprecision(6)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << irrNormsFine[iindex];
              pout() << " \\\\ " << endl <<   "\\hline"  <<  endl;
            }
          pout() << "\\end{tabular}" << endl;
          pout() << "\\end{center}" << endl;
          pout() << "\\caption{Solution error convergence rates using L-" << inorm << " norm. " << endl;
          pout() << "$h_c = 2 h_f$, $D = " << SpaceDim << "$ }"  << endl;
          pout() << "\\end{table}" << endl;
          pout() << endl << endl;
        }
    }

  delete[] normsCoar;
  delete[] normsFine;
  delete[] regNormsCoar;
  delete[] regNormsFine;
  delete[] irrNormsCoar;
  delete[] irrNormsFine;
}

/******/
Real
MFArith::norm(const Vector< LevelData<MFCellFAB>* >& a_src,
              const Vector< DisjointBoxLayout >&     a_grids,
              const Vector< Vector<EBISLayout> >&    a_ebislv,
              const Vector<int>&                     a_refRatio,
              const int&                             a_comp,
              const int&                             a_pval,
              EBNormType::NormMode a_mode)
{
  CH_TIME("MFArith::norm");
  CH_assert(a_comp >= 0);
  CH_assert(a_pval >= -2);

  Real sum = 0;
  Real volume = 0;
  RealVect dxLev = RealVect::Unit;
  int nlevels = a_src.size();

  Vector< LevelData<EBCellFAB>* > src0(nlevels,NULL),
                                  src1(nlevels,NULL);
  for (int ilev=0; ilev<nlevels; ilev++)
    {
      src0[ilev] = new LevelData<EBCellFAB>();
      src1[ilev] = new LevelData<EBCellFAB>();
    }
  aliasMF(src0, 0, a_src);
  aliasMF(src1, 1, a_src);

  for (int ilev  = 0; ilev < nlevels; ilev++)
    {
      //don't count stuff covered by finer levels
      IntVectSet ivsExclude;
      if (ilev < (nlevels-1))
        {
          //put next finer grids into an IVS
          for (LayoutIterator lit = a_grids[ilev+1].layoutIterator(); lit.ok(); ++lit)
            {
              ivsExclude |= a_grids[ilev+1].get(lit());
            }
          //coarsen back down to this level.
          //ivs will now hold the image of the next finer level
          ivsExclude.coarsen(a_refRatio[ilev]);
        }

      Real sumLev, volLev;
      // Phase 0
      EBArith::volWeightedSum(sumLev, volLev, *(src0[ilev]),  a_grids[ilev],
                              a_ebislv[0][ilev], dxLev, ivsExclude,
                              a_comp, a_pval, a_mode);
      if (a_pval == 0)
        {
          if (volLev > 0)
            {
              sum = Max(sum, sumLev);
            }
        }
      else
        {
          sum += sumLev;
          volume += volLev;
        }
      // Phase 1
      EBArith::volWeightedSum(sumLev, volLev, *(src1[ilev]),  a_grids[ilev],
                              a_ebislv[1][ilev], dxLev, ivsExclude,
                              a_comp, a_pval, a_mode);
      if (a_pval == 0)
        {
          if (volLev > 0)
            {
              sum = Max(sum, sumLev);
            }
        }
      else
        {
          sum += sumLev;
          volume += volLev;
        }

      dxLev /= a_refRatio[ilev];
    }

  for (int ilev=0; ilev<nlevels; ilev++)
    {
      delete src0[ilev];
      delete src1[ilev];
    }

  Real normval;
  if (a_pval == 0)
    {
      normval = sum;
    }
  else if ( (a_pval == 1) || (a_pval == -1) || (a_pval == -2))
    {
      if (volume > 0.0)
        {
          normval = sum/volume;
        }
      else
        {
          normval = 0.0;
        }
    }
  else
    {
      Real denom = a_pval;
      Real exponent = 1.0/denom;
      if (volume > 0.0)
        {
          sum /= volume;
          normval = pow(sum, exponent);
        }
      else
        {
          normval = 0.0;
        }
    }
  return normval;
}

/******/
Real
MFArith::normBoundaryExcluded(const Vector< LevelData<MFCellFAB>* >& a_src,
                              const Vector< DisjointBoxLayout >&     a_grids,
                              const Vector< Vector<EBISLayout> >&    a_ebislv,
                              const ProblemDomain&                   a_coarDomain,
                              const Vector<int>&                     a_refRatio,
                              const int&                             a_coarBoundaryExclusion,
                              const int&                             a_comp,
                              const int&                             a_pval,
                              EBNormType::NormMode a_mode)
{
  CH_TIME("MFArith::normBoundaryExcluded");
  CH_assert(a_comp >= 0);
  CH_assert(a_pval >= -2);
  CH_assert(a_coarBoundaryExclusion >= 0);

  Real sum = 0;
  Real volume = 0;
  RealVect dxLev = RealVect::Unit;
  ProblemDomain domLev(a_coarDomain);
  int nlevels = a_src.size();

  // create ivs for boundary exclusion at each level
  Box coarBox = a_coarDomain.domainBox();
  IntVectSet boundaryIvsLev(coarBox);
  boundaryIvsLev -= coarBox.grow(-a_coarBoundaryExclusion);

  Vector< LevelData<EBCellFAB>* > src0(nlevels,NULL),
                                  src1(nlevels,NULL);
  for (int ilev=0; ilev<nlevels; ilev++)
    {
      src0[ilev] = new LevelData<EBCellFAB>();
      src1[ilev] = new LevelData<EBCellFAB>();
    }
  aliasMF(src0, 0, a_src);
  aliasMF(src1, 1, a_src);

  for (int ilev  = 0; ilev < nlevels; ilev++)
    {
      //don't count stuff covered by finer levels
      IntVectSet ivsExclude;
      if (ilev < (nlevels-1))
        {
          //put next finer grids into an IVS
          for (LayoutIterator lit = a_grids[ilev+1].layoutIterator(); lit.ok(); ++lit)
            {
              ivsExclude |= a_grids[ilev+1].get(lit());
            }
          //coarsen back down to this level.
          //ivs will now hold the image of the next finer level
          ivsExclude.coarsen(a_refRatio[ilev]);
        }

      // add on boundary exclusion
      ivsExclude |= boundaryIvsLev;

      Real sumLev, volLev;
      // Phase 0
      EBArith::volWeightedSum(sumLev, volLev, *(src0[ilev]),  a_grids[ilev],
                              a_ebislv[0][ilev], dxLev, ivsExclude,
                              a_comp, a_pval, a_mode);
      if (a_pval == 0)
        {
          if (volLev > 0)
            {
              sum = Max(sum, sumLev);
            }
        }
      else
        {
          sum += sumLev;
          volume += volLev;
        }
      // Phase 1
      EBArith::volWeightedSum(sumLev, volLev, *(src1[ilev]),  a_grids[ilev],
                              a_ebislv[1][ilev], dxLev, ivsExclude,
                              a_comp, a_pval, a_mode);
      if (a_pval == 0)
        {
          if (volLev > 0)
            {
              sum = Max(sum, sumLev);
            }
        }
      else
        {
          sum += sumLev;
          volume += volLev;
        }

      dxLev /= a_refRatio[ilev];
      domLev.refine(a_refRatio[ilev]);
      boundaryIvsLev.refine(a_refRatio[ilev]);
    }

  for (int ilev=0; ilev<nlevels; ilev++)
    {
      delete src0[ilev];
      delete src1[ilev];
    }

  Real normval;
  if (a_pval == 0)
    {
      normval = sum;
    }
  else if ( (a_pval == 1) || (a_pval == -1) || (a_pval == -2))
    {
      if (volume > 0.0)
        {
          normval = sum/volume;
        }
      else
        {
          normval = 0.0;
        }
    }
  else
    {
      Real denom = a_pval;
      Real exponent = 1.0/denom;
      if (volume > 0.0)
        {
          sum /= volume;
          normval = pow(sum, exponent);
        }
      else
        {
          normval = 0.0;
        }
    }
  return normval;
}

/******/
void
MFArith::normByPhase(Vector<Real>&                          a_phaseNorm,
                     const Vector< LevelData<MFCellFAB>* >& a_src,
                     const Vector< DisjointBoxLayout >&     a_grids,
                     const Vector< Vector<EBISLayout> >&    a_ebislv,
                     const Vector<int>&                     a_refRatio,
                     const int&                             a_comp,
                     const int&                             a_pval,
                     EBNormType::NormMode a_mode)
{
  CH_TIME("MFArith::normByPhase");
  CH_assert(a_comp >= 0);
  CH_assert(a_pval >= -2);

  int numPhase = 2;
  a_phaseNorm.resize(numPhase, 0.);

  Vector<Real> sum, volume;
  sum.resize(numPhase, 0.);
  volume.resize(numPhase, 0.);

  RealVect dxLev = RealVect::Unit;
  int nlevels = a_src.size();

  Vector< Vector< LevelData<EBCellFAB>* > > src;
  src.resize(numPhase);
  for (int iphase=0; iphase<numPhase; iphase++)
    {
      src[iphase].resize(nlevels,NULL);
      for (int ilev=0; ilev<nlevels; ilev++)
        {
          src[iphase][ilev] = new LevelData<EBCellFAB>();
        }
      aliasMF(src[iphase], iphase, a_src);
    }

  for (int ilev  = 0; ilev < nlevels; ilev++)
    {
      //don't count stuff covered by finer levels
      IntVectSet ivsExclude;
      if (ilev < (nlevels-1))
        {
          //put next finer grids into an IVS
          for (LayoutIterator lit = a_grids[ilev+1].layoutIterator(); lit.ok(); ++lit)
            {
              ivsExclude |= a_grids[ilev+1].get(lit());
            }
          //coarsen back down to this level.
          //ivs will now hold the image of the next finer level
          ivsExclude.coarsen(a_refRatio[ilev]);
        }

      for (int iphase=0; iphase<numPhase; iphase++)
        {
          Real sumLev = 0.;
          Real volLev = 0.;
          EBArith::volWeightedSum(sumLev, volLev, *(src[iphase][ilev]),
                                  a_grids[ilev], a_ebislv[iphase][ilev], dxLev,
                                  ivsExclude, a_comp, a_pval, a_mode);
          if (a_pval == 0)
            {
              if (volLev > 0)
                {
                  sum[iphase] = Max(sum[iphase], sumLev);
                }
            }
          else
            {
              sum[iphase] += sumLev;
              volume[iphase] += volLev;
            }
          dxLev /= a_refRatio[ilev];
        }
    }

  for (int iphase=0; iphase<numPhase; iphase++)
    {
      for (int ilev=0; ilev<nlevels; ilev++)
        {
          delete src[iphase][ilev];
        }
    }

  for (int iphase=0; iphase<numPhase; iphase++)
    {
      if (a_pval == 0)
        {
          a_phaseNorm[iphase] = sum[iphase];
        }
      else if ( (a_pval == 1) || (a_pval == -1) || (a_pval == -2))
        {
          if (volume[iphase] > 0.0)
            {
              a_phaseNorm[iphase] = sum[iphase]/volume[iphase];
            }
          else
            {
              a_phaseNorm[iphase] = 0.0;
            }
        }
      else
        {
          Real denom = a_pval;
          Real exponent = 1.0/denom;
          if (volume[iphase] > 0.0)
            {
              sum[iphase] /= volume[iphase];
              a_phaseNorm[iphase] = pow(sum[iphase], exponent);
            }
          else
            {
              a_phaseNorm[iphase] = 0.0;
            }
        }
    }
}
#include "NamespaceFooter.H"
