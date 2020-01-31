#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "EBAMRCNSParams.H"
#include "ParmParse.H"
#include "PolyGeom.H"

#include "NamespaceHeader.H"


ostream&
operator<< (ostream&              a_os,
            const EBAMRCNSParams& a_p)
{
  a_os << "AMRCNSParams:" << endl;

  a_os << "doDiffusion                =" << a_p.m_doDiffusion            << endl;
  a_os << "slipBoundaries             =" << a_p.m_slipBoundaries         << endl;
  a_os << "useMassRedist              =" << a_p.m_useMassRedist          << endl;
  a_os << "doSmushing                 =" << a_p.m_doSmushing             << endl;
  a_os << "variable transport coeffs  =" << a_p.m_variableTransportCoeffs<< endl;
  a_os << "reference temperature      =" << a_p.m_referenceTemperature   << endl;
  a_os << "redistRad                  =" << a_p.m_redistRad              << endl;
  a_os << "tagBufferSize              =" << a_p.m_tagBufferSize          << endl;
  a_os << "verbosity                  =" << a_p.m_verbosity              << endl;
  a_os << "domainLength               =" << a_p.m_domainLength           << endl;
  a_os << "initialDtMultiplier        =" << a_p.m_initialDtMultiplier    << endl;
  a_os << "cfl                        =" << a_p.m_cfl                    << endl;
  a_os << "refineThresh               =" << a_p.m_refineThresh           << endl;
  a_os << "specHeatCv                 =" << a_p.m_specHeatCv             << endl;
  a_os << "thermalCond                =" << a_p.m_thermalCond            << endl;
  a_os << "viscosityMu                =" << a_p.m_viscosityMu            << endl;
  a_os << "viscosityLa                =" << a_p.m_viscosityLa            << endl;      
  a_os << "load balance               =" << a_p.m_loadBalanceType        << endl;      

  return a_os;
}

#include "NamespaceFooter.H"
