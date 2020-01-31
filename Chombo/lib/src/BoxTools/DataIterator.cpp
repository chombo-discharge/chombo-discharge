#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ClockTicks.H"
#include "DataIterator.H"
#include "memtrack.H"
#include "NamespaceHeader.H"

#ifdef CH_MPI

bool stupidVerbose = false;

inline unsigned long long  getTimeTDC()
{
  return ch_ticks();
}

#ifdef CH_USE_MEMORY_TRACKING
inline unsigned long long  getPeakTDC()
{
  if(stupidVerbose)
    {
      long long current, peak;
      overallMemoryUsage(current, peak);
      pout() << "ch_memcount = " << ch_memcount << ", current = " << current << ", peak = " << peak << endl;
    }
  return ch_memcount;
}
#else
inline unsigned long long  getPeakTDC()
{
  MayDay::Error("calling memory function when none is available");
  return 0;
}
#endif
DataIterator::
DataIterator(const BoxLayout& plan,
             const int* layoutID)
  :m_layout(plan), m_indices(plan.m_dataIndex), m_current(0)
{
  m_timeEnabled = false;
  m_timeDefined = false;
  m_peakEnabled = false;
  m_peakDefined = false;
}

void 
DataIterator::
clearTime()
{
  //only defines if not defined before
  defineTime();
  for (int ibox = 0; ibox < m_time.size(); ibox++)
    {
      m_time[ibox] = 0;
    }
}

void 
DataIterator::
defineTime()
{
  if (!m_timeDefined)
    {
      m_timeDefined = true;
      m_time.resize(m_layout.size(), 0);
    }
}

void 
DataIterator::
clearPeak()
{
  //only defines if not defined before
  definePeak();
  for (int ibox = 0; ibox < m_time.size(); ibox++)
    {
      m_peak[ibox] = 0;
    }
}

void 
DataIterator::
definePeak()
{
  if (!m_peakDefined)
    {
      m_peakDefined = true;
      m_peak.resize(m_layout.size(), 0);
    }
}

void 
DataIterator::
end()
{
  m_current = m_indices->size();
}

void 
DataIterator::
mergeTime()
{
  int count = m_time.size();
  Vector<unsigned long long> tmp(count);
  MPI_Allreduce(&(m_time[0]),&(tmp[0]), count, MPI_LONG_LONG_INT, MPI_SUM, Chombo_MPI::comm);
  m_time = tmp;
}

void 
DataIterator::
mergePeak()
{
  int count = m_peak.size();
  Vector<unsigned long long> tmp(count);
  MPI_Allreduce(&(m_peak[0]),&(tmp[0]), count, MPI_LONG_LONG_INT, MPI_SUM, Chombo_MPI::comm);
  m_peak = tmp;
}

DataIterator::
DataIterator()
  : m_current(-1)
{
  m_timeEnabled = false;
  m_timeDefined = false;
  m_peakEnabled = false;
  m_peakDefined = false;
}


const DataIndex& 
DataIterator::
operator()() const
{
  CH_assert(ok());
  return m_indices->operator[](m_current);
}

void 
DataIterator::
operator++()
{
  if (m_timeEnabled)
    {
      unsigned long long tdiff = getTimeTDC();
      tdiff -= m_startTime;
      const DataIndex& current = this->operator()();
      m_time[current.m_index] += tdiff;
    }
  if(m_peakEnabled)
    {
      unsigned long long pdiff = getPeakTDC();
      pdiff -= m_startPeak;
      const DataIndex& current = this->operator()();
      m_peak[current.m_index] += pdiff;
    }
  ++m_current;
}

bool 
DataIterator::
ok() const
{
  if (m_timeEnabled)
    {
      DataIterator& nonConstDI = const_cast<DataIterator&>(*this);
      nonConstDI.m_startTime = getTimeTDC();
    }

  if (m_peakEnabled)
    {
      DataIterator& nonConstDI = const_cast<DataIterator&>(*this);
      nonConstDI.m_startPeak = getTimeTDC();
    }
  return m_current < m_indices->size();
}

void 
DataIterator::
reset()
{
  begin();
}

void 
DataIterator::
begin()
{
  m_current=0;
}

#endif

#include "NamespaceFooter.H"
