#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _DATAITERATOR_H_
#define _DATAITERATOR_H_

#include "Vector.H"
#include "DataIndex.H"
#include "BoxLayout.H"
#include "SPMD.H"
#include "LayoutIterator.H"
#include "NamespaceHeader.H"


#ifdef CH_MPI
///An Iterator based on a BoxLayout object.
/**
   An Iterator based on a BoxLayout object.  It does not
   support a dereferencing operation, since it is intended
   to work with all of BoxLayouts, DisjointBoxLayouts, BoxLayoutDatas
   LevelDatas, and any object that is built on top
   of a BoxLayout object.

   DataIterator accesses the data in a BoxLayout-based
   object in a <b> data-parallel </b> manner.
   This means that it skips over entries for Boxes not
   assigned to this processor.  The order of access is not defined.

   In serial execution mode, there is no difference between
   DataIterator and LayoutIterator.
*/
class DataIterator
{
  //Note: parallel version doesn't inherit LayoutIterator, serial version does
public:

  /// a null constructed DataIterator will return false on ok()
  DataIterator();

  DataIterator(const BoxLayout& a_layout)
  {
    *this = a_layout.dataIterator();
  }

  virtual ~DataIterator()
  {}

  /// return the index that this iterator is at
  /** Aborts if the iterator is not ok() */
  const DataIndex& operator()() const ;

  /// return a copy of the index that this iterator is at
  /** Aborts if the iterator is not ok() */
  DataIndex i() const
  {
    return this->operator()();
  }

  DataIndex operator[](int ivec) const
  {
    return (DataIndex)((*m_indices)[ivec]);
  }


  /// move the iterator to the next index in the layout
  virtual void operator++();

  /// move the iterator to the next index in the layout
  void incr()
  {
    ++(*this);
  }

  /// return true if this iterator is still in the layout
  virtual bool ok() const;

  /// initialize this iterator to the first index in the layout
  void begin();

  /// same as begin()
  void reset();

  /// move this iterator to after the last index in the layout
  /** The iterator will be !ok() afterwards. */
  void end();

  int size() const
  {
    return m_indices->size();
  }

  ///functions for using measurements of time and memory for load balancing.
  
  ///sets  m_time values to zero
  void clearTime();

  //sets memory measurement to zero
  void clearPeak();

  ///gets current time data
  Vector<unsigned long long> getTime() const
  {
    return m_time;
  }

  ///gets current memory data
  Vector<unsigned long long> getPeak() const
  {
    return m_peak;
  }

  Vector<Box> getBoxes()
  {
    return m_layout.boxArray();
  }

  ///enables timing.  does not set to zero.  use clear time for that
  void enableTime()
  {
    //only defines if not defined before
    defineTime();
    m_timeEnabled = true;
  }

  ///enables measurements based on how much peak memory moves
  void enablePeak()
  {
    //only defines if not defined before
    definePeak();
    m_peakEnabled = true;
  }

  ///turns off timing
  void disableTime()
  {
    m_timeEnabled = false;
  }

  ///turns off memory measurement
  void disablePeak()
  {
    m_peakEnabled = false;
  }

  /// After you are finished timing your local elements, call mergeTimes to fill-in off-processor times.
  void mergeTime();

  /// After you are finished measuring the memory your local elements, call mergeTimes to fill-in off-processor times.
  void mergePeak();

private:

  friend class BoxLayout;
  friend class DisjointBoxLayout;
  friend class TimedDataIterator;

  DataIterator(const BoxLayout& boxlayout, const int* layoutID);

  BoxLayout  m_layout;
  //  RefCountedPtr<Vector<DataIndex> > m_indices;
protected:
  const Vector<DataIndex>* m_indices;

  int m_current;

  void defineTime();
  void definePeak();

  Vector<unsigned long long> m_time;
  Vector<unsigned long long> m_peak;

  bool m_timeEnabled;
  bool m_timeDefined;

  bool m_peakEnabled;
  bool m_peakDefined;
  unsigned long long m_startTime;
  unsigned long long m_startPeak;
};
#else
// serial version

class DataIterator : public LayoutIterator
{
public:
  virtual ~DataIterator()
  {}

  DataIterator()
  {}

  DataIterator(const BoxLayout& a_layout)
  {
    *this = a_layout.dataIterator();
  }

  /// return the index that this iterator is at
  /** Aborts if the iterator is not ok() */
  const DataIndex& operator()() const
  {
    return (const DataIndex&)(LayoutIterator::operator()());
  };

  /// return a copy of the index that this iterator is at
  /** Aborts if the iterator is not ok() */
  DataIndex i() const
  {
    return this->operator()();
  }

  int size() const
  {
    return this->m_layout.size();
  }

  DataIndex operator[](int ivec) const
  {
    return (DataIndex)((*m_indicies)[ivec]);
  }


  ///load balancing function.  no-op in serial
  void clearTime()
  {
  }

  ///load balancing function.  no-op in serial
  void clearPeak()
  {
  }

  ///load balancing function.  no-op in serial
  Vector<unsigned long long> getTime() const
  {
    Vector<unsigned long long> retval(m_layout.size(), 1.);
    return retval;
  }

  Vector<unsigned long long> getPeak() const
  {
    Vector<unsigned long long> retval(m_layout.size(), 1.);
    return retval;
  }

  Vector<Box> getBoxes()
  {
    return m_layout.boxArray();
  }

  ///load balancing function.  no-op in serial
  void enableTime()
  {
  }

  void enablePeak()
  {
  }

  ///load balancing function.  no-op in serial
  void disableTime()
  {
  }

  void disablePeak()
  {
  }

  ///load balancing function.  no-op in serial
  void mergeTime()
  {
  }

  void mergePeak()
  {
  }
private:
  friend class BoxLayout;
  friend class DisjointBoxLayout;
  friend class TimedDataIterator;

protected:
  DataIterator(const BoxLayout& boxlayout, const int* layoutID)
    :LayoutIterator(boxlayout, layoutID)
  {}
};



#define DATAITERATOR(CLASS, BOXLAYOUT)                \
        DataIterator dit = BOXLAYOUT .dataIterator(); \
        for (dit.begin(); dit.ok(); ++dit)            \
          {                                           \
            DataIndex di = dit();                     \
            MT_BEGIN1(CLASS, DataIndex, di)

#define ENDITERATOR(CLASS)                            \
            MT_END1(CLASS, DataIndex, di)             \
          }

#define DATAITERATOR1(CLASS, BOXLAYOUT, TYPE1, VAL1)     \
        DataIterator dit = BOXLAYOUT .dataIterator();    \
        for (dit.begin(); dit.ok(); ++dit)               \
          {                                              \
            DataIndex di = dit();                        \
            MT_BEGIN2(CLASS, TYPE1, VAL1, DataIndex, di)

#define ENDITERATOR1(CLASS, TYPE1, VAL1)                 \
            MT_END2(CLASS, TYPE1, VAL1, DataIndex, di)   \
          }

#define DATAITERATOR2(CLASS, BOXLAYOUT, TYPE1, VAL1, TYPE2, VAL2)     \
        DataIterator dit = BOXLAYOUT .dataIterator();                 \
        for (dit.begin(); dit.ok(); ++dit)                            \
          {                                                           \
            DataIndex di = dit();                                     \
            MT_BEGIN3(CLASS, TYPE1, VAL1, TYPE2, VAL2, DataIndex, di)

#define ENDITERATOR2(CLASS, TYPE1, VAL1, TYPE2, VAL2)                 \
            MT_END3(CLASS, TYPE1, VAL1, TYPE2, VAL2, DataIndex, di)   \
          }

#endif

#include "NamespaceFooter.H"
#endif
