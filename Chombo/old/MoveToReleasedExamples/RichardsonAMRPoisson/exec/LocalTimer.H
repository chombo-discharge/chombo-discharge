#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LOCALTIMER_H_
#define _LOCALTIMER_H_

// if TIMER is not defined, then we just have a stub class.
#ifdef TIMER

extern "C"
{
  // these are for gettimeofday() wall-clock timing
#include <unistd.h>
#include <sys/time.h>


  // obviously, if you don't have PAPI installed on your system
  // (which requires kernel modification) then the counters wont
  // work, so don't compile with -DPAPI
#ifdef PAPI
#include <papi.h>
#endif

}

#ifndef TIMER_COUNTER
#define TIMER_COUNTER 0
#endif

#include <list>
#include <string>
#include "UsingNamespace.H"
using namespace std;


#ifdef CH_MPI
#include "mpi.h"
#endif

#define DIAGNOSTIC false

class LocalTimer
{

public:

  LocalTimer();                             // unmanaged
  LocalTimer(const string&);                // root
  LocalTimer(const string&, LocalTimer&);        // managed
  LocalTimer(const string&, LocalTimer&, bool);  // diagnostic
  ~LocalTimer();

  void setup();
  void start();
  void stop();
  void clear();

  inline double getTimeStampWC();
  double wc_time()
  {
    return accumulated_WCtime;
  }

#ifdef PAPI
  long long int papi_counter1()
  {
    return accumulated_counter1;
  }
  long long int papi_counter2()
  {
    return accumulated_counter2;
  }
#endif

  string Name()
  {
    return timer_name;
  }

  long int Count()
  {
    return count;
  }

  static list<LocalTimer*> LocalTimerList;

  static void LocalTimerSummary(void);

private:

  //bool timer_on;    // State of timer, either on(true) or off(false)
  bool diagnostic;

  double accumulated_WCtime, last_WCtime_stamp;

  // wall-clock timer data
#ifndef CH_MPI
  struct timeval tv;   //  Values from call to gettimeofday
  struct timezone tz;
#endif

  int ID;
  string timer_name;
  LocalTimer& Parent;

  double avgWC,  minWC,  maxWC;

#ifdef PAPI
  long long int values[2];
  long long int accumulated_counter1;
  long long int accumulated_counter2;

#endif

  long long int count;
};


#else


// stub LocalTimer class
#include <string>
#include <iostream>
using namespace std;

#define DIAGNOSTIC false

class LocalTimer
{
public:

  LocalTimer()
  {
  }
  LocalTimer(const string&)
  {
  }
  //LocalTimer(const string&, int)
  {
  }
  LocalTimer(const string&, LocalTimer&)
  {
  }
  LocalTimer(const string&, LocalTimer&, bool)
  {
  }
  ~LocalTimer()
  {
  }

  //static list<LocalTimer*> LocalTimerList;

  inline void setup()
  {
  }
  inline void start()
  {
  }
  inline void stop()
  {
  }
  inline void clear()
  {
  }

  inline double getTimeStampWC()
  {
    return 0.0;
  }
  inline double wc_time()
  {
    return 0.0;
  }

#ifdef PAPI
  inline long long int papi_counter1()
  {
    return 0;
  }
  inline long long int papi_counter2()
  {
    return 0;
  }
#endif

  inline string Name()
  {
    return NULL;
  }

  inline long int Count()
  {
    return 0;
  }

  static void LocalTimerSummary(void)
  {
    int rank;
#ifdef CH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank=0;
#endif
    if (rank==0)
    {
      cout << " LocalTimers Disabled.  No timing output written."  << endl;
    }
  }

};

#endif // TIMER

#endif // TIMER_H
