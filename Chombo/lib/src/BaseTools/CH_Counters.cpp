#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CH_Counters.H"
#include <cstddef>
#include <cstdio>
#include <iostream>
#include "MayDay.H"
#include "ClockTicks.H"
#include "parstream.H"
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif


#include "FortranNameMacro.H"

#include "BaseNamespaceHeader.H"
extern "C"{
extern void FORTRAN_BASENAME(FLOPINIT,flopinit)();
}

long long int ch_counters[NCOUNTERS];

int ch_eventset=0, cacheLevels, cacheSize[3], lineSize[3];


#define BASEWORKINGSET 16
#define MINTRIALS      10
#define MAXTRIALS      250
#define SAMPLES        135
#define ERT_ALIGN      32

unsigned long long int kernelTicks[SAMPLES];
unsigned  int memorySize[SAMPLES];

bool initialized=false;
bool countersAdded = true;

unsigned long long int spTicks;
double spTime;

void streamDump(std::ostream& os)
{
  const unsigned long long int GB = 1<<30;
  
  unsigned int trials = MAXTRIALS;
  for(unsigned int i=0 ; i<SAMPLES ; i++)
    {
      os<<memorySize[i]<<" "<< memorySize[i]*trials*spTicks/((GB*kernelTicks[i])*spTime)<<"\n";
      
      if(trials>MINTRIALS) trials-=2;
    }
}

void CountersInit()
{

  if(initialized) return; 
  initialized=true;
  
  for(int i=0; i<NCOUNTERS; ++i)
    {
      ch_counters[i]=0;
    }

  FORTRAN_BASENAME(FLOPINIT,flopinit)() ;


  /*  
  // measure bandwidth 
  spTicks=ch_ticks();
  spTime = TimerGetTimeStampWC();

  unsigned int W=BASEWORKINGSET;
  unsigned int offsets[SAMPLES];
  for(unsigned int i=0; i<SAMPLES; i++,W=W*1.102)
    {
      kernelTicks[i]=0;
      offsets[i] = 1200*rand()/RAND_MAX;
      //kernelTicks[i] = ~kernelTicks[i];
    }
  double* __restrict A __attribute__((aligned(ERT_ALIGN))); 
  double* __restrict B __attribute__((aligned(ERT_ALIGN))); 
  double* __restrict C __attribute__((aligned(ERT_ALIGN)));
  posix_memalign((void**)&A, ERT_ALIGN, W*sizeof(double)); // N is 1.1 times larger than we need, giving us some array padding.
  posix_memalign((void**)&B, ERT_ALIGN, W*sizeof(double));
  posix_memalign((void**)&C, ERT_ALIGN, W*sizeof(double));
 
  unsigned int tcount=1;
  unsigned int tid=0;
#pragma omp parallel
  { 
#ifdef _OPENMP
    tid=omp_get_thread_num();
#pragma omp single
    {
      tcount=omp_get_num_threads();
      pout() <<"OMP THREADS "<< tcount <<std::endl;
    }
#endif

#pragma omp for schedule(dynamic,10) // spread the "first touch" around the node aggressively
    for(unsigned int p=0; p<W; p++)
      {
        B[p]=1.43;
        C[p]=p*1e-4;
      }

#pragma omp barrier 
  unsigned int trials = MAXTRIALS;
  unsigned int N = BASEWORKINGSET;
  unsigned long long int kticks[SAMPLES];
  for(unsigned int i=0 ; i<SAMPLES ; i++, N=N*1.102)
    {
      unsigned int subN = N/tcount;  // compute loop bounds manually
      unsigned int start=offsets[i] + tid*subN;
      unsigned int end  = start + subN;
      if(tid==tcount-1)
        {
          memorySize[i] = 3*8*(N/tcount)*tcount;
        }
      kticks[i]=ch_ticks();  
      for(unsigned int t=0; t<trials; t++)
        {
          double beta = 1.0000045;
          for(unsigned int p = start; p<end; p++)
            {
              A[p] = beta*B[p]+C[p];
            }
          beta *= (1.0 - 1e-8);
        }
      kticks[i] = ch_ticks() - kticks[i];
#pragma omp barrier

      if(trials>MINTRIALS) trials-=2;// this IF has to match the IF in streamDump below.
    }

  // now need to reduce kticks into kernelTicks
#pragma omp critical
  {
    for(unsigned int i=0; i<SAMPLES; i++)
      {
        //pout()<<" "<<kticks[i];
        kernelTicks[i]+=kticks[i];
      }
    // pout()<<"\n";
  }
  }

  for(int i=0; i<SAMPLES; i++)
    {
      kernelTicks[i]/=tcount;
      //   pout()<<"\n memorySize:"<<memorySize[i]<<" kernelTicks:"<<kernelTicks[i];
    }
 
 

  free(A);  free(B); free(C);
  
  spTicks = ch_ticks() - spTicks;
  spTime  = TimerGetTimeStampWC() - spTime;
  streamDump(pout());
  pout()<<std::endl;
  */
  
#ifdef CH_PAPI

  int retval=0;
  retval = PAPI_library_init(PAPI_VER_CURRENT);

  const PAPI_hw_info_t *hwinfo = NULL;     
  if ((hwinfo = PAPI_get_hardware_info()) == NULL)
    MayDay::Error("could not start PAPI\n");

  cacheLevels = hwinfo->mem_hierarchy.levels;

  printf("%d cache levels \n", cacheLevels);
  for(int i=0; i<cacheLevels; i++)
   {
      cacheSize[i]=hwinfo->mem_hierarchy.level[i].cache[0].size;
      lineSize[i]=hwinfo->mem_hierarchy.level[i].cache[0].line_size;
      printf("  L%d size %d line_size %d \n",i+1, cacheSize[i], lineSize[i]);
   }
    

  int EventSet=PAPI_NULL;
  PAPI_create_eventset(&EventSet);
  //if ((retval = PAPI_add_event( EventSet, PAPI_TOT_INS))  != PAPI_OK )  PFAIL("add event");
  if ((retval = PAPI_add_named_event( EventSet, "PAPI_L1_DCM" )) != PAPI_OK )  countersAdded=false; 
  if ((retval = PAPI_add_named_event( EventSet, "PAPI_L2_DCM" )) != PAPI_OK )  countersAdded=false;   
  if ((retval = PAPI_add_named_event( EventSet, "PAPI_L3_TCM" )) != PAPI_OK )  countersAdded=false;

  PAPI_reset(EventSet);
  PAPI_start(EventSet);

#endif
}

#include "BaseNamespaceFooter.H"
