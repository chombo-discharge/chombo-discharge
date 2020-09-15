
// This defines a function used to scatter particles to the proper process

#ifndef _MPI_UTIL_H_
#define _MPI_UTIL_H_

#ifdef CH_MPI
#include <mpi.h>
#include <vector>
#include <typeinfo>
#include <iomanip>
#include <map>
using std::map;

#include "SPMD.H"
#include "REAL.H"
#include "List.H"
#include "RealVect.H"
#include "parstream.H"
#include "MayDay.H"

#include "NamespaceHeader.H"

// a_p is a local container in which to collect the P-objects passed by other
// processes. a_pp is a vector numProc long, whose elements are the containers
// of P-objects to be sent to each other process.
template <class P> void mpi_scatter_part(map<unsigned,List<P> >& a_p,
                                         vector<map<unsigned,List<P> > >& a_pp)
{
  CH_assert(a_pp.size()==numProc());

  const size_t psize = P().size();
  vector<unsigned int> snd_sizes(numProc(),0), rcv_sizes(numProc(),0);

  typename map<unsigned int,List<P> >::iterator vi;
  for (int p=0; p<numProc(); p++)
    {
      // total number of ps
      for (vi=a_pp[p].begin(); vi!=a_pp[p].end(); ++vi)
        snd_sizes[p] += vi->second.length();
      // times p-size
      snd_sizes[p] *= psize;
      // +map keys + list lengths
      snd_sizes[p] += a_pp[p].size() * (sizeof(unsigned int)+sizeof(size_t));
    }

  int err = MPI_Alltoall(&snd_sizes[0],1,MPI_UNSIGNED,
                         &rcv_sizes[0],1,MPI_UNSIGNED,
                         Chombo_MPI::comm);
  if (err != MPI_SUCCESS)
    {
      MayDay::Error("MPI communcation error in mpi_assign_part");
    }
  // sanity check
  CH_assert(rcv_sizes[procID()]==0);

  vector<MPI_Request> rcv_req(numProc()-1);
  vector<char*> rcv_buf(numProc(),NULL);
  int rcv_req_cnt=0;

  for (int p=0; p<numProc(); p++)
    {
      if (rcv_sizes[p] > 0)
        {
          rcv_buf[p] = new char[rcv_sizes[p]];
          MPI_Irecv(rcv_buf[p], rcv_sizes[p], MPI_CHAR, p, 1,
                    MPI_COMM_WORLD, &rcv_req[rcv_req_cnt] );
          rcv_req_cnt++;
        }
    }

  vector<MPI_Request> snd_req(numProc()-1);
  vector<char*> snd_buf(numProc(),NULL);
  int snd_req_cnt=0;

  // pack up data
  for (int p=0; p<numProc(); p++)
    {
      if (snd_sizes[p] > 0)
        {
          snd_buf[p] = new char[snd_sizes[p]];

          char* data = snd_buf[p];
          for (vi=a_pp[p].begin(); vi!=a_pp[p].end(); ++vi)
            {
              *((unsigned int*)data)=vi->first;
              data += sizeof(unsigned int);
              *((size_t*)data)=vi->second.length();
              data += sizeof(size_t);
              for (ListIterator<P> li(vi->second); li.ok(); ++li)
                {
                  li().linearOut((void*)data);
                  data += psize;
                }
            }

          MPI_Isend(snd_buf[p], snd_sizes[p], MPI_CHAR, p, 1,
                    MPI_COMM_WORLD, &snd_req[snd_req_cnt] );
          snd_req_cnt++;
        }
    }

  vector<MPI_Status> rcv_stat(numProc()-1),snd_stat(numProc()-1);
  if (snd_req_cnt>0)
    {
      int err = MPI_Waitall(snd_req_cnt, &snd_req[0], &snd_stat[0]);
      if (err != MPI_SUCCESS)
        {
          MayDay::Error("mpi_scatter_part: send communication failed");
        }
    }

  if (rcv_req_cnt>0)
    {
      int err = MPI_Waitall(rcv_req_cnt, &rcv_req[0], &rcv_stat[0]);
      if (err != MPI_SUCCESS)
        {
          MayDay::Error("mpi_scatter_part: receive communication failed");
        }
    }

  // unpack buffers
  for (int p=0; p<numProc(); ++p)
    {
      P q;
      if (rcv_sizes[p] > 0)
        {
          char* data = rcv_buf[p];
          unsigned int in=0;
          while (in < rcv_sizes[p])
            {
              unsigned int idx=*((unsigned int*)data);
              data += sizeof(unsigned int);
              size_t nps = *((size_t*)data);
              data += sizeof(size_t);
              unsigned int np=0;
              while (np++<nps)
                {
                  q.linearIn((void*)data);
                  a_p[idx].add(q);
                  data += psize;
                }
              in += sizeof(unsigned int)+sizeof(size_t)+psize*nps;
            }
        }
    }

  // delete buffers
  for (int p=0; p<numProc(); ++p)
    {
      a_pp[p].clear();
      if (snd_sizes[p] > 0) delete snd_buf[p];
      if (rcv_sizes[p] > 0) delete rcv_buf[p];
    }
}

#include "NamespaceFooter.H"

#endif // if CH_MPI

#endif
