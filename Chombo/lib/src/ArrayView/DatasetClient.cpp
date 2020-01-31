#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef CH_SPACEDIM

// -------------------------------------------------------------------
//  DatasetClient.cpp
// -------------------------------------------------------------------

#include <cstring>
#include <csignal>
#include <cstdlib>
#include <cstdio>
#include <iostream>

using namespace std;

#include <strings.h>
//XXX// this is not std C++ and should be replaced with <sstream>
//XXX#include <strstream>
#include <sstream>
#include <fcntl.h>
#include <sys/socket.h>
// NOTE: it doesn`t look like errno.h is actually used.
// #include <sys/errno.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
//XXX -- <dbs> why is this here when it is already up above?
//XXX#if (defined(BL_Solaris) || defined(BL_IRIX) || defined(BL_IRIX64))
//XXX#  include <strings.h>
//XXX#endif

#include "Box.H"
#include "BaseFab.H"
#include "SPMD.H"
#include "REAL.H"
#include "FArrayBox.H"
#include "BoxLayout.H"
#include "LayoutData.H"
#include "BoxIterator.H"
#include "LayoutIterator.H"

#include "DatasetClient.H"
#include "NamespaceHeader.H"

const int MAXBUFSIZE  = 1024;
const int PORTOFFSET  = 14159 + SpaceDim;
const char *defaultFormat = "%7.5e";
const char *defaultLabel = " ";

// This allows the user to specify a different port to connect to
int arrayview_user_port_offset = 0;

bool ArrayViewInt(BaseFab<int>* a_intFabPtr)
{
  BaseFab<int>& intfab = *a_intFabPtr;
  BaseFab<Real> realfab(intfab.box(), intfab.nComp());

  for (BoxIterator bit(intfab.box()); bit.ok(); ++bit)
    {
      for (int icomp = 0; icomp < intfab.nComp(); icomp++)
        {
          realfab(bit(), icomp) = Real(intfab(bit(), icomp));
        }
    }

  return ArrayView(&realfab);
}

bool CreateSocket(int& a_newsocket)
{
  int                   a_sockfd;
  struct sockaddr_in    serveraddr;
  char                 *serverhost = "localhost";
  struct hostent       *serverhostp;

  // use to contact the server
  int GETUID_SERVER_PORT = 3*getuid() + PORTOFFSET + arrayview_user_port_offset;

  // create socket
  if ((a_sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
  {
    perror("Bad client socket create");
    return false;
  }
  // cout << "=== after opening socket." << endl;

  // set up the socket structures
  bzero((char *) &serveraddr, sizeof(struct sockaddr_in));
  serveraddr.sin_family = AF_INET;

  if ((serverhostp = gethostbyname(serverhost)) == (struct hostent *) NULL)
  {
    std::cerr << "gethostbyname on " << serverhost << " failed" << std::endl;
    return false;
  }

#ifdef CH_CRAY
  unsigned int saddr = serveraddr.sin_addr.s_addr ;
  bcopy(serverhostp->h_addr, (char *)&saddr,
        serverhostp->h_length);
#else
  bcopy(serverhostp->h_addr, (char *)&(serveraddr.sin_addr.s_addr),
        serverhostp->h_length);
#endif
  serveraddr.sin_port = htons(GETUID_SERVER_PORT);

  // connect to the server
  if (connect(a_sockfd, (sockaddr *)&serveraddr, sizeof(serveraddr)) < 0)
  {
    perror ("Bad client connect");
    return false;
  }
  // cout << "=== connection successful." << endl;

  a_newsocket = a_sockfd;

  return true;
}

bool SendString(int         a_sockfd,
                const char* a_sendstring)
{
  int count;
  char ptrbuffer[MAXBUFSIZE];

  if (send(a_sockfd, a_sendstring, strlen(a_sendstring), 0) < 0)
  {
    perror("Bad client a_sendstring send");
    return false;
  }

  // wait for acknowledgment
  if ((count = recv(a_sockfd, ptrbuffer, MAXBUFSIZE, 0)) < 0)
  {
    perror("Bad a_sendstring ack.");
    return false;
  }

  ptrbuffer[count] = '\0';
  // cout << "<<< received a_sendstring ack:  " << ptrbuffer << endl;

  return true;
}

bool SendRealArray(int        a_sockfd,
                   Real*      a_data[],
                   int        a_nvar,
                   const int* a_lodim,
                   const int* a_hidim)
{
  int  count;
  char buffer[MAXBUFSIZE];
//XXX  char ptrbuffer[MAXBUFSIZE];

  IntVect ivlo(a_lodim);
  IntVect ivhi(a_hidim);
  Box dataBox(ivlo, ivhi);

  // --------------------------------------------------- send the box
  // cout << ">>> sending box." << endl;
//XXX  std::ostrstream bufferstream(buffer, sizeof(buffer));
  std::ostringstream bufferstream ;

  bufferstream << dataBox ;

//XXX  if (send(a_sockfd, buffer, strlen(buffer), 0) < 0)
  if (send(a_sockfd, bufferstream.str().c_str(), bufferstream.str().length(), 0) < 0)
  {
    perror("Bad client box send");
    return false;
  }

  // wait for acknowledgment
  if ((count = recv(a_sockfd, buffer, MAXBUFSIZE, 0)) < 0)
  {
    perror("Bad box ack.");
    return false;
  }
  buffer[count] = '\0';
  // cout << "<<< received box ack:  " << buffer << endl;

  // --------------------------------------------------- send nVar
  // cout << ">>> sending nVar." << endl;
//XXX  sprintf(buffer, "%d", a_nvar);
  bufferstream.str("");
  bufferstream << a_nvar ;

//XXX  if (send(a_sockfd, buffer, strlen(buffer), 0) < 0)
  if (send(a_sockfd, bufferstream.str().c_str(), bufferstream.str().length(), 0) < 0)
  {
    perror("Bad client nVar send");
    return false;
  }

  // wait for acknowledgment
  if ((count = recv(a_sockfd, buffer, MAXBUFSIZE, 0)) < 0)
  {
    perror("Bad nVar ack.");
    return false;
  }
  buffer[count] = '\0';
  // cout << "<<< received nVar ack:  " << buffer << endl;

  // --------------------------------------------------- send the data.
  // cout << ">>> sending data." << endl;

  int totalDataBytes = sizeof(Real) * dataBox.numPts();
  int totalBytesSent, dataBytesRemaining;
  int dataBufferSize;
  char *getDataHere, *dataComponentStartingAddress;

  for (int dataComponent = 0; dataComponent < a_nvar; dataComponent++)
  {
    // cout << "dataComponent = " << dataComponent << endl;
    totalBytesSent = 0;
    dataBytesRemaining = totalDataBytes;
    dataComponentStartingAddress = (char *) (a_data[dataComponent]);

    // send a chunk of data
    while (totalBytesSent < totalDataBytes)
    {
      getDataHere = dataComponentStartingAddress + totalBytesSent;
      dataBufferSize = dataBytesRemaining;

      if ((count = write(a_sockfd, getDataHere, dataBufferSize)) < 0)
      {
        perror("Bad client data send");
        return false;
      }
      // cout << "  bytes sent = " << count << endl;
      totalBytesSent     += count;
      dataBytesRemaining -= count;
    }  // end while
  }  // end for

  // --------------------------------------------------- send the pointer

//XXX  std::ostrstream ptrbufferstream(ptrbuffer, sizeof(ptrbuffer));
  std::ostringstream ptrbufferstream ;
  ptrbufferstream << a_data[0] ;

//XXX  if (send(a_sockfd, ptrbuffer, strlen(ptrbuffer), 0) < 0)
  if (send(a_sockfd, ptrbufferstream.str().c_str(), ptrbufferstream.str().length(),  0) < 0)
  {
    perror("Bad client data ptr send");
    return false;
  }

  // wait for acknowledgment
//XXX  if ((count = recv(a_sockfd, ptrbuffer, MAXBUFSIZE, 0)) < 0)
  if ((count = recv(a_sockfd, buffer, MAXBUFSIZE, 0)) < 0)
  {
    perror("Bad data ptr ack.");
    return false;
  }
//XXX  ptrbuffer[count] = '\0';
  buffer[count] = '\0';
  // cout << "<<< received data ptr ack:  " << ptrbuffer << endl;

  // --------------------------------------------------- done sending data

  return true;
} // end SendRealArray

// -------------------------------------------------------------------
// pointer to BaseFab interface
// -------------------------------------------------------------------
bool ArrayView(BaseFab<Real>* a_debugFab)
{
  return (ArrayViewFabFormatLabel(a_debugFab, defaultFormat, "Fab"));
}

bool ArrayViewFabFormatLabel(BaseFab<Real>* a_debugFab,
                             const char*    a_format,
                             const char*    a_label)
{
  bool returnValue;
  int nvar = a_debugFab->nComp();

  if (nvar < 1)
  {
    std::cerr << "Error in ArrayView:  fab nVar < 1:  fab->nVar = " << nvar << std::endl;
    return false;
  }

  if (a_debugFab->box().isEmpty())
  {
    std::cerr << "Error in ArrayView:  bad fab box = " << a_debugFab->box() << std::endl;
    return false;
  }

  Real** dataArray = new Real* [nvar];
  // build the array of real pointers
  for (int d = 0; d < nvar; d++)
  {
    dataArray[d] = a_debugFab->dataPtr(d);  // don't assume contiguous
  }

  returnValue = ArrayViewRealPtrArrayNVarDims(dataArray, nvar,
                                              a_debugFab->box().smallEnd().getVect(),
                                              a_debugFab->box().bigEnd().getVect(), a_format, a_label);

  delete [] dataArray;
  return returnValue;
}

// -------------------------------------------------------------------
// pointer to LayoutData interface
// -------------------------------------------------------------------
bool ArrayViewLDF(LevelData<FArrayBox>* a_debugLevelData)
{
  return (ArrayViewLDFFormatLabel(a_debugLevelData,
                                  defaultFormat,
                                  "LevelData<FArrayBox>"));
}

bool MultiArrayViewFab(LayoutData<FArrayBox>* a_debugLayoutData)
{
  return (MultiArrayViewFormatLabel((LayoutData<BaseFab<Real> >*) a_debugLayoutData,
                                    defaultFormat,
                                    "LayoutData<FArrayBox>"));
}

bool MultiArrayView(LayoutData<FArrayBox>* a_debugLayoutData)
{
  return (MultiArrayViewFab(a_debugLayoutData));
}

bool MultiArrayViewRealBaseFab(LayoutData<BaseFab<Real> >* a_debugLayoutData)
{
  return (MultiArrayViewFormatLabel(a_debugLayoutData,
                                    defaultFormat,
                                    "LayoutData<BaseFab<Real>>"));
}

// -------------------------------------------------------------------
// pointer to real interface
// -------------------------------------------------------------------
bool ArrayViewRealPtrArrayNVarDims(Real*       a_data[],
                                   int         a_nvar,
                                   const int*  a_lodim,
                                   const int*  a_hidim,
                                   const char* a_format,
                                   const char* a_label)
{
  int sockfd;

  if (!CreateSocket(sockfd))
  {
    return false;
  }

  // --------------------------------------------------- send data label
  if (!SendString(sockfd, a_label))
  {
    return false;
  }

  // --------------------------------------------------- send format
  if (!SendString(sockfd, a_format))
  {
    return false;
  }

  // --------------------------------------------------- send isIrregular
  if (!SendString(sockfd, "false")) // is not irregular
  {
    return false;
  }

  // --------------------------------------------------- send isMultiFab
  if (!SendString(sockfd, "false")) // not a MultiFab
  {
    return false;
  }

  // --------------------------------------------------- send nElements
  // don't send nElements

  // --------------------------------------------------- send the data
  return (SendRealArray(sockfd, a_data, a_nvar, a_lodim, a_hidim));

} // end of function

bool MultiArrayViewFormatLabel(LayoutData<BaseFab<Real> >* a_layoutdata,
                               const char*                 a_format,
                               const char*                 a_label)
{
  int sockfd;
  char buffer[MAXBUFSIZE];

  if (!CreateSocket(sockfd))
  {
    return false;
  }

  // --------------------------------------------------- send data label
  if (! SendString(sockfd, a_label))
  {
    return false;
  }

  // --------------------------------------------------- send format
  if (a_format == NULL)
  {
    if (!SendString(sockfd, defaultFormat))
    {
      return false;
    }
  }
  else
  {
    if (!SendString(sockfd, a_format))
    {
      return false;
    }
  }

  // --------------------------------------------------- send isIrregular
  if (!SendString(sockfd, "false")) // is not irregular
  {
    return false;
  }

  // --------------------------------------------------- send isMulti
  if (!SendString(sockfd, "true")) // this has multiple grids
  {
    return false;
  }

  // --------------------------------------------------- send nElements
  int num_elems = a_layoutdata->boxLayout().numBoxes(procID());
  // cout << ">>> sending nElements = " << num_elems << endl;
  sprintf(buffer, "%d", num_elems);

  if (!SendString(sockfd, buffer))
  {
    return false;
  }

  // ArrayViewData data(a_layoutdata);

  // --------------------------------------------------- send the data
  for (DataIterator it = a_layoutdata->dataIterator(); it.ok(); ++it)
  {
    // construct dataArray for this element
    BaseFab<Real>& fab = a_layoutdata->operator[](it());
    int nvar = fab.nComp();
    Real** dataArray = new Real * [nvar];

    for (int d = 0; d < nvar; d++) // build the array of Real *
    {
      dataArray[d] = fab.dataPtr(d);  // don't assume contiguous
    }

    int lo_vect[SpaceDim];
    int hi_vect[SpaceDim];

    for (int d = 0; d < SpaceDim; ++d)
    {
      lo_vect[d] = fab.box().smallEnd(d);
      hi_vect[d] = fab.box().bigEnd(d);
    }

    /*
      if (! SendRealArray(sockfd, dataArray, nvar,
                          (fab.box()).loVect(), (fab.box()).hiVect()))
    */
    if (! SendRealArray(sockfd, dataArray, nvar, lo_vect, hi_vect) )
    {
      return false;
    }

    delete [] dataArray;
  }

  return true;
} // end of function

bool ArrayViewLDFFormatLabel(LevelData<FArrayBox>* a_debugLevelData,
                             const char*           a_format,
                             const char*           a_label)
{
  DisjointBoxLayout inputDBL = a_debugLevelData->getBoxes();
  Vector<Box> boxes;

  LayoutIterator lit = inputDBL.layoutIterator();
  for (lit.reset(); lit.ok(); ++lit)
  {
    boxes.push_back(inputDBL.get(lit()));
  }

  //assign all boxes to proc 0
  Vector<int> assign(boxes.size(), 0);
  DisjointBoxLayout tempDBL(boxes, assign);

  tempDBL.close();
  LevelData<FArrayBox> tempLDF(tempDBL,
                               a_debugLevelData->nComp(),
                               a_debugLevelData->ghostVect());

  DataIterator dit = tempLDF.dataIterator();
  for (dit.reset(); dit.ok(); ++dit)
  {
    tempLDF[dit()].setVal(0.);
  }

  a_debugLevelData->copyTo(a_debugLevelData->interval(),
                         tempLDF,
                         tempLDF.interval());

  LevelData<FArrayBox>* leveldata = &tempLDF;

  if (procID() == 0)
  {
    int sockfd;
    char buffer[MAXBUFSIZE];

    if (!CreateSocket(sockfd))
    {
      return false;
    }

    // --------------------------------------------------- send data label
    if (! SendString(sockfd, a_label))
    {
      return false;
    }

    // --------------------------------------------------- send format
    if (a_format == NULL)
    {
      if ( ! SendString(sockfd, defaultFormat))
      {
        return false;
      }
    }
    else
    {
      if ( ! SendString(sockfd, a_format))
      {
        return false;
      }
    }

    // --------------------------------------------------- send isIrregular
    if (!SendString(sockfd, "false")) // is not irregular
    {
      return false;
    }

    // --------------------------------------------------- send isMulti
    if (!SendString(sockfd, "true")) // this has multiple grids
    {
      return false;
    }

    // --------------------------------------------------- send nElements
    int num_elems = leveldata->boxLayout().numBoxes(procID());
    // cout << ">>> sending nElements = " << num_elems << endl;
    sprintf(buffer, "%d", num_elems);

    if (!SendString(sockfd, buffer))
    {
      return false;
    }

    // ArrayViewData data(leveldata);

    // --------------------------------------------------- send the data
    for (DataIterator it = leveldata->dataIterator(); it.ok(); ++it)
    {
      // construct dataArray for this element
      BaseFab<Real>& fab = leveldata->operator[](it());
      int nvar = fab.nComp();
      Real** dataArray = new Real * [nvar];

      for (int d = 0; d < nvar; d++)    // build the array of Real *
      {
        dataArray[d] = fab.dataPtr(d);  // don't assume contiguous
      }

      int lo_vect[SpaceDim];
      int hi_vect[SpaceDim];

      for (int d = 0; d < SpaceDim; ++d)
      {
        lo_vect[d] = fab.box().smallEnd(d);
        hi_vect[d] = fab.box().bigEnd(d);
      }

      /*
        if (!SendRealArray(sockfd, dataArray, nvar,
                           (fab.box()).loVect(), (fab.box()).hiVect()))
      */
      if (!SendRealArray(sockfd, dataArray, nvar, lo_vect, hi_vect))
      {
        return false;
      }

      delete [] dataArray;
    }
  }

  return true;
} // end of function
#include "NamespaceFooter.H"

#endif // CH_SPACEDIM
