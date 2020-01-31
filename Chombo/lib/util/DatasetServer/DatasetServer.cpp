#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// -------------------------------------------------------------
// DatasetServer.cpp
// -------------------------------------------------------------

#include <csignal>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cerrno>

#include <iostream>
using std::endl;
using std::cerr;
using std::cout;
using std::cin;
using std::flush;

// important -- xlC did not like the .h on <strstream.h>
#include <strstream>
using std::istrstream;

extern "C"
{
#include <Xm/Xm.h>
#include <Xm/MainW.h>
#include <Xm/PushB.h>
#include <Xm/PushBG.h>
#include <Xm/Form.h>

// eww. yuck.
#include <netinet/in.h>
#include <netdb.h>
#if (defined(CH_Solaris) || defined(CH_IRIX) || defined(CH_IRIX64))
#  include <strings.h>
#endif
#include <fcntl.h>
#include <sys/socket.h>

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
}

#include "Box.H"
#include "FArrayBox.H"
#include "LayoutData.H"
#include "LayoutIterator.H"
// ParmParse is included for the PP_Array and PP_String classes
#include "ParmParse.H"
#include "ArrayViewData.H"

#include "Dataset.H"

const int MAXBUFSIZE  = 1024;
const int PORTOFFSET  = 14159 + SpaceDim;

void
CBSocket(XtPointer client_data, int* source, XtInputId* iid);

bool
ReceiveData(int newsockfd,
            bool& isMulti, bool& isIrregular,
            BaseFab<Real>** newfab,
            LayoutData<BaseFab<Real> >** newfultifab,
            char* formatString, char* dataLabel);

void
DrawDatasheet(int argc, char* argv[], int newsockfd, bool isMulti,
              bool isIrregular,
              BaseFab<Real>* newFab,
              LayoutData<BaseFab<Real> >* newLayoutData,
              const char* formatString, const char* dataLabel);

bool verbose;

// -------------------------------------------------------------
void
PrintUsage(char* arg0)
{
  cerr << "Usage:  " << arg0 << " [-v]" << endl;
  exit(-1);
}

// -------------------------------------------------------------
void
TestCount(int testCount)
{
  if (testCount == 0)
  {
    exit(-3);
  }
}

// -------------------------------------------------------------
void
SendAck(int sockfd, const char* ackString)
{
  if (send(sockfd, ackString, strlen(ackString), 0) < 0)
  {
    perror("Bad server ack send.");
    cerr << "  **** ackString = " << ackString << endl;
    exit(-1);
  }
}

// -------------------------------------------------------------
int
main(int argc, char* argv[])
{
  BaseFab<Real>* newFab = NULL;
  LayoutData<BaseFab<Real> >* newLayoutData = NULL;

  bool loop, receivedData;
  verbose = false;

// valid command-line options:
// -v: verbose mode
// -p n: set socket port offset to n
  int user_port_base = 0;
  for (int count = 1; count < argc; ++count)
  {
    const char* arg = argv[count];
    if (strcmp(arg, "-v") == 0)
    {
      verbose = true;
    }
    else if (strcmp(arg, "-p") == 0)
    {
      user_port_base = atoi(argv[++count]);
    }
  }
  if (verbose)
  {
    cout << "user port base = " << user_port_base << endl;
  }

  // How to contact the server
  // not the best way to find a port, services might be a better method
  int GETUID_SERVER_PORT = 3*getuid() + PORTOFFSET + user_port_base;

  struct sockaddr_in serveraddr;
  struct sockaddr_in clientaddr;
  char formatString[MAXBUFSIZE + 1];
  char dataLabel[MAXBUFSIZE + 1];
  int sockfd, newsockfd, pid, pid2;
  int clientaddrlen = sizeof(clientaddr);
  bool isMulti;
  bool isIrregular;

  signal(SIGCHLD, SIG_IGN);

  // Create a socket
  if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
  {
    perror("Bad server socket create");
    exit(1);
  }
  if (verbose)
  {
    cout << "--- socket created:  sockfd = " << sockfd << endl;
  }

  bzero((char *) &serveraddr, sizeof(struct sockaddr_in));
  serveraddr.sin_family      = AF_INET;
  serveraddr.sin_addr.s_addr = htonl(INADDR_ANY);
  serveraddr.sin_port        = htons(GETUID_SERVER_PORT);

  // Bind the address to this socket to make it known
  if (bind(sockfd, (struct sockaddr *)&serveraddr, sizeof(struct sockaddr_in)) < 0)
  {
    perror("Bad server bind");
    exit(-2);
  }

  listen(sockfd, 8);

  loop = true;
  while (loop)
  {
    int* addrptr = (&clientaddrlen);
    bool extraBool;
    // Allow AIX to use the same call as Linux is using.
    // The bottom one did not compile for me on AIX, while the
    // top one did.  That's all i know!  (ndk)
#if (defined(CH_Linux) || defined(CH_AIX))
    extraBool = ((newsockfd = accept(sockfd, (struct sockaddr *) &clientaddr,
                                     (socklen_t *) addrptr)) < 0);
#else
    extraBool =((newsockfd = accept(sockfd, (struct sockaddr *) &clientaddr,
                                    addrptr)) < 0);
#endif
    if (extraBool)
      {
        if (errno == EINTR)
        {
          cout << "Exiting." << endl;
          exit(0);
        }
        else
        {
          perror("Bad server accept");
          exit(-3);
        }
      }

    pid = fork();

    switch(pid)
    {
    case -1:    // error
      perror("Bad server fork");
      break;

    default:            // parent
      // here the pid is the pid of the child
      close(newsockfd);

      // wait for the child to exit (should be immediate)
      int statloc;
      waitpid(pid, &statloc, 0);

      break;

    case 0:               // child
    {
      close(sockfd);

      // ignore signals to kill child if parent dies
      (void) signal(SIGHUP, SIG_IGN);

      // now fork another child and orphan it
      pid2 = fork();

      switch(pid2)
      {
      case -1:  // error
        perror("Bad server fork (for pid2)");
        break;

      default:            // parent (child)
        // orphan the grandchild
        exit(0);
        break;

      case 0:             // child (grandchild)

        // get data from the dataset client
        receivedData = ReceiveData(newsockfd, isMulti, isIrregular,
                                   &newFab, &newLayoutData,
                                   formatString, dataLabel);

        if ( ! receivedData)
        {
          cerr << "Bad data receive:  exiting." << endl;
          exit(-3);
        }

        // ------------------------------------------------ close the socket
        close(newsockfd);
        loop = false;

        break;

      }  // end switch(pid2)

      break;

    }  // end case 0:

    }  // end switch(pid)
  }  // end while (loop)

  DrawDatasheet(argc, argv, newsockfd, isMulti, isIrregular,
                newFab, newLayoutData,      // only one of these is used
                formatString, dataLabel);
  return 0 ;
}  // end main(...)

// -------------------------------------------------------------
bool
ReceiveData(int newsockfd,
            bool& isMulti,
            bool& isIrregular,
            BaseFab<Real>** newfab,
            LayoutData<BaseFab<Real> >** newlayoutdata,
            char* formatString,
            char* dataLabel)
{
  BaseFab<Real>* newFab = NULL;
  LayoutData<BaseFab<Real> >*  newLayoutData = NULL;
  Vector<Box> newboxes;
//  BoxLayout newBoxLayout;
  PP_Array<BaseFab<Real> *> fabPtrArray;
  char datapipe[MAXBUFSIZE + 1];
  int nElements, count;

  // ----------------------------------------- receive the window label
  if ((count = recv(newsockfd, datapipe, MAXBUFSIZE, 0)) < 0)
  {
    perror("Bad server label recv.");
    return false;
  }
  TestCount(count);
  datapipe[count] = '\0';
  SendAck(newsockfd, "label ack");
  strcpy(dataLabel, datapipe);

  // ---------------------------------------------- receive the format
  if ((count = recv(newsockfd, datapipe, MAXBUFSIZE, 0)) < 0)
  {
    perror("Bad server format recv.");
    return false;
  }
  TestCount(count);
  datapipe[count] = '\0';
  SendAck(newsockfd, "format ack");
  strcpy(formatString, datapipe);

  // ---------------------------------------------- receive isIrregular
  if ((count = recv(newsockfd, datapipe, MAXBUFSIZE, 0)) < 0)
  {
    perror("Bad server isIrregular recv.");
    return false;
  }
  TestCount(count);
  datapipe[count] = '\0';
  SendAck(newsockfd, "isIrregular ack");

  char isIrregularString[8];
  strcpy(isIrregularString, datapipe);
  if (strcmp(isIrregularString, "true") == 0)
  {
    isIrregular = true;
    if (verbose)
    {
      cout << " +++ isIrregular" << endl;
    }
  }
  else
  {
    isIrregular = false;
    if (verbose)
    {
      cout << " +++ ! isIrregular" << endl;
    }
  }

  // ---------------------------------------------- receive isMulti
  if ((count = recv(newsockfd, datapipe, MAXBUFSIZE, 0)) < 0)
  {
    perror("Bad server isMulti recv.");
    return false;
  }
  TestCount(count);
  datapipe[count] = '\0';
  SendAck(newsockfd, "isMulti ack");

  char isMultiString[8];
  strcpy(isMultiString, datapipe);
  if (strcmp(isMultiString, "true") == 0)
  {
    isMulti = true;
    if (verbose)
    {
      cout << " +++ isMulti" << endl;
    }
  }
  else
  {
    isMulti = false;
    if (verbose)
    {
      cout << " +++ ! isMulti" << endl;
    }
  }

  if (isMulti)
  {
    // -------------------------------------------------- receive nElements
    if ((count = recv(newsockfd, datapipe, MAXBUFSIZE, 0)) < 0)
    {
      perror("Bad server fab nElements recv.");
      return false;
    }
    TestCount(count);
    datapipe[count] = '\0';
    SendAck(newsockfd, "nElements ack");
    nElements = atoi(datapipe);

    if (nElements < 1)
    {
      cerr << "Error in datasheet server:  received bad nElements." << endl;
      cerr << "  received nElements = " << nElements << endl;
      return false;
    }

    if (verbose)
    {
      cout << "+++ creating Multiple grid data with nElements = " << nElements << endl;
    }
    newLayoutData = new LayoutData<BaseFab<Real> >();
    fabPtrArray.resize(nElements);
  }
  else
  {
    nElements = 1;
  }

  for (int ne = 0; ne < nElements; ne++)
  {

    if (isMulti)
    {
      if (verbose)
      {
        cout << "+++ receiving data for Multiple grid data element " << ne << endl;
      }
    }

    // -------------------------------------------------- receive the box
    if ((count = recv(newsockfd, datapipe, MAXBUFSIZE, 0)) < 0)
    {
      perror("Bad server fab box recv.");
      return false;
    }
    TestCount(count);
    datapipe[count] = '\0';
    SendAck(newsockfd, "box ack");

    istrstream bufferstream(datapipe, sizeof(datapipe));
    Box receivedBox;
    bufferstream >> receivedBox;
    if ( receivedBox.isEmpty())
    {
      cerr << "Error in datasheet server:  received bad box." << endl;
      cerr << "  received box = " << receivedBox << endl;
      return false;
    }
    if (verbose)
    {
      cout << "receivedBox = " << receivedBox << endl;
    }

    // -------------------------------------------------- receive nVar
    if ((count = recv(newsockfd, datapipe, MAXBUFSIZE, 0)) < 0)
    {
      perror("Bad server fab nVar recv.");
      return false;
    }
    TestCount(count);
    SendAck(newsockfd, "nVar ack");

    datapipe[count] = '\0';
    int nFabComponents = atoi(datapipe);
    if (nFabComponents < 1)
    {
      cerr << "Error in datasheet server:  received bad nVar." << endl;
      cerr << "  received nVar = " << nFabComponents << endl;
      return false;
    }
    if (verbose)
    {
      cout << "nVar = " << nFabComponents << endl;
    }

    // -------------------------------------------------- create the fab
    newFab = new BaseFab<Real>(receivedBox, nFabComponents);

    // -------------------------------------------------- receive fab data
    int totalFabDataBytes = sizeof(Real) * receivedBox.numPts();
    int totalBytesReceived, fabDataBytesRemaining;
    int fabDataBufferSize;
    char *putDataHere, *fabComponentStartingAddress;

    for (int fabComponent = 0; fabComponent < nFabComponents; fabComponent++)
    {
      fabDataBytesRemaining = totalFabDataBytes;
      totalBytesReceived = 0;
      fabComponentStartingAddress = (char *) (newFab->dataPtr(fabComponent));

      while (totalBytesReceived < totalFabDataBytes)
      {
        putDataHere = fabComponentStartingAddress + totalBytesReceived;
        fabDataBufferSize = fabDataBytesRemaining;
        if ((count = read(newsockfd, putDataHere, fabDataBufferSize)) < 0)
        {
          perror("Bad server fab data recv");
          break;
        }
        TestCount(count);

        totalBytesReceived    += count;
        fabDataBytesRemaining -= count;
      } // end while
    }  // end for
    if (verbose)
    {
      cout << "received fab data." << endl;
    }

    // -------------------------------------------------- receive the data ptr
    if ((count = recv(newsockfd, datapipe, MAXBUFSIZE, 0)) < 0)
    {
      perror("Bad server fab data ptr recv.");
      return false;
    }
    TestCount(count);
    datapipe[count] = '\0';
    SendAck(newsockfd, "data ptr ack");

    if (verbose)
    {
      cout << "receivedDataPtr = " << datapipe << endl;
    }

    if (isMulti)
    {
//      newBoxLayout.addBox(receivedBox);
      newboxes.push_back (receivedBox);
      fabPtrArray.set(ne, newFab);
    }
  } // end for (ne...)  (nElements loop)

  // ------------------------------------------------ done receiving data

  Vector<int> assign(newboxes.size(), 0);
  BoxLayout newBoxLayout (newboxes, assign);
  newBoxLayout.close();

  if (isMulti)
  {
    newLayoutData->define(newBoxLayout);
    ArrayViewData data(newLayoutData);
    for (int ne = 0; ne < nElements; ne++)
    {
      data.set(ne, fabPtrArray[ne]);
    }
    if (verbose)
    {
      cout << "+++ BoxLayout = " << endl << newBoxLayout.size() << ":" << endl;
      LayoutIterator lit = newBoxLayout.layoutIterator();
      for ( ; lit.ok(); ++lit)
      {
        LayoutIndex li = lit();
        Box bob = newBoxLayout.get(li);
        cout << "  " << li.intCode() << ": " << bob << endl;
      }
    }
  }

  *newfab = newFab;
  *newlayoutdata = newLayoutData;

  return true;
}  // end ReceiveData

// -------------------------------------------------------------
void
DrawDatasheet(int argc, char* argv[], int newsockfd, bool isMulti,
              bool isIrregular,
              BaseFab<Real>* newFab,
              LayoutData<BaseFab<Real> >* newLayoutData,
              const char* formatString, const char* dataLabel)
{
  XtAppContext  app;
  Widget        wTopLevel;
  Dataset      *dataSheet;
// XtInputId code was commented out in 0.2
  //XtInputId     inputId;

  //inputId = XtAppAddInput(app, newsockfd,
  //(XtPointer) XtInputReadMask, CBSocket, NULL);
  //inputId = XtAppAddInput(app, sockfd,
  //(XtPointer) XtInputReadMask, CBSocket, NULL);

  wTopLevel = XtVaAppInitialize(&app, "DataView", NULL, 0,
                                (int *) &argc, argv, NULL,
                                NULL);

  if (verbose)
  {
    cout << "Creating a new dataset." << endl;
  }
  if (isMulti)
  {
    dataSheet = new Dataset(wTopLevel, newLayoutData, formatString, dataLabel);
  }
  else
  {
    dataSheet = new Dataset(wTopLevel, newFab, formatString, dataLabel,
                            isIrregular);
  }
  XtAppMainLoop(app);
}

// -------------------------------------------------------------
void
CBSocket(XtPointer client_data, int* source, XtInputId* iid)
{
  int intin;
  cout << "_in CBSocket." << endl;
  cout << "source = " << *source << endl;
  cout << "iid    = " << *iid << endl;
  cout << "enter:  " << endl << flush;
  cin >> intin;

  int i = 0;
  PP_Array<int> gridsInLevel;

#if (CH_SPACEDIM == 2)
  IntVect ivlo(8, 3);
  IntVect ivhi(28, 34);
#elif (CH_SPACEDIM == 3)
  IntVect ivlo(8, 3, 4);
  IntVect ivhi(28, 34, 42);
#else
  bogus spacedim;
#endif
  Box fbox(ivlo, ivhi);
  BaseFab<Real> *newFab = new BaseFab<Real>(fbox, 1);
  Real *newFabPtr = newFab->dataPtr();
  for (i = 0; i < newFab->box().numPts(); i++)
  {
    newFabPtr[i] = i;
  }
}
// -------------------------------------------------------------
// -------------------------------------------------------------
