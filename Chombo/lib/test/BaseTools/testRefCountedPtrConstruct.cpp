#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include <cstring>

// Define this to get more output from RefCountedPtr
//#define DEBUGRCP

#include "REAL.H"
#include "RefCountedPtr.H"
#include "parstream.H"
#include "SPMD.H"

#include "UsingBaseNamespace.H"

using std::endl;

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testRCPConstruct();

/// Global variables for handling output:
static const char *pgmname = "testRefCountedPtrConstruct" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;

  if ( verbose )
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int ret = testRCPConstruct() ;
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}

class A
{
public:
  A(int a_ai)
  :
  ai(a_ai)
    { }
  virtual ~A()
    { }
  virtual int geti() const
    {
      return ai;
    }
  virtual void print() const
    {
      pout() << "A: " << ai << endl;
    }
  int ai;
};

class B : public A
{ 
public:
  B(int a_bi)
  :
  A(a_bi+1),
  bi(a_bi)
    { }
  virtual ~B()
    { }
  virtual int geti() const
    {
      return bi;
    }
  virtual void print() const
    {
      pout() << "B: " << bi << ' ' << ai << endl;
    }
  int bi;
};

// The tests include constructions that should work and also some that should
// fail.  The failures are of course commented out but you can uncomment them
// to verify that they do in fact fail (some at compiler time, others at
// runtime)
int
testRCPConstruct()
{
  int status = 0;


//--Construction (same pointer, tests qualifiers)

  if (verbose)
    {
      pout() << "\n-1 - Regular class pointer (ok)\n";
    }
  RefCountedPtr<A> a(new A(1));
  if (a->geti() != 1) ++status;
  // Note - the following should not be ambiguous
  if (verbose)
    {
      pout() << "\n-2 - Construct non-const from non-cost pointer (ok)\n";
    }
  RefCountedPtr<A> a2(a);
  if (a2->geti() != 1) ++status;
  if (verbose) 
    {
      pout() << "\n-3 - Construct const from non-const (ok)\n";
    }
  RefCountedPtr<const A> a3(a);
  if (a3->geti() != 1) ++status;
  // if (verbose)
  //   {
  //     pout() << "\n-4- Construct non-const from const (bad - fails ";
  //     "compilation)\n";
  //   }
  // RefCountedPtr<A> a4(a3);
  // a4->print();
  if (verbose) 
    {
      pout() << "\n-5 - Derived class held in base class pointer (ok)\n";
    }
  RefCountedPtr<A> b(new B(2));
  if (b->geti() != 2) ++status;
  if (b->A::geti() != 3) ++status;
  if (verbose) 
    {
      pout() << "\n-6 - Construct const from non-const (ok)\n";
    }
  RefCountedPtr<const A> b2(b);
  if (b2->geti() != 2) ++status;
  if (b2->A::geti() != 3) ++status;

//--Assignment (same pointer, tests qualifiers)

  if (verbose)
    {
      pout() << "\n-7 - Assigning const from non-const (ok)\n";
    }
  b2 = b;
  if (b2->geti() != 2) ++status;
  if (b2->A::geti() != 3) ++status;
  // if (verbose)
  //   {
  //     pout() << "\n-8 - Assigning non-const from const (bad - fails "
  //       "complilation)\n";
  //   }
  // b = b2;
  // b->print();

//--Conversion to POD pointer (tests qualifiers)

  if (verbose)
    {
      pout() << "\n-9 - Auto-conversion to POD pointer (ok)\n";
    }
  A* a5 = b;
  if (a5->geti() != 2) ++status;
  if (verbose)
    {
      pout() << "\n-10- Auto-conversion to const POD pointer (ok)\n";
    }
  const A* a6 = b;
  if (a6->geti() != 2) ++status;
  if (verbose)
    {
      pout() << "\n-11- Auto-conversion of const RC to const POD pointer "
        "(ok)\n";
    }
  const A* a7 = b2;
  if (a7->geti() != 2) ++status;
  // if (verbose)
  //   {
  //     pout() << "\n-12- Auto-conversion of non-const RC to const POD "
  //       "pointer (bad - fails compilation)\n";
  //   }
  // A* a8 = b2;
  // a8->print();

//--Conversions of RC ptrs with related POD pointers through construction

  if (verbose)
    {
      pout() << "\n-13- Conversion to derived through constructor (ok)\n";
    }
  RefCountedPtr<B> bb3(b);
  if (bb3->geti() != 2) ++status;
  if (verbose)
    {
      pout() << "\n-14- Conversion to const derived through constructor (ok)\n";
    }
  RefCountedPtr<const B> bb4(b);
  if (bb4->geti() != 2) ++status;
  // if (verbose)
  //   {
  //     pout() << "\n-15- Conversion of const base to non-const derived "
  //       "through constructor (bad - fails compilation)\n";
  //   }
  // RefCountedPtr<B> bb5(b2);
  // bb5->print();
  // if (verbose)
  //   {
  //     pout() << "\n-16- Conversion of actual base to derived through "
  //       "constructor (bad - fails runtime assertion)\n";
  //   }
  // RefCountedPtr<B> bb6(a);
  // bb6->print();
  // if (verbose)
  //   {
  //     pout() << "\n-17- Conversion of actual const base to non-const "
  //       "derived through constructor (bad - fails compilation)\n";
  //   }
  // RefCountedPtr<B> bb7(a3);
  // bb7->print();
  if (verbose)
    {
      pout() << "\n-18- Conversion to base through constructor (ok)\n";
    }
  RefCountedPtr<A> b4(bb3);
  if (b4->geti() != 2) ++status;
  static_cast<RefCountedPtr<A> >(bb3);  // Casting also works...
  if (verbose)
    {
      pout() << "\n-19- Conversion to const base through constructor (ok)\n";
    }
  RefCountedPtr<const A> b5(bb3);
  if (b5->geti() != 2) ++status;
  // if (verbose)
  //   {
  //     pout() << "\n-20- Conversion of const derived to non-const base "
  //       "through constructor (bad - fails compilation)\n";
  //   }
  // RefCountedPtr<A> b6(bb4);
  // b6->print();

//--Conversions of RC ptrs with related POD pointers during assignment

  if (verbose)
    {
      pout() << "\n-21- Assignment of base to derived (ok)\n";
    }
  RefCountedPtr<B> bb5(new B(3));
  RefCountedPtr<B> bb6(bb5);
  bb6 = b;
  if (bb6->geti() != 2) ++status;
  // if (verbose)
  //   {
  //     pout() << "\n-22- Assignment to constant RC (bad - fails compilation "
  //       "but very noisy)\n";
  //   }
  // const RefCountedPtr<B> bb7(bb5);
  // bb7 = b;
  // bb7->print();
  if (verbose)
    {
      pout() << "\n-23- Assignment of non-const base to const derived (ok)\n";
    }
  RefCountedPtr<const B> bb8(bb5);
  bb8 = b;
  if (bb8->geti() != 2) ++status;
  // if (verbose)
  //   {
  //     pout() << "\n-24- Assignment of const base to non-const derived "
  //       "(bad - fails compilation)\n";
  //   }
  // RefCountedPtr<B> bb9(bb5);
  // bb9 = b2;
  // bb9->print();
  // if (verbose)
  //   {
  //     pout() << "\n-25- Assignment of actual base to derived (bad - fails "
  //       "runtime assertion)\n";
  //   }
  // RefCountedPtr<B> bb10(bb5);
  // bb10 = a;
  // bb10->print();
  // if (verbose)
  //   {
  //     pout() << "\n-26- Assignment of actual const base to non-const "
  //       "derived (bad - fails compilation)\n";
  //   }
  // RefCountedPtr<B> bb11(bb5);
  // bb11 = a3;
  // bb11->print();
  if (verbose)
    {
      pout() << "\n-27- Assignment of derived to base (ok)\n";
    }
  RefCountedPtr<A> b12(bb5);
  b12 = bb3;
  if (b12->geti() != 2) ++status;
  if (verbose)
    {
      pout() << "\n-28- Assignment of non-const derived to const base (ok)\n";
    }
  RefCountedPtr<const A> b13(bb5);
  b13 = bb3;
  if (b13->geti() != 2) ++status;
  // if (verbose)
  //   {
  //     pout() << "\n-29- Assignment of const derived to non-const base "
  //       "(bad - fails compilation)\n";
  //   }
  // RefCountedPtr<A> b14(bb5);
  // b14 = bb4;
  // b14->print();

//--Array construction

  if (verbose)
    {
      pout() << "\n-30- Array (ok)\n";
    }
  RefCountedPtr<int, RCPArrayPolicy> r1(new int[4]);
  for (int i = 0; i != 4; ++i) r1[i] = 0;
  r1[1] = 1;
  if (r1[1] != 1) ++status;
  if (verbose)
    {
      pout() << "\n-31- Construct const array from non-const (ok)\n";
    }
  RefCountedPtr<const int, RCPArrayPolicy> r2(r1);
  // This would fail since r2.operator[]() returns type 'const int&'
  // r2[1] = 2;
  if (r1[1] != 1) ++status;
  // if (verbose)
  //   {
  //     pout() << "\n-32- Construct non-const array from const (bad - fails "
  //       "compilation)\n";
  //   }
  // RefCountedPtr<int, RCPArrayPolicy> r3(r2);
  // if (verbose)
  //   {
  //     pout() << "\n-33- Construct RC pointer from RC array (bad - fails "
  //       "compilation)\n";
  //   }
  // RefCountedPtr<int> r4(r1);
  // if (verbose)
  //   {
  //     pout() << "\n-34- Construct const RC pointer from RC array (bad - "
  //       "fails compilation)\n";
  //   }
  // RefCountedPtr<const int> r5(r1);

//--Array assignment

  if (verbose)
    {
      pout() << "\n-35- Array assignment (ok)\n";
    }
  RefCountedPtr<int, RCPArrayPolicy> r6(new int[4]);
  for (int i = 0; i != 4; ++i) r6[i] = -1;
  r6 = r1;
  if (r6[1] != 1) ++status;
  if (verbose)
    {
      pout() << "\n-36- Array assign non-const to const (ok)\n";
    }
  RefCountedPtr<const int, RCPArrayPolicy> r7(new int[4]);
  r7 = r1;
  if (r7[1] != 1) ++status;
  // if (verbose)
  //   {
  //     pout() << "\n-37- Array assign const to non-const (bad - fails "
  //       "compilation)\n";
  //   }
  // r6 = r2;
  // if (verbose)
  //   {
  //     pout() << "\n-38- Assign RC array to RC pointer (bad - fails "
  //       "compilation\n";
  //   }
  // RefCountedPtr<int> r8(new int);
  // r8 = r1;
  // if (verbose)
  //   {
  //     pout() << "\n-39- Assign non-const RC array to const RC pointer "
  //       "(bad - fails compilation\n";
  //   }
  // RefCountedPtr<const int> r9(new int);
  // r9 = r1;
  // if (verbose)
  //   {
  //     pout() << "\n-40- Assign const RC array to non-const RC pointer "
  //       "(bad - fails compilation\n";
  //   }
  // RefCountedPtr<int> r10(new int);
  // r10 = r2;

//--Test malloc and free

  int* mem = (int*) malloc(2*sizeof(int));
  mem[0] = -1;
  mem[1] = -2;
  if (verbose)
    {
      pout() << "\n-41- Malloc (ok)\n";
    }
  RefCountedPtr<int, RCPFreePolicy> m1(mem);
  if (m1[0] != -1 || m1[1] != -2) ++status;
  if (verbose)
    {
      pout() << "\n-42- Construct const mem from non-const (ok)\n";
    }
  RefCountedPtr<const int, RCPFreePolicy> m2(m1);
  // This would fail since m2.operator[]() returns type 'const int&'
  // m2[1] = 2;
  if (m2[0] != -1 || m2[1] != -2) ++status;
  // if (verbose)
  //   {
  //     pout() << "\n-43- Construct non-const mem from const (bad - fails "
  //       "compilation)\n";
  //   }
  // RefCountedPtr<int, RCPFreePolicy> m3(m2);
  // if (verbose)
  //   {
  //     pout() << "\n-44- Construct RC pointer from RC mem (bad - fails "
  //       "compilation)\n";
  //   }
  // RefCountedPtr<int> m4(m1);
  
  return status;
}

///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
  {
    if ( argv[i][0] == '-' ) //if it is an option
    {
      // compare 3 chars to differentiate -x from -xx
      if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
      {
        verbose = true ;
        // argv[i] = "" ;
      }
      else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
      {
        verbose = false ;
        // argv[i] = "" ;
      }
      else
      {
        break ;
      }
    }
  }
  return ;
}
