#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Tues, Oct 5 1999

#include <cmath>
#include "BoxLayout.H"
#include "Misc.H"
#include "DebugDump.H"
#include "UsingNamespace.H"

/**
   mortonTest returns:
    0: all tests passed
 */
int
mortonTest()
{
  int retval = 0;
  //make two sets of the same boxes but in different order
  Vector<Box> boxesOne(16);
  Vector<Box> boxesTwo(16);
  boxesOne[0 ] = Box(IntVect(D_DECL6(0  , 0 ,0,0,0,0)) , IntVect(D_DECL6(31 ,31 ,0,0,0,0)));      
  boxesOne[1 ] = Box(IntVect(D_DECL6(0  ,32 ,0,0,0,0)) , IntVect(D_DECL6(31 ,63 ,0,0,0,0)));      
  boxesOne[2 ] = Box(IntVect(D_DECL6(0  ,64 ,0,0,0,0)) , IntVect(D_DECL6(31 ,95 ,0,0,0,0)));      
  boxesOne[3 ] = Box(IntVect(D_DECL6(0  ,96 ,0,0,0,0)) , IntVect(D_DECL6(31 ,127,0,0,0,0)));
  boxesOne[4 ] = Box(IntVect(D_DECL6(32 , 0 ,0,0,0,0)) , IntVect(D_DECL6(63 ,31 ,0,0,0,0)));    
  boxesOne[5 ] = Box(IntVect(D_DECL6(32 ,32 ,0,0,0,0)) , IntVect(D_DECL6(63 ,63 ,0,0,0,0)));    
  boxesOne[6 ] = Box(IntVect(D_DECL6(32 ,64 ,0,0,0,0)) , IntVect(D_DECL6(63 ,95 ,0,0,0,0)));    
  boxesOne[7 ] = Box(IntVect(D_DECL6(32 ,96 ,0,0,0,0)) , IntVect(D_DECL6(63 ,127,0,0,0,0)));
  boxesOne[8 ] = Box(IntVect(D_DECL6(64 ,0  ,0,0,0,0)) , IntVect(D_DECL6(95 ,31 ,0,0,0,0)));    
  boxesOne[9 ] = Box(IntVect(D_DECL6(64 ,32 ,0,0,0,0)) , IntVect(D_DECL6(95 ,63 ,0,0,0,0)));    
  boxesOne[10] = Box(IntVect(D_DECL6(64 ,64 ,0,0,0,0)) , IntVect(D_DECL6(95 ,95 ,0,0,0,0)));    
  boxesOne[11] = Box(IntVect(D_DECL6(64 ,96 ,0,0,0,0)) , IntVect(D_DECL6(95 ,127,0,0,0,0)));
  boxesOne[12] = Box(IntVect(D_DECL6(96 ,0  ,0,0,0,0)) , IntVect(D_DECL6(127,31 ,0,0,0,0)));  
  boxesOne[13] = Box(IntVect(D_DECL6(96 ,32 ,0,0,0,0)) , IntVect(D_DECL6(127,63 ,0,0,0,0)));  
  boxesOne[14] = Box(IntVect(D_DECL6(96 ,64 ,0,0,0,0)) , IntVect(D_DECL6(127,95 ,0,0,0,0)));  
  boxesOne[15] = Box(IntVect(D_DECL6(96, 96 ,0,0,0,0)) , IntVect(D_DECL6(127,127,0,0,0,0)));
  boxesTwo[0 ] = Box(IntVect(D_DECL6(0  , 0 ,0,0,0,0)) , IntVect(D_DECL6(31 ,31 ,0,0,0,0)));      
  boxesTwo[14] = Box(IntVect(D_DECL6(0  ,32 ,0,0,0,0)) , IntVect(D_DECL6(31 ,63 ,0,0,0,0)));      
  boxesTwo[3 ] = Box(IntVect(D_DECL6(0  ,64 ,0,0,0,0)) , IntVect(D_DECL6(31 ,95 ,0,0,0,0)));      
  boxesTwo[2 ] = Box(IntVect(D_DECL6(0  ,96 ,0,0,0,0)) , IntVect(D_DECL6(31 ,127,0,0,0,0)));
  boxesTwo[5 ] = Box(IntVect(D_DECL6(32 , 0 ,0,0,0,0)) , IntVect(D_DECL6(63 ,31 ,0,0,0,0)));    
  boxesTwo[4 ] = Box(IntVect(D_DECL6(32 ,32 ,0,0,0,0)) , IntVect(D_DECL6(63 ,63 ,0,0,0,0)));    
  boxesTwo[7 ] = Box(IntVect(D_DECL6(32 ,64 ,0,0,0,0)) , IntVect(D_DECL6(63 ,95 ,0,0,0,0)));    
  boxesTwo[6 ] = Box(IntVect(D_DECL6(32 ,96 ,0,0,0,0)) , IntVect(D_DECL6(63 ,127,0,0,0,0)));
  boxesTwo[9 ] = Box(IntVect(D_DECL6(64 ,0  ,0,0,0,0)) , IntVect(D_DECL6(95 ,31 ,0,0,0,0)));    
  boxesTwo[8 ] = Box(IntVect(D_DECL6(64 ,32 ,0,0,0,0)) , IntVect(D_DECL6(95 ,63 ,0,0,0,0)));    
  boxesTwo[11] = Box(IntVect(D_DECL6(64 ,64 ,0,0,0,0)) , IntVect(D_DECL6(95 ,95 ,0,0,0,0)));    
  boxesTwo[10] = Box(IntVect(D_DECL6(64 ,96 ,0,0,0,0)) , IntVect(D_DECL6(95 ,127,0,0,0,0)));
  boxesTwo[13] = Box(IntVect(D_DECL6(96 ,0  ,0,0,0,0)) , IntVect(D_DECL6(127,31 ,0,0,0,0)));  
  boxesTwo[12] = Box(IntVect(D_DECL6(96 ,32 ,0,0,0,0)) , IntVect(D_DECL6(127,63 ,0,0,0,0)));  
  boxesTwo[15] = Box(IntVect(D_DECL6(96 ,64 ,0,0,0,0)) , IntVect(D_DECL6(127,95 ,0,0,0,0)));  
  boxesTwo[1 ] = Box(IntVect(D_DECL6(96, 96 ,0,0,0,0)) , IntVect(D_DECL6(127,127,0,0,0,0)));

  mortonOrdering(boxesOne);
  mortonOrdering(boxesTwo);
  for(int ibox = 0; ibox < 16; ibox++)
    {
      if(boxesOne[ibox] != boxesTwo[ibox])
        {
          pout() << "mismatch at ibox = " << ibox << ", boxesone[ibox] = " << boxesOne[ibox] << ", boxestwo[ibox] = " << boxesTwo[ibox] << endl;
          retval = ibox+1;
        }
    }
  pout() <<"post ordering vbox one"  << endl;
  dumpVBox(&boxesOne);
  pout() <<"post ordering vbox two"  << endl;
  dumpVBox(&boxesTwo);
  return retval;
}

/// Code:
int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  int icode = mortonTest();
  if (icode != 0)
    {
      pout() << " mortonOrdering failed with error code " << icode   << endl;
    }
  else
    {
      pout() << "mortonOrdering passed all tests"   << endl;
    }
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return icode;
}


