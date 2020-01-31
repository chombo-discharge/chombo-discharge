#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//
// It's ok to #include Chombo headers here, but not CCSE headers (except for
// ccse facade headers like ccse/neutralFacade.H).
//

#include <iostream>

#include "CCSE_all/neutralFacade.H"
#include "IntVect.H" // Chombo header -- as are all others below here.
#include "Box.H"
#include "Vector.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"


Chombo::Box envelope( const std::vector<Chombo::Box>& boxes )
{
  Chombo::IntVect lo( boxes[0].smallEnd() ), hi( boxes[0].bigEnd() );
  for ( unsigned b=1;b<boxes.size();++b )
  {
    lo = boxes[b].smallEnd() < lo  ?  boxes[b].smallEnd() : lo;
    hi = boxes[b].bigEnd()   > hi  ?  boxes[b].bigEnd()   : hi;
  }
  return Chombo::Box( lo, hi );
}



int main(int argc, char** argv)
{
  BLfacade::initializeBoxLib(argc,argv);

  //
  // IntVect
  //
  Chombo::IntVect chIV( D_DECL(17, 31, 23) );
  std::cout << "chIV=" << chIV << '\n';

  BLfacade::IntVect blIV( chIV );
  std::cout << "blIV=" << blIV << '\n';

  Chombo::IntVect backagainIV( blIV );
  std::cout << "backagainIV=" << backagainIV << '\n';


  //
  // Box
  //
  Chombo::Box chB( Chombo::IntVect(D_DECL(1,2,3)),
                   Chombo::IntVect(D_DECL(11,12,13)) );
  std::cout << "chB=" << chB << '\n';

  BLfacade::Box blB( chB );
  std::cout << "blB=" << blB << '\n';

  Chombo::Box backagainB( blB );
  std::cout << "backagainB=" << backagainB << '\n';

  //
  // BoxArray
  //
  int const nboxes=3;
  Chombo::Vector<Chombo::Box> boxes(nboxes);
  for ( int ib=0; ib<nboxes; ++ib )
  {
    boxes[ib].define(Chombo::IntVect(D_DECL(1+ib*10,2,3)),
                     Chombo::IntVect(D_DECL(1+ib*10,12,13)));
  }
  Chombo::DisjointBoxLayout dbl;
  Chombo::Vector<int> procIDs;
  dbl.defineAndLoadBalance( boxes, &procIDs );

  BLfacade::BoxArray fBoxArray( dbl );
  std::cout << fBoxArray << '\n';

  //
  // MultiFab
  //
  int const ncomps = 2;
  int const nghosts = 1;

  Chombo::Box boundingBox( envelope( boxes.stdVector() ) );

  BLfacade::MultiFab multifab1( fBoxArray, ncomps, nghosts, 0 );
  Chombo::LevelData<Chombo::FArrayBox>* pLD1 =
    multifab1.pLevelData( boundingBox, Chombo::IntVect::Unit );
  CH_assert( pLD1 );

  for ( Chombo::DataIterator dit(pLD1->dataIterator()); dit.ok(); ++dit )
  {
    (*pLD1)[dit()].setVal(1.1);
    CH_assert( (*pLD1)[dit()].dataPtr(0)[2] == 1.1 );
  }


  BLfacade::BoxArray fBoxArray2( dbl );
  BLfacade::MultiFab multifab2( fBoxArray2, ncomps, nghosts,
                                &(procIDs.stdVector()) );
  Chombo::LevelData<Chombo::FArrayBox>* pLD2 =
    multifab2.pLevelData( boundingBox, Chombo::IntVect::Unit );

  CH_assert( pLD2 );

  for ( Chombo::DataIterator dit(pLD2->dataIterator()); dit.ok(); ++dit )
  {
    (*pLD2)[dit()].setVal(2.2);
    CH_assert( (*pLD2)[dit()].dataPtr(0)[2] == 2.2 );
  }

  // Construct another MultiFab just like multifab2, and check that their DBLs
  // are equivalent.
  BLfacade::MultiFab multifab3( fBoxArray2, ncomps, nghosts,
                                &(fBoxArray2.dbl()->procIDs().stdVector()) );
  Chombo::LevelData<Chombo::FArrayBox>* pLD3 =
    multifab3.pLevelData( boundingBox, Chombo::IntVect::Unit );
  CH_assert( pLD2->disjointBoxLayout() == pLD3->disjointBoxLayout() );

  BLfacade::finalizeBoxLib();
}
