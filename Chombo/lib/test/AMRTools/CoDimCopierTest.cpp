#include "LevelData.H"
#include "LoadBalance.H"
#include "CornerCopier.H"
#include "CoDimCopier.H"
#include "CoDimCopyManager.H"
#include "FABView.H"
#include "parstream.H"

//#include "CoDimCopyManager.H"

#ifdef CH_USE_PETSC
#include "petsc.h"
#endif 

#include "UsingNamespace.H"

int n;
int verbosity;
Real epsilon = 1.0e-8;


int checkValues(const LevelData<FArrayBox>& ld)
{
  int status = 0;
  const DisjointBoxLayout& dbl = ld.getBoxes();
  // now loop over all ghost cells and verify that they are what they should be
  DataIterator dit = dbl.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FArrayBox& thisData = ld[dit];
      const Box& thisBox = dbl[dit];
      const IntVect& ghostVect = ld.ghostVect();
      for (int faceDir=0; faceDir<SpaceDim; faceDir++)
        {
          // lo-side
          Box loGhost = adjCellLo(thisBox,faceDir,1);
          loGhost.grow(ghostVect);
          loGhost.grow(faceDir, -1*ghostVect[faceDir]);
          // hi-side
          Box hiGhost = adjCellHi(thisBox,faceDir,1);
          hiGhost.grow(ghostVect);
          hiGhost.grow(faceDir, -1*ghostVect[faceDir]);           
          
          BoxIterator loBit(loGhost);
          for (loBit.begin(); loBit.ok(); ++loBit)
            {
              IntVect iv = loBit();
              int sumIV = iv.sum();
              // account for periodicity
              if (iv[SpaceDim-1] == -1)
                {
                  sumIV += 2*n;
                }
              if (iv[SpaceDim-1] == (2*n))
                {
                  sumIV -= 2*n;
                }
              
              if (abs(thisData(iv,0)-sumIV) > epsilon)
                {
                  // fail!
                  if (verbosity > 1)
                    {
                      pout() << "CoDimCopierTest Failed: Value at "
                             << iv << " should be " << sumIV << " but is actually "
                                  << thisData(iv,0)
                             << " lo side, faceDir = " << faceDir << endl;
                    }
                  ++status;
                }
            }
          
          // could do this with a SideIterator, but I'm being lazy
          BoxIterator hiBit(hiGhost);
          for (hiBit.begin(); hiBit.ok(); ++hiBit)
            {
              IntVect iv = hiBit();
              int sumIV = iv.sum();
              // account for periodicity
              if (iv[SpaceDim-1] == -1)
                {
                  sumIV += 2*n;
                }
              if (iv[SpaceDim-1] == (2*n))
                {
                  sumIV -= 2*n;
                }
              
              if (abs(thisData(iv,0)-sumIV) > epsilon)
                {
                  // fail!
                  if (verbosity > 1)
                    {                       
                      pout() << "CoDimCopierTest Failed: Value at "
                             << iv << " should be " << sumIV << " but is "
                             << thisData(iv,0)
                             << " hi side, faceDir = " << faceDir << endl;
                    }
                  ++status;
                }
            }
        } // end loop over face dir
    } // end loop over boxes           
  return status;
}


void initData(LevelData<FArrayBox>& ld)
{
  const DisjointBoxLayout& dbl = ld.getBoxes();
  
  DataIterator dit = dbl.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& gridBox = dbl[dit];
      FArrayBox& thisData = ld[dit];
      
      thisData.setVal(10000000.0);
      // set the interior (valid-region) values to zero
      thisData.setVal(0.0, dbl[dit],0);
      
      // now loop over face directions and initialize coDim1 ghost cells, which are the starting
      // point for any CoDim copy chains.
      for (int dir=0; dir<SpaceDim; dir++)
        {
          Box coDim1Box;
          SideIterator sit;
          for (sit.begin(); sit.ok(); ++sit)
            {
              if (sit() == Side::Lo)
                {
                  coDim1Box = adjCellLo(gridBox,dir,1);
                }
              else
                {
                  coDim1Box = adjCellHi(gridBox,dir,1);
                }
              
              // need to grow in directions transverse to the box-coupled direction
              if (dir != (SpaceDim-1) )
                {
                  coDim1Box.grow(ld.ghostVect());
                  coDim1Box.grow(dir, -1);
                  coDim1Box.grow((SpaceDim-1),-1);
                }
              
              // Set the values in box1_codim2 to 1 and
              // the values in box2_codim2 to 2.
              BoxIterator bit(coDim1Box);
              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  int sumIV = iv.sum();
                  // account for periodicity
                  if (iv[SpaceDim-1] == -1)
                    {
                      sumIV += 2*n;
                    }
                  if (iv[SpaceDim-1] == (2*n))
                    {
                      sumIV -= 2*n;
                    }
                  
                  thisData(iv,0) = sumIV;
                } // end loop over coDim1 set
            } //end loop over high-lo
        } // end loop over face directions
    } // end loop over boxes  
}


int main( int argc, char* argv[] )
{
#ifdef CH_USE_PETSC
  int ierr;
  ierr = PetscInitialize(&argc, &argv,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
#else
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif 
#endif // end petsc conditional

   int status = 0;   
   verbosity = 4;

   // this doesn't make any sense in 1D
   if (SpaceDim > 1)
     {
       // Make two 4x4x4 boxes that are adjacent in the SpaceDim-1 direction
       
       n = 4;
       
       Box box1(IntVect::Zero,(n-1)*IntVect::Unit);
       Box box2(box1);
       box2.shift(SpaceDim-1,n);
       
       // Make a LevelData on these boxes with one ghost cell
       
       Vector<Box> boxes;
       boxes.push_back(box1);
       boxes.push_back(box2);
       
       Vector<int> pids;
       LoadBalance(pids, boxes);
       
       Box domainBox(boxes[0]);
       domainBox.growHi(SpaceDim-1, n);   
       ProblemDomain probDomain(domainBox);
       probDomain.setPeriodic(SpaceDim-1,true);
       
       DisjointBoxLayout dbl(boxes, pids, probDomain);
       
       
       LevelData<FArrayBox> ld(dbl, 1, IntVect::Unit);

       initData(ld);
       
       // step through to enable viz
       DataIterator dit = dbl.dataIterator();
       for (dit.begin(); dit.ok(); ++dit)
         {
           FArrayBox& thisData = ld[dit];
           const Box& thisBox = dbl[dit];
         }
       
       // Use CoDim copiers to make the ghost cells consistent
     
       if (verbosity > 4) pout() << "manually organizing coDim copies..." << endl;
       
       //   ld.exchange();  // Just for the heck of it
       //CornerCopier copier(dbl, dbl, dbl.physDomain(), IntVect::Unit, true);
       for (int coDim=2; coDim<=SpaceDim; coDim++)
         {
           CoDimCopier coDimCopy(dbl, dbl, dbl.physDomain(), IntVect::Unit, coDim, true);
           
           ld.exchange(coDimCopy);
           
           // step through to enable viz
           for (dit.begin(); dit.ok(); ++dit)
             {
               FArrayBox& thisData = ld[dit];
               const Box& thisBox = dbl[dit];
             }
         }
       
       int checkVal = checkValues(ld);
       status += checkVal;

       // now do this again using the CoDimCopyManager...
       LevelData<FArrayBox> ld2(dbl, 1, IntVect::Unit);

       initData(ld2);
       
       // step through to enable viz
          for (dit.begin(); dit.ok(); ++dit)
         {
           FArrayBox& thisData = ld2[dit];
           const Box& thisBox = dbl[dit];
         }
       
       // Use CoDimCopyManager to orchestrate it all
       if (verbosity > 4) pout() << "now using coDimCopyManager..." << endl;
       
       CoDimCopyManager<FArrayBox> manager(dbl, dbl, ld2.ghostVect(), true);
     
       manager.manageExchanges(ld2);
       
       // step through to enable viz
       for (dit.begin(); dit.ok(); ++dit)
         {
           FArrayBox& thisData = ld2[dit];
           const Box& thisBox = dbl[dit];
         }
     
       
       checkVal = checkValues(ld2);
       status += checkVal;
       
     }  // end if DIM > 1
   
   if (verbosity > 1)
     {
       if (status == 0)
         {
           pout() << "CoDimCopierTest PASSED all tests" << endl;
         }
       else
         {
           pout() << "CoDimCopierTest FAILED at least one test" << endl; 
           pout() << "Status = " << status << endl;
         }      
     }
   
#ifdef CH_USE_PETSC
   ierr = PetscFinalize(); CHKERRQ(ierr);
#else
#ifdef CH_MPI
   MPI_Finalize();
#endif // mpi conditional
#endif // petsc conditional  
   
   
   return status;
}

