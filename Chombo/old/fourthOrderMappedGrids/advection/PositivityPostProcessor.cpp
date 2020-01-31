#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PositivityPostProcessor.H"

#include "Box.H"
#include "BoxIterator.H"
#include "BoxLayoutData.H"
#include "Copier.H"
#include "DataIterator.H"
#include "DisjointBoxLayout.H"
#include "Interval.H"
#include "PositivityPostProcessorF_F.H"

#include "NamespaceHeader.H"


PositivityPostProcessor::PositivityPostProcessor(
   const IntVect& a_halo,
   const int&     a_maximum_number_of_passes,
   const bool&    a_verbose)
{
   define( a_halo, a_maximum_number_of_passes, a_verbose);
}


void PositivityPostProcessor::define(
   const IntVect& a_halo,
   const int&     a_maximum_number_of_passes,
   const bool&    a_verbose)
{
   CH_assert( a_halo>=IntVect::Unit );
   m_neighborhood.define( -a_halo, a_halo );
   CH_assert( a_maximum_number_of_passes>0 );
   m_maximum_number_of_passes = a_maximum_number_of_passes;
   m_verbose = a_verbose;
}


inline bool notDone( const int& a_count )
{
   return ( a_count>0 ? true : false );
}


void PositivityPostProcessor::enforce( LevelData<FArrayBox>& a_phi,
                                       const Real& a_ref_val ) const
{
   checkIfDefined();
   CH_assert( a_ref_val>0.0 );

   const IntVect& nghosts( a_phi.ghostVect() );
   checkNeighborhood( nghosts );

   const DisjointBoxLayout& boxes( a_phi.getBoxes() );
   Copier copier;
   copier.define( boxes, boxes, nghosts, true );
   copier.reverse();

   LevelData<FArrayBox> deltaPhi( boxes, a_phi.nComp(), nghosts );

   int count( numberOfNegativeValues( a_phi ) );
   int pass_number(0);
   while ( notDone( count ) && (pass_number<m_maximum_number_of_passes) )
   {
      if (m_verbose)
      {
         pout() << "  Minimum Value Redistribution Pass " << pass_number << endl;
         pout() << "    Number of values less than minVal:\t" << count << std::endl;
         if (procID()==0)
         {
            cout << "  Minimum Value Redistribution Pass " << pass_number << endl;
            cout << "    Number of values less than minVal:\t" << count << std::endl;
         }
      }

      computeRedistribution( deltaPhi, a_phi, a_ref_val );
      accumulateRedistribution( deltaPhi, copier );
      applyRedistribution( a_phi, deltaPhi );
      count = numberOfNegativeValues( a_phi );
      pass_number++;

   } // end loop over passes

   if (m_verbose)
   {
      pout() << "    Number of values less than minVal:\t" << count << std::endl;
      if (procID()==0)
      {
         cout << "    Number of values less than minVal:\t" << count << std::endl;
      }
   }

   if (notDone( count ))
   {
      if (procID()==0)
         MayDay::Warning("PositivityPostProcessor: Minimum value not enforced on entire domain!");
   }
}


void
PositivityPostProcessor::computeRedistribution(
   LevelData<FArrayBox>& a_deltaPhi,
   const LevelData<FArrayBox>& a_phi,
   const Real& a_ref_val ) const
{
   const DisjointBoxLayout& grids( a_phi.getBoxes() );

   DataIterator dit( a_phi.dataIterator() );
   int unable_to_redistribute(0);
   for (dit.begin(); dit.ok(); ++dit)
   {

      Box gridBox( grids[dit] );
      const FArrayBox& thisPhi( a_phi[dit] );
      FArrayBox& thisDeltaPhi( a_deltaPhi[dit] );
      thisDeltaPhi.copy( thisPhi );

      // redistribute on interior of box
      FORT_REDISTRIBUTENEGATIVES(CHF_FRA(thisDeltaPhi),
                                 CHF_BOX(gridBox),
                                 CHF_BOX(m_neighborhood),
                                 CHF_CONST_REAL(a_ref_val),
                                 CHF_INT(unable_to_redistribute));

      // save to fix boundaries later
      thisDeltaPhi -= thisPhi;
   } // end loop over grid boxes

   int count(unable_to_redistribute);
#ifdef CH_MPI
   MPI_Allreduce( &unable_to_redistribute, &count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
#endif
   if (m_verbose)
   {
      pout() << "    Number of values unable to redistribute:\t" << count << std::endl;
      if (procID()==0)
         cout << "    Number of values unable to redistribute:\t" << count << std::endl;
   }
}


class AddOp : public LDOperator<FArrayBox>
{
   public:
      Real scale;
      AddOp()
        :
        scale(1.0)
      {
      }

      virtual void linearIn( FArrayBox& arg,
                             void* buf,
                             const Box& R,
                             const Interval& comps ) const
      {
         Real* buffer = (Real*)buf;
         if (scale != 1.0)
         {
            ForAllXBNNnoindx( Real, arg, R, comps.begin(), comps.size() )
            {
               argR += (*buffer) * scale;
               buffer++;
            } EndFor
                 }
         else
         {
            ForAllXBNNnoindx( Real, arg, R, comps.begin(), comps.size() )
            {
               argR += (*buffer);
               buffer++;
            } EndFor
                 }
      }

      void op( FArrayBox& dest,
               const Box& RegionFrom,
               const Interval& Cdest,
               const Box& RegionTo,
               const FArrayBox& src,
               const Interval& Csrc ) const
      {
         if (scale != 1.0)
            dest.plus(src, RegionFrom, RegionTo, scale, Csrc.begin(), Cdest.begin(), Cdest.size());
         else
            dest.plus(src, RegionFrom, RegionTo, Csrc.begin(), Cdest.begin(), Cdest.size());
      }
};


void PositivityPostProcessor::accumulateRedistribution(
   LevelData<FArrayBox>& a_deltaPhi,
   const Copier& a_reverseCopier ) const
{
   AddOp accumOp;
   a_deltaPhi.copyTo( a_deltaPhi, a_reverseCopier, accumOp );
}


void PositivityPostProcessor::applyRedistribution(
   LevelData<FArrayBox>& a_phi,
   const LevelData<FArrayBox>& a_deltaPhi ) const
{
   const DisjointBoxLayout& grids( a_phi.getBoxes() );
   DataIterator dit( a_phi.dataIterator() );
   for (dit.begin(); dit.ok(); ++dit)
   {
      Box gridBox( grids[dit] );
      FArrayBox& thisPhi( a_phi[dit] );
      const FArrayBox& thisDeltaPhi( a_deltaPhi[dit] );
      thisPhi += thisDeltaPhi;
   } // end loop over grid boxes
   a_phi.exchange();
}


int PositivityPostProcessor::numberOfNegativeValues(
   const LevelData<FArrayBox>& a_phi ) const
{
   const DisjointBoxLayout& grids( a_phi.getBoxes() );
   DataIterator dit( a_phi.dataIterator() );

   int count_lcl(0);
   for (dit.begin(); dit.ok(); ++dit)
   {
      Box gridBox( grids[dit] );
      const FArrayBox& thisPhi( a_phi[dit] );
      if (m_verbose)
      {
         FORT_FINDALLNEGATIVES(CHF_INT(count_lcl),
                               CHF_CONST_FRA(thisPhi),
                               CHF_BOX(gridBox));
      }
      else
      {
         FORT_FINDANYNEGATIVES(CHF_INT(count_lcl),
                               CHF_CONST_FRA(thisPhi),
                               CHF_BOX(gridBox));
      }
   } // end loop over grid boxes

   int count(count_lcl);
#ifdef CH_MPI
   MPI_Allreduce( &count_lcl, &count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
#endif

   return count;
}


void PositivityPostProcessor::checkNeighborhood(
   const IntVect& a_nghosts ) const
{
   if (m_neighborhood.bigEnd()>a_nghosts)
   {
      if (procID()==0)
         MayDay::Error("PositivityPostProcessor: Redistribution neighborhood is larger than the number of avialable ghost cells!");
   }
}


void PositivityPostProcessor::checkIfDefined() const
{
   if (m_maximum_number_of_passes<0)
   {
      if (procID()==0)
         MayDay::Error("PositivityPostProcessor: Attempt to use undefined object!");
   }
}


#include "NamespaceFooter.H"
