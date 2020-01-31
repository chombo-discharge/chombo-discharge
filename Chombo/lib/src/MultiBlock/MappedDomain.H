#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// bvs  06/06/06

#ifndef _MAPPEDDOMAIN_H_
#define _MAPPEDDOMAIN_H_


#include <iostream>
#include "MappedBlock.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "Box.H"
#include "Misc.H"
#include "SPACE.H"
#include <map>
#include <string>

// class MappedBlock;

#include "NamespaceHeader.H"

// ---------------------------------------------------------
class MappedDomain
{
public:

  // Constructors

  /// The default constructor.  The constructed as an empty domain.
  /**
   */
  MappedDomain ();

  /// Construct MappedDomain with a single \em a_domBox as computational domain
  /**
     This constructor defaults to non-periodic domain
  */
  MappedDomain(const Box& a_domBox,
               int a_defaultBC); // void* a_defaultBC


  /// useful constructor where you give it all it needs.
  /**
      constructor for multi-block domain.  During construction the input
    {\tt a\_blocks} are checked to see that mated surface specifications
     are complementary, that conformal mappings are mutual inverses,
      and that non-conformal matings map into one another.  Surfaces that
      are not explicitly set in the MappedBlock objects can be assigned a
      default boundary condition value.  If the default boudnary condition
      value is NULL, then MappedDomain construction aborts if an UNKOWN boundary
      does not mate with a known boundary.

  */
  MappedDomain(const Vector<MappedBlock>& a_blocks,
               int a_defaultBC = -1); // void* a_defaultBC = NULL


  /// Construct MappedDomain with \em a_domBox as computational domain
  /**
     This constructor defaults to non-periodic domain
  */
  void define(const Box& a_domBox);


  /// useful define for strong construction
  /**

  */
  void define(const Vector<MappedBlock>& a_blocks,
              int a_defaultBC = -1); // void* a_defaultBC = NULL

  //  Accessors


  /// Returns true if this MappedDomain is empty or undefined.
  /**
  */
  bool isEmpty () const;

  ///
  /** Returns number of blocks.
   */
  int numBlocks() const;

  ///  MappedBlock corresponding to BoxLayout::blockID call.
  const MappedBlock& block(int index) const
  {
    return m_blocks[index];
  }

  ///
  /** Returns whole domain
   */
  const Box& domain() const
  {
    return m_domain;
  }

  ///
  /** Returns index of block that contains given box.
   */
  int findBlock(const Box&   a_bx) const;

  ///
  /** Returns index of block that contains given point in physical space.
   */
  int findBlock(const RealVect&   a_physPoint) const;

  ///
  /** Returns index of block that contains valid cell that is
      a_cell in a_block's index space.
      If a_cell is in one of the domains of m_blocks, then that is the
      right block, so return it.
      Otherwise, consider a_cell as a ghost cell of a_block, and return
      the index of the neighboring block to a_block that contains it.
      If that also fails, then return NO_BLOCK.
   */
  int blockValid(const IntVect&   a_cell,
                 int              a_block) const;

  ///
  /**
     In a_ivs, return indices in block a_block of
     neighborhood of radius a_radius around cell a_baseCell.
     The new center of the neighborhood is a_centerCell,
     which is usually a_baseCell but may be shifted
     to get away from the external boundary.
  */
  void getNeighborhood(DenseIntVectSet&   a_ivs,
                       IntVect&           a_centerCell,
                       const IntVect&     a_baseCell,
                       int                a_block,
                       int                a_radius);

  ///
  /** Returns indices in a_blockNew corresponding to
      a_indicesOld in a_blockOld.
   */
  IntVect convertBetweenBlocks(int              a_blockOld,
                               int              a_blockNew,
                               const IntVect&   a_indicesOld);

  ///
  /** Returns a_indicesNew in a_blockNew corresponding to
      a_indicesOld in a_blockOld,
      but avoiding corner ghost cells of a_blockNew.
   */
  bool convertGhostBlocks(IntVect&         a_indicesNew,
                          int              a_blockNew,
                          const IntVect&   a_indicesOld,
                          int              a_blockOld);

  // refinement

  /// Refine this problem domain.
  /**
     Modifies this MappedDomain by refining it by given (positive) refinement
     ratio.  The Empty MappedDomain is not modified by this function.
  */
  MappedDomain& refine (int a_refinement_ratio);

  /// Return a MappedDomain which is a refinement of \em a_probdomain
  /**
     Returns a MappedDomain that is the argument MappedDomain refined by
     given (positive) refinement ratio.  The Empty MappedDomain is not
     modified by this function.
  */
  friend MappedDomain refine (const MappedDomain& a_probdomain,
                              int   a_refinement_ratio);

  /// Refine this MappedDomain
  /**
     Modifies this MappedDomain by refining it by given (positive) refinement
     ratio.  The Empty MappedDomain is not modified by this function.
  */
  MappedDomain& refine (const IntVect& a_refinement_ratio);

  /// Refinement function
  /**
     Returns a MappedDomain that is the argument MappedDomain refined
     by given (positive) refinement ratio.  The Empty MappedDomain is
     not modified by this function.

  */
  friend MappedDomain refine (const MappedDomain&     a_probdomain,
                              const IntVect& a_refinement_ratio);

  // coarsening

  /// Coarsen this MappedDomain
  /**
     Modifies this MappedDomain by coarsening it by given (positive)
     refinement ratio.  The Empty MappedDomain is not modified by
     this function.
  */

  MappedDomain& coarsen (int a_refinement_ratio);

  /// Coarsening function
  /**
     Returns a MappedDomain that is the argument MappedDomain coarsened
     by given (positive) refinement ratio.  The Empty MappedDomain is not
     modified by this function.
  */
  friend MappedDomain coarsen (const MappedDomain& a_probdomain,
                               int        a_refinement_ratio);

  /// Coarsen this MappedDomain
  /**
     Modifies this MappedDomain by coarsening by given (positive) refinement
     ratio.  The Empty MappedDomain is not modified by this function.
  */
  MappedDomain& coarsen (const IntVect& refinement_ratio);

  /// Coarsening function
  /**
     Returns a MappedDomain that is the argument MappedDomain coarsened
     by given (positive) refinement ratio.  The Empty MappedDomain is
     not modified by this function.

  */
  friend MappedDomain coarsen (const MappedDomain&  a_probdomain,
                               const IntVect& a_refinement_ratio);

  // I/O Functions

  /// Write an ASCII representation to the ostream.
  /**
  */
  void write(std::ostream&   os) const;

  /// Read from istream.
  /**
  */
  void define(std::istream& is);

  const static int NO_BLOCK = -1;

protected:

  friend class MappedDomainIO;

  Vector<MappedBlock> m_blocks;

  // smallest box containing all boxes of m_blocks
  Box m_domain;
};

#include "NamespaceFooter.H"

#endif // include guard
