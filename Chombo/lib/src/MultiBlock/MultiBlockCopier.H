#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKCOPIER_H_
#define _MULTIBLOCKCOPIER_H_

#include "DisjointBoxLayout.H"
#include "Pool.H"
#include "Vector.H"
#include "ProblemDomain.H"
#include "Copier.H"
#include "NamespaceHeader.H"

/// An even stranger (than Copier) thing to copy from valid cells of one block to ghost cells of another block.
/**
   This is a specialized Copier to fill ghost cells of one block using
   valid cells of another block.

   Each "block" is specified by a box.  A box in the layout is said to
   be in the block if it is contained in the block box.

   This is the same as a regular Copier if the source block box is the
   entire source domain, and the destination block box is the entire
   destination domain.
*/

class MultiBlockCopier: public Copier
{
public:

  ///null constructor, copy constructor and operator= can be compiler defined.
  MultiBlockCopier()
  {
    // Initialize as empty boxes.
    m_levelBlockDomain = Box();
    m_destBlockDomain = Box();
  }

  /// Basic defining constructor
  /** Note that this constructor doesn't take domains of
      source and dest blocks, so it's useless for MultiBlockCopier.
      Provided for compatibility with Copier base class.
  */
  MultiBlockCopier(const DisjointBoxLayout& a_level,
                   const BoxLayout& a_dest,
                   bool a_exchange = false);

  /// another defining constructor
  /** Note that this constructor doesn't take domains of
      source and dest blocks, so it's useless for MultiBlockCopier.
      Provided for compatibility with Copier base class.
  */
  MultiBlockCopier(const DisjointBoxLayout& a_level,
                   const BoxLayout& a_dest,
                   const ProblemDomain& a_domain,
                   bool  a_exchange  = false);

  /// another defining constructor
  /** Note that this constructor doesn't take domains of
      source and dest blocks, so it's useless for MultiBlockCopier.
      Provided for compatibility with Copier base class.
  */
  MultiBlockCopier(const DisjointBoxLayout& a_level,
                   const BoxLayout& a_dest,
                   const IntVect& a_ghost,
                   bool  a_exchange = false);

  /// another defining constructor
  /** Note that this constructor doesn't take domains of
      source and dest blocks, so it's useless for MultiBlockCopier.
      Provided for compatibility with Copier base class.
  */
  MultiBlockCopier(const DisjointBoxLayout& a_level,
                   const BoxLayout& a_dest,
                   const ProblemDomain& a_domain,
                   const IntVect& a_ghost,
                   bool  a_exchange = false);

  /// defining constructor with domains of source and dest blocks
  /** Constructor taking domains of source and destination blocks.
  */
  MultiBlockCopier(const DisjointBoxLayout& a_level,
                   const BoxLayout& a_dest,
                   const IntVect& a_ghost,
                   const Box& a_levelBlockDomain,
                   const Box& a_destBlockDomain,
                   bool  a_exchange = false);

  /// defining constructor with domains of source and dest blocks
  /** Constructor taking domains of source and destination blocks.
  */
  MultiBlockCopier(const DisjointBoxLayout& a_level,
                   const BoxLayout& a_dest,
                   const ProblemDomain& a_domain,
                   const IntVect& a_ghost,
                   const Box& a_levelBlockDomain,
                   const Box& a_destBlockDomain,
                   bool  a_exchange = false);

  ///
  virtual ~MultiBlockCopier();

  /// define to match constructor (also doesn't do anything)
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      bool  a_exchange = false);

  /// define to match constructor (also doesn't do anything)
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      bool  a_exchange = false);

  /// define to match constructor (also doesn't do anything)
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const IntVect& a_ghost,
                      bool  a_exchange = false);

  /// define to match original (also doesn't do anything)
  virtual void define(const BoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      const IntVect& a_ghost,
                      bool  a_exchange = false);

  /// define to match constructor that actually does something
  virtual void define(const DisjointBoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const IntVect& a_ghost,
                      const Box& a_levelBlockDomain,
                      const Box& a_destBlockDomain,
                      bool  a_exchange = false);

  /// define to match constructor that actually does something
  virtual void define(const BoxLayout& a_level,
                      const BoxLayout& a_dest,
                      const ProblemDomain& a_domain,
                      const IntVect& a_ghost,
                      const Box& a_levelBlockDomain,
                      const Box& a_destBlockDomain,
                      bool  a_exchange = false);

  ///
  virtual void clear();

  const Box& levelBlockDomain()
  {
    return m_levelBlockDomain;
  }

  const Box& destBlockDomain()
  {
    return m_destBlockDomain;
  }

protected:

private:

  /// domain of source block
  Box m_levelBlockDomain;

  /// domain of destination block
  Box m_destBlockDomain;
};

std::ostream& operator<<(std::ostream& os, const MultiBlockCopier& copier);


#include "NamespaceFooter.H"
#endif
