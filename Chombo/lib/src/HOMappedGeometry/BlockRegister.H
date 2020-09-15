#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// jafh, 5 Jan 2011

#ifndef _BLOCKREGISTER_H_
#define _BLOCKREGISTER_H_

#include <map>

#include "MultiBlockCoordSys.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LoHiSide.H"
#include "BlockBaseRegister.H"

#include "NamespaceHeader.H"

/// Class to hold data on a mapped multiblock domain block boundary.
class BlockRegister
{
public:

  /**
     \name Constructor, destructor, and defines.
  */
  /*@{*/

  /// Empty constructor.
  /**
   */
   BlockRegister()
     :
     m_verbosity(0),
     m_isOpen(true),
     m_fluxRegister(),
     m_auxRegister()
   {}

  /// Constructor.
  /**
     – a_mblock: Pointer to multiblock coordinate system.
     – a_grids: DisjointBoxLayout over which multiblock grid is defined.
     – a_ghost: Extension of register domain for auxiliary data.
   */
  BlockRegister(const MultiBlockCoordSys*  a_mblock,
                /// CELL-centered grids
                const DisjointBoxLayout&   a_grids,
                /// ghost width for auxiliary data
                const int&                 a_ghost = 0 );

  /// virtual destructor
  virtual ~BlockRegister();

  /// Define function that matches constructor
  /**
     – a_mblock: Pointer to multiblock coordinate system.
     – a_grids: DisjointBoxLayout over which multiblock grid is defined.
     – a_ghost: Extension of register domain for auxiliary data.
  */
  void define(const MultiBlockCoordSys*  a_mblock,
              /// CELL-centered grids
              const DisjointBoxLayout&   a_grids,
              /// ghost width for auxiliary data
              const int&                 a_ghost = 0 );

  /// Set verbosity.
  /**
     Default is 0.  Positive values output increasing amounts of information.
   */
  virtual void setVerbose(int   a_verbose);

  /// Store fluxes on block-boundary faces.
  /**
     – a_flux: Face-centered flux to be stored.
     – a_dataIndex: Index of box where flux has been computed.
     – a_dir: Direction of the face.
     – a_side: Low or high side of the box where the block-boundary
     face lies.
  */
  void storeFlux(const FArrayBox& a_flux,
                 const DataIndex& a_dataIndex,
                 int a_dir,
                 Side::LoHiSide a_side);
 
  /// increment fluxes on block-boundary faces.
  /**
     – a_flux: Face-centered flux to be stored.
     – a_dataIndex: Index of box where flux has been computed.
     – a_dir: Direction of the face.
     – a_side: Low or high side of the box where the block-boundary
     face lies.
  */
  void incrementFlux(const FArrayBox& a_flux,
                     const DataIndex& a_dataIndex,
                     int a_dir,
                     Side::LoHiSide a_side);

  /// Store auxiliary data on block-boundary faces including layers of ghost faces.
  /**
     – a_data: Face-centered data to be stored.
     – a_dataIndex: Index of box containing data to be stored.
     – a_dir: Direction of the face.
     – a_side: Low or high side of the box where the block-boundary
     face lies.
  */
  void storeAux(const FArrayBox& a_data,
                const DataIndex& a_dataIndex,
                int a_dir,
                Side::LoHiSide a_side);

/// Store auxiliary data on block-boundary faces including layers of ghost faces.
  /**
     – a_data: Face-centered data to be stored.
     – a_dataIndex: Index of box containing data to be stored.
     – a_dir: Direction of the face.
     – a_side: Low or high side of the box where the block-boundary
     face lies.
  */
  void incrementAux(const FArrayBox& a_data,
                    const DataIndex& a_dataIndex,
                    int a_dir,
                    Side::LoHiSide a_side);

  ///
  /**
     Perform data manipulation of flux and auxiliary data, after it has been
     stored, so that accessor functions have access to it.
  */
  void close();

  ///
  /**
     Sets the contents of the register to zero, and sets the register to be
     open.  If data is not already defined, argument a_n_comp must be defined
     to be greater than zero. 
     – a_n_comp: Number of components in register.
  */
  void setToZero( const int a_n_comp=-1 );

  /// Returns whether or not a box has a block-boundary interface on one of its faces.
  /**
     – a_dataIndex: Index of box.
     – a_dir: Direction of the face of the box.
     – a_side: Low or high side of the box.
  */
  bool hasInterface(const DataIndex& a_dataIndex,
                    int a_dir,
                    Side::LoHiSide a_side) const;

  /// Return fluxes on a block-boundary interface of a box, if hasInterface(a_dataIndex, a_dir, a_side) is true.
  /**
     – a_flux: Face-centered flux returned.
     – a_dataIndex: Index of box.
     – a_dir: Direction of the face.
     – a_side: Low or high side of the box.
     – a_sideData: Low or high side of the block-boundary interface,
     from which flux will be taken.
  */
  void getFlux(FArrayBox& a_flux,
               const DataIndex& a_dataIndex,
               int a_dir,
               Side::LoHiSide a_side,
               Side::LoHiSide a_sideData) const;

  /// Return auxiliary data on a block-boundary interface of a box, if hasInterface(a_dataIndex, a_dir, a_side) is true.
  /**
     – a_flux: Face-centered data returned.
     – a_dataIndex: Index of box.
     – a_dir: Direction of the face.
     – a_side: Low or high side of the box.
     – a_sideData: Low or high side of the block-boundary interface,
     from which data will be taken.
  */
  void getAux(FArrayBox& a_data,
              const DataIndex& a_dataIndex,
              int a_dir,
              Side::LoHiSide a_side,
              Side::LoHiSide a_sideData) const ;

protected:

  /// Returns true if register is open, false otherwise.
  inline bool isNotClosed() const
  {
    return m_isOpen;
  }

  /// Returns false if register is open, true otherwise.
  inline bool isClosed() const
  {
    return !m_isOpen;
  }

  int m_verbosity;
  bool m_isOpen;

  BlockBaseRegister<FArrayBox> m_fluxRegister;
  BlockBaseRegister<FArrayBox> m_auxRegister;
};

#include "NamespaceFooter.H"

#endif // include guard
