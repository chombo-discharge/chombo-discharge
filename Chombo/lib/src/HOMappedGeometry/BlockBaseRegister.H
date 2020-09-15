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

#ifndef _BLOCKBASEREGISTER_H_
#define _BLOCKBASEREGISTER_H_

#include <map>

#include "MultiBlockCoordSys.H"
#include "BaseFab.H"
#include "LevelData.H"
#include "LoHiSide.H"

#include "NamespaceHeader.H"

/// Class to implement a single data register a mapped multiblock domain block boundary.

template <class T> class BlockBaseRegister
{
public:

  /**
     \name Constructor, destructor, and defines.
  */
  /*@{*/

  /// Empty constructor.
  /**
   */
   BlockBaseRegister()
      : m_ghost(0), m_verbosity(0), m_isDefined(false), m_hasData(false), m_hasDefinedDataLayouts(false)
   {}

  /// Constructor.
  /**
     – a_mblock: Pointer to multiblock coordinate system.
     – a_grids: DisjointBoxLayout over which multiblock grid is defined.
     – a_ghost: Extension of register domain for auxiliary data.
   */
   BlockBaseRegister(const MultiBlockCoordSys*  a_mblock,
                        /// CELL-centered grids
                        const DisjointBoxLayout&   a_grids,
                        /// ghost width for state variables
                        const int&                 a_ghost = 0 );

  /// virtual destructor
  virtual ~BlockBaseRegister();

  /// Define function that matches constructor
  /**
     – a_mblock: Pointer to multiblock coordinate system.
     – a_grids: DisjointBoxLayout over which multiblock grid is defined.
     – a_ghost: Extension of register domain for auxiliary data.
  */
  void define(const MultiBlockCoordSys*  a_mblock,
              /// CELL-centered grids
              const DisjointBoxLayout&   a_grids,
              /// ghost width for state variables
              const int&                 a_ghost = 0 );

  /// Returns true if class has been defined, false otherwise.
  bool isDefined() const
  {
    return m_isDefined;
  }

  /// Set verbosity.
  /**
     Default is 0.  Positive values output increasing amounts of information.
   */
  virtual void setVerbose(int a_verbose)
  {
    m_verbosity = a_verbose;
  }

  /// Store data on block-boundary faces.
  /**
     – a_flux: Face-centered data to be stored.
     – a_dataIndex: Index of box where flux has been computed.
     – a_dir: Direction of the face.
     – a_side: Low or high side of the box where the block-boundary
     face lies.
  */
  void store(const T& a_data,
             const DataIndex& a_dataIndex,
             int a_dir,
             Side::LoHiSide a_side);

  /// Increment data on block-boundary faces.
  /**
     – a_flux: Face-centered flux to be stored.
     – a_dataIndex: Index of box where flux has been computed.
     – a_dir: Direction of the face.
     – a_side: Low or high side of the box where the block-boundary
     face lies.
  */
  void increment(const T& a_flux,
                 const DataIndex& a_dataIndex,
                 int a_dir,
                 Side::LoHiSide a_side);

  /// Exchange data to fill exterior register regions.
  /**
     Perform data manipulation so that both interior and exterior register
     data is current.
  */
  void exchange();

  ///
  /**
     Sets the contents of the register to zero.
  */
  void zeroRegister( const int a_n_comp=-1 );

  /// Returns whether or not a box has a block-boundary interface on one of its faces.
  /**
     – a_dataIndex: Index of box.
     – a_dir: Direction of the face of the box.
     – a_side: Low or high side of the box.
  */
  bool hasInterface(const DataIndex& a_dataIndex,
                    int a_dir,
                    Side::LoHiSide a_side) const;

  /// Return data on a block-boundary interface of a box, if hasInterface(a_dataIndex, a_dir, a_side) is true.
  /**
     – a_data: Face-centered data returned.
     – a_dataIndex: Index of box.
     – a_dir: Direction of the face.
     – a_side: Low or high side of the box.
     – a_sideData: Low or high side of the block-boundary interface,
     from which data will be taken.
  */
  void fill(T& a_data,
            const DataIndex& a_dataIndex,
            int a_dir,
            Side::LoHiSide a_side,
            Side::LoHiSide a_sideData) const;

  /// Set the name of the register.
  /**
     String identifier for this register, useful for debugging output.
     Default is UNDEFINED.
   */
  virtual void setName(const std::string& a_name)
  {
    m_name = a_name;
  }

protected:

  /// Returns true if the block-boundary in this direction on this side is an interface between blocks, false otherwise.
  /**
     – a_block: Index of the block.
     – a_dir: Direction of the face.
     – a_side: Low or high side of the box.
  */
  inline bool isInterface( const int& a_block,
                           const int& a_dir,
                           const Side::LoHiSide& a_side ) const;

  /// Returns true if data has been stored in the register, false otherwise.
  inline bool hasData() const
  {
    return m_hasData;
  }

  /// Returns the integer index of the block in which the box corresponding to the DataIndex resides.
  /**
     – a_dataIndex: Index of box.
  */
  inline int getBlockID( const DataIndex& a_dataIndex ) const;

  /// Returns the box corresponding to the block specified by the block integer index.
  /**
     – a_block: Index of the block.
  */
  inline const Box& getBlockBox( const int& a_block ) const;

  /// Returns the box corresponding to the block opposite to the block specified by the block integer index through the block-boundary specified by the integer index.
  /**
     – a_block: Index of the block.
     – a_block: Index of the block boundary.
  */
  inline const Box& getNeighborBlockBox( const int& a_block,
                                         const int& a_blockBdryIndex ) const;

  /// Returns the index transformation box corresponding to the Blockboundary specified by the block integer index and the block-boundary index.
  /**
     – a_block: Index of the block.
     – a_block: Index of the block boundary.
  */
  inline IndicesTransformation getTransformation(
     const int& a_block,
     const int& a_blockBdryIndex ) const;

  /// Returns the index of the neighboring block opposite the block in the direction and side spceified.
  /**
     – a_block: Index of the block.
     – a_dir: Direction of the face.
     – a_side: Low or high side of the box.
  */
  inline int getSrcBlockBdryIndex( const int& a_block,
                                   const int& a_dir,
                                   const Side::LoHiSide& a_side );

  /// Builds an array of vectors of boxes, one for each direction-side combination, that corresponds to either the full register, the src for external register regions, or the image of the src external regions.
  /**
     – a_boxes: Array of vector of boxes.
     – a_boxProcMap: Array of vector of processor mappings.
     – a_dataIndex: Array of vector of DataIndices that relate each box to a box in the underlying DisjointBoxLayout.
     – a_type: Type of register box to build [DST, SRC, or EXC].
  */
  enum UseType
  {
    Invalid=-1,
    DST=0,
    SRC,
    EXC,
    NUM_USE_TYPES
  };

  void buildBoxVectors( Vector<Box> a_boxes[2*CH_SPACEDIM],
                        Vector<int> a_boxProcMap[2*CH_SPACEDIM],
                        Vector<DataIndex> a_dataIndex[2*CH_SPACEDIM],
                        const UseType a_type );

  /// Builds a mapping from the block-boundary index of neighboring block to the local block-boundary index.
  /**
     – a_map: Array of vector relating block-boundary indices.
  */
  void buildInverseBlockBdryIndexMap( Vector<int> a_map[2*CH_SPACEDIM] );

  /// Build all box layouts for the data register
  void buildRegisterLayout();

  /// Build all box layouts for the source boxes
  void buildSrcRegisterLayout();

  /// Build all disjoint box layouts for the data exchange
  void buildExchangeRegisterLayout();

  /// Define each of the three internally used data layouts
  void defineDataLayouts( const int a_ncomp );

  /// Undefine each of the three internally used data layouts
  void clearDataLayouts();

  /// Build a box layout for the data register given a vector of boxes with processor mappings and related DataIndices from the underlying disjoint box layout.
  /**
     – a_layout: Resulting box layout.
     – a_map: Mapping between DataIndexes of the underlying disjoint box layout and DataIndexes of the resulting box layout.
     – a_boxes: Vector of boxes.
     – a_boxProcMap: Vector of processor mappings.
     – a_indices: Vector of DataIndices that relate each box to a box in the underlying DisjointBoxLayout.
  */
  inline void createDstBoxLayout( DisjointBoxLayout& a_layout,
                                  LayoutData<DataIndex>& a_map,
                                  const Vector<Box>& a_boxes,
                                  const Vector<int>& a_procs,
                                  const Vector<DataIndex>& a_indices );

  /// Build the box layout and associated box mappings for the source boxes
  /**
     – a_layout: Resulting box layout.
     – a_dstToSrcMap: Array of mappings between DataIndexes of the destination box layout and DataIndexes of the resulting box layout.
     – a_boxes: Vector of boxes.
     – a_boxProcMap: Vector of processor mappings.
     – a_indices: Vector of DataIndices that relate each box to a box in the underlying DisjointBoxLayout.
     – a_dstBlockBdryIndices: Vector of DataIndices that relate each box to a block boundry of the destination block.
     – a_grdToDstMap: Array of vector of DataIndices that relate each in the underlying DisjointBoxLayout to a box in the destination register.
  */
  inline void createSrcBoxLayout(
     DisjointBoxLayout& a_layout,
     LayoutData<DataIndex> a_dstToSrcMap[2*CH_SPACEDIM],
     const Vector<Box>& a_boxes,
     const Vector<int>& a_procs,
     const Vector<DataIndex>& a_indices,
     const Vector<int>& a_dstBlockBdryIndices,
     const LayoutData<DataIndex> a_grdToDstMap[2*CH_SPACEDIM] );

  /// Build the box layout and associated box mappings for the data exchange
  /**
     – a_layout: Resulting box layout.
     – a_dstToExcMap: Array of mappings between DataIndexes of the destination box layout and DataIndexes of the resulting box layout.
     – a_boxes: Vector of boxes.
     – a_boxProcMap: Vector of processor mappings.
     – a_indices: Vector of DataIndices that relate each box to a box in the underlying DisjointBoxLayout.
     – a_grdToDstMap: Array of vector of DataIndices that relate each in the underlying DisjointBoxLayout to a box in the destination register.
  */
  inline void createExchangeBoxLayout(
     DisjointBoxLayout& a_layout,
     LayoutData<DataIndex>& a_dstToExcMap,
     const Vector<Box>& a_boxes,
     const Vector<int>& a_procs,
     const Vector<DataIndex>& a_indices,
     const LayoutData<DataIndex>& a_grdToDstMap );

  /// Copy, including parallel exchange, of data from interior exchange registers to exterior source registers.
  inline void copyToSrcRegister();

  /// Local copy of data from exterior source registers to exterior full register.
  inline void copyToDstRegister();

  /// Print all registers to pout().
  void printRegisters() const;

  /// Print the specified register.
  /**
     – a_register: Array of registers to print.
  */
  inline void print(
     const LevelData<T> a_register[2*CH_SPACEDIM] ) const;
  inline void print(
     const BoxLayoutData<T> a_register[2*CH_SPACEDIM] ) const;

  inline bool hasStorage() const { return m_hasDefinedDataLayouts; }
      
  DisjointBoxLayout m_grids;
  ProblemDomain m_domain;
  MultiBlockCoordSys* m_coordsys;
  int m_ghost;
  int m_verbosity;
  bool m_isDefined;
  bool m_hasData;
  bool m_hasDefinedDataLayouts;
  std::string m_name;
  IntVect m_ghostVect;

 typedef std::map<int,int> BlockBdryIndexMap;
  Vector<BlockBdryIndexMap> m_blockBdryIndexMap;

  // Data register
//  BoxLayout m_registerLayout[2*CH_SPACEDIM];
//  BoxLayoutData<T> m_register[2*CH_SPACEDIM];
  DisjointBoxLayout m_registerLayout[2*CH_SPACEDIM];
  LevelData<T> m_register[2*CH_SPACEDIM];

  DisjointBoxLayout m_excRegisterLayout[2*CH_SPACEDIM];
  LevelData<T> m_excRegister[2*CH_SPACEDIM];

  DisjointBoxLayout m_srcRegisterLayout[2*CH_SPACEDIM];
  LevelData<T> m_srcRegister[2*CH_SPACEDIM];
//  BoxLayout m_srcRegisterLayout[2*CH_SPACEDIM];
//  BoxLayoutData<T> m_srcRegister[2*CH_SPACEDIM];

  LayoutData<DataIndex> m_grdToDstMap[2*CH_SPACEDIM];
  LayoutData<DataIndex> m_dstToExcMap[2*CH_SPACEDIM];
  LayoutData<DataIndex> m_dstToSrcMap[2*CH_SPACEDIM];
};

#include "NamespaceFooter.H"

#include "BlockBaseRegisterImplem.H"

#endif // include guard
