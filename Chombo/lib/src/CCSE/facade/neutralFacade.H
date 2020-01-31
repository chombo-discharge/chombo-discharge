#ifndef _NEUTRALFACADE_H_
#define _NEUTRALFACADE_H_

//
// This is a file that can be #included from both Chombo and CCSE files.
// So don't #include any Chombo or CCSE files *here*.
//

#ifdef CH_USE_DOUBLE
typedef double Real;
#else
typedef float Real;
#endif

//
// Forward declarations
//
namespace Chombo
{
  class IntVect;
  class Box;
  class BoxLayout;
  class DisjointBoxLayout;
  template<typename T> class LayoutData;
  template<typename T> class LevelData;
  class FArrayBox;
  class DataIterator;
}

class IntVect;  // CCSE
class Box;      // CCSE
class BoxArray; // CCSE
class MultiFab; // CCSE
class MFIter;   // CCSE

// There's no such namespace in the CCSE code, and that's precisely the point;
// we're creating one for exclusive use in the facade code, just for clarity.
// Note no actual CCSE code was harmed in the process.
namespace CCSE
{
  typedef ::IntVect IntVect;
  typedef ::Box Box;
  typedef ::BoxArray BoxArray;
  typedef ::MultiFab MultiFab;
  typedef ::MFIter MFIter;
}

#include <iosfwd>
#include <vector>

//
// The external interface
//
namespace BLfacade
{
class Box;
class BoxArray;
class MultiFab;


struct Streamable
{
  virtual ~Streamable() {}
  virtual std::ostream& streamOut( std::ostream& out ) const;
};

/**
  Chombo code that would like to use a CCSE IntVect should use this instead.
*/
class IntVect : public Streamable
{
public:
  IntVect( const Chombo::IntVect& );
  IntVect( const CCSE::IntVect& );
  IntVect( const IntVect& );
  IntVect& operator=( const IntVect& );
  virtual ~IntVect();
  operator Chombo::IntVect() const;
  operator CCSE::IntVect() const;

  virtual std::ostream& streamOut( std::ostream& ) const;

private:
  void initPimpl();

  int* m_serialization;
  CCSE::IntVect* m_pimpl;
};


class Box : public Streamable
{
public:
  Box( const Chombo::Box& chBox );
  Box( const CCSE::Box& ccseBox );
  Box( const Box& );
  Box& operator=( const Box& b );
  virtual ~Box();
  operator Chombo::Box() const;  
  operator CCSE::Box() const;
  virtual std::ostream& streamOut( std::ostream& ) const;
  bool operator<( Box const& that ) const;
  bool operator==( Box const& that ) const;

private:
  void initPimpl();

  IntVect m_lo;
  IntVect m_hi;
  CCSE::Box* m_pimpl;
};


class BoxArray : public Streamable
{
public:
  BoxArray( const Chombo::BoxLayout& chBL );
  virtual ~BoxArray();
  operator CCSE::BoxArray const &() const;
  virtual std::ostream& streamOut( std::ostream& ) const;

  int numBoxes() const { return m_numBoxes; }
  Chombo::Box operator[](int b) const;

  // For use by MultiFab only:
  Chombo::DisjointBoxLayout* dbl() const { return m_dbl; }
  // For use by MultiFab only:
  void dbl( Chombo::DisjointBoxLayout* a_dbl ) { m_dbl = a_dbl; }

private:
  void initPimpl();

  // Deliberately unimplemented.  Would be nice if it can last like this.
  BoxArray( const BoxArray& );
  BoxArray& operator=( const BoxArray& );

  int                        m_numBoxes;
  Box**                      m_boxes;
  CCSE::BoxArray*            m_pimpl;
  Chombo::DisjointBoxLayout* m_dbl;
};

// Used inside BLfacade::MultiFab.
struct BoxNDptr
{
  BoxNDptr( Box const& a_box, Real* a_dptr )
   : box(a_box), dptr(a_dptr) {}
  Box   box;
  Real* dptr;
  bool operator<( BoxNDptr const& ) const;
};

struct BoxNDit;


class MultiFab
{
public:
  MultiFab( BoxArray&, int ncomps, int nghosts,
            const std::vector<int>* procIDs=0);
  ~MultiFab();
  operator CCSE::MultiFab&() const;
  Chombo::LevelData<Chombo::FArrayBox>* pLevelData( const Chombo::Box& domain,
                                                    const Chombo::IntVect& prd);

private:
  void initPimpl( const std::vector<int>* procIDs );
  void chomboDelete();
  std::vector<int> procIDs() const;

  // Utility functions for constructing a Chombo::LevelData.
  std::vector<BoxNDptr> getBoxNDptrs() const;
  std::vector<BoxNDit>  getBoxNDits(const Chombo::LayoutData<Real*>&) const;

  // Deliberately unimplemented.  Would be nice if it can last like this.
  MultiFab( const MultiFab& );
  MultiFab& operator=( const MultiFab& );

  BoxArray&          m_boxarray; // Net of ghost cells
  std::vector<Real*> m_farrays;  // FAB data
  CCSE::MultiFab*    m_pimpl;
  int                m_ncomps;
  int                m_nghosts;
  Chombo::LevelData<Chombo::FArrayBox>* m_pLevelData; // cache
};

void initializeBoxLib(int argc, char** argv);
void finalizeBoxLib();

} // namespace BLfacade

std::ostream& operator<<(std::ostream&, const BLfacade::Streamable&);

#endif // include guard
