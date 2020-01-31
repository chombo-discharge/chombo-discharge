#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CHOMBOPRINTER_H_
#define _CHOMBOPRINTER_H_

//
// The ChomboPrinter is a simple TestListener that
// just prints "OK" if everything goes well, otherwise
// reports the error in the format of compiler messages.
// The ChomboPrinter uses parstream.H's pout
//

#include <cxxtest/Flags.h>

#ifndef _CXXTEST_HAVE_STD
#   define _CXXTEST_HAVE_STD
#endif // _CXXTEST_HAVE_STD

#include <cxxtest/ErrorFormatter.h>
#include <cxxtest/StdValueTraits.h>

#include "parstream.H"

namespace CxxTest
{
  class ChomboPrinter : public ErrorFormatter
  {
  public:
    ChomboPrinter( CXXTEST_STD(ostream) &o = pout(), const char *preLine = ":", const char *postLine = "" )
      :
      ErrorFormatter( new Adapter(o), preLine, postLine )
    {
    }

    virtual ~ChomboPrinter()
    {
      delete outputStream();
    }

  private:
    class Adapter : public OutputStream
    {
      CXXTEST_STD(ostream) &_o;

    public:
      Adapter( CXXTEST_STD(ostream) &o )
        :
        _o(o)
      {
      }

      void flush()
      {
        _o.flush();
      }

      OutputStream &operator<<( const char *s )
      {
        _o << s;
        return *this;
      }

      OutputStream &operator<<( Manipulator m )
      {
        return OutputStream::operator<<( m );
      }

      OutputStream &operator<<( unsigned i )
      {
        char s[1 + 3 * sizeof(unsigned)];
        numberToString( i, s );
        _o << s;
        return *this;
      }
    };
  };
}

#endif // __cxxtest__ChomboPrinter_h__
