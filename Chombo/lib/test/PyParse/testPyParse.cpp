#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PyParse.H"
#include "CH_assert.H"
#include <iostream>
using namespace std;

//-----------------------------------------------------------------------
int
main(int a_argc, char* a_argv[])
{
  vector<string> stuffWeWant;
  stuffWeWant.push_back("foo");
  stuffWeWant.push_back("goo");
  stuffWeWant.push_back("bar");
  stuffWeWant.push_back("F");
  PyParse parser("testPyParse.inputs", stuffWeWant);

  CH_assert(parser.parsed());
  CH_assert(parser.contains("foo"));
  CH_assert(parser.contains("goo"));
//  CH_assert(parser.contains("goo", PyParse::STRING));
  CH_assert(parser.contains("bar"));

  int foo;
  parser.get("foo", foo);
  Real foo2;
  parser.get("foo", foo2);
  string goo;
  parser.get("goo", goo);
  vector<int> bar;
  parser.get("bar", bar);
  IntVect bar2;
  parser.get("bar", bar2);

  // See if Boxes work properly.
  Box box1;
  parser.get("box1", box1);
  vector<Box> boxes;
  parser.get("boxes", boxes);
  vector<vector<Box> > boxHierarchy;
  parser.get("boxHierarchy", boxHierarchy);
  vector<vector<vector<Box> > > boxBlockHierarchy;
  parser.get("boxBlockHierarchy", boxBlockHierarchy);

  // Check functions.
  RefCountedPtr<PyUnaryFunction> F;
  parser.get("F", F);
  RefCountedPtr<PyBinaryFunction> G;
  parser.get("G", G);

  // Report our results.
  cout << "foo(int) = " << foo << endl;
  cout << "foo(Real) = " << foo2 << endl;
  cout << "goo = " << goo << endl;
  cout << "bar(vector<int>) = [" << bar[0];
  for (int i = 1; i < bar.size(); ++i)
    cout << ", " << bar[i];
  cout << "]" << endl;
  cout << "bar(IntVect) = " << bar2 << endl;
  cout << "F(2) = " << (*F)(2.0) << endl;
  cout << "G(1, 2) = " << (*G)(1.0, 2.0) << endl;
  cout << "box1 = " << box1 << endl;
  cout << "boxes = ";
  for (int i = 0; i < boxes.size(); ++i)
    cout << boxes[i] << " ";
  cout << endl;

  cout << "boxHierarchy = " << endl;
  for (int i = 0; i < boxHierarchy.size(); ++i)
  {
    cout << " " << i << ": ";
    for (int j = 0; j < boxHierarchy[i].size(); ++j)
      cout << boxHierarchy[i][j] << " ";
    cout << endl;
  }
  cout << endl;

  cout << "boxBlockHierarchy = " << endl;
  for (int i = 0; i < boxBlockHierarchy.size(); ++i)
  {
    cout << "Block " << i << ": " << endl;
    for (int j = 0; j < boxBlockHierarchy[i].size(); ++j)
    {
      cout << " Level " << j << ": ";
      for (int k = 0; k < boxBlockHierarchy[i][j].size(); ++k)
        cout << boxBlockHierarchy[i][j][k] << " ";
      cout << endl;
    }
    cout << endl;
  }
  cout << endl;

  PyParse::finalize(); // <- Needs to be called somewhere.

  return 0;
}
//-----------------------------------------------------------------------
