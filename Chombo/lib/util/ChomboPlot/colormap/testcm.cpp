#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <stream.h>
#include <fstream.h>
int main(int argc, char* argv[])
{
  fstream is;
  is.open("colormap.rgb", ios::out);

  for (int i =0; i < 256; i++)
    {
      int ir, ig, ib;
      ir =(256*i)/256;
      ig = (256*(256-i))/256;
      if (i < 128)
        ib = (128*(128-i))/128;
      else
        ib = (128*(i-128))/128;
      is  << ir << "   "
          << ig << "   "
          << ib << endl;
    }
  is.close();
}
