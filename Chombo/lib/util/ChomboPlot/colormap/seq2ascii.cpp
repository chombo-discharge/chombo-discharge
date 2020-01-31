#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//
// seq2ascii
//
// converts colormap file from seq format, produced by icol, to ascii rgb
// triplet format.  there must be num_colorslots color entries.
//
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

int
ReadSeqPalette (std::vector <unsigned char>& red,
                std::vector <unsigned char>& green,
                std::vector <unsigned char>& blue,
                const std::string& fileName
  );

void
WriteASCIIPalette (const std::vector <unsigned char>& red,
                   const std::vector <unsigned char>& green,
                   const std::vector <unsigned char>& blue,
                   const std::string& fileName
  );

int num_colorslots;

// -------------------------------------------------------------------
int main(int argc, char** argv)
{
  if (argc != 3)
  {
    cerr << "usage: seq2ascii seqfile asciifile" << endl;
    return (-1);
  }
  num_colorslots = 256;

  std::vector <unsigned char> red   (num_colorslots);
  std::vector <unsigned char> green (num_colorslots);
  std::vector <unsigned char> blue  (num_colorslots);

  std::string seqfile (argv[1]);
  ReadSeqPalette (red, green, blue, seqfile);

  std::string asciifile (argv[2]);
  WriteASCIIPalette (red, green, blue, asciifile);
}

// -------------------------------------------------------------------
int
ReadSeqPalette (std::vector <unsigned char>& red,
                std::vector <unsigned char>& green,
                std::vector <unsigned char>& blue,
                const std::string& fileName
  )
{
  unsigned char rbuff[num_colorslots];
  unsigned char gbuff[num_colorslots];
  unsigned char bbuff[num_colorslots];
  int   i, fd;          /* file descriptor */

  if ((fd = open(fileName.c_str(), O_RDONLY, NULL)) < 0)
  {
    cout << "Can't open colormap file:  " << fileName << endl;
    return (-1);
  }

  if (read(fd, rbuff, num_colorslots) != num_colorslots)
  {
    cout << "palette is not a seq colormap." << endl;
    return(-2);
  }
  if (read(fd, gbuff, num_colorslots) != num_colorslots)
  {
    cout << "file is not a seq colormap." << endl;
    return(-2);
  }
  if (read(fd, bbuff, num_colorslots) != num_colorslots)
  {
    cout << "palette is not a seq colormap." << endl;
    return(-2);
  }
  (void) close(fd);

  for (int i = 0; i < num_colorslots; ++i)
  {
    red[i]   = rbuff[i];
    green[i] = gbuff[i];
    blue[i]  = bbuff[i];
  }
  return(0);

}  // end ReadSeqPalette()

// -------------------------------------------------------------------
void
WriteASCIIPalette (const std::vector <unsigned char>& red,
                   const std::vector <unsigned char>& green,
                   const std::vector <unsigned char>& blue,
                   const std::string& fileName
  )
{
 CH_assert (red.size () == green.size () && red.size () == blue.size ());
  ofstream strm (fileName.c_str ());
  for (int i = 0; i < red.size (); ++i)
  {
    strm << int(red[i]) << "  " << int(green[i]) << "  " << int(blue[i]) << endl;
  }
  strm.close ();
}
