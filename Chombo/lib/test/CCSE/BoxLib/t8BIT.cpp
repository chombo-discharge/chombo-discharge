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
// A simple program to read in a MultiFab and write out in 8BIT format.
//

#include <string>

#include <VisMF.H>

int
main (int argc, char** argv)
{
    argc--; argv++;

    FArrayBox::setFormat(FABio::FAB_8BIT);

    for (int i = 0; i < argc; i++)
    {
        std::cout << "Transforming " << argv[i] << " ... " << std::flush;

        std::string name = argv[i];

        MultiFab mf;

        VisMF::Read(mf, name);

        VisMF::Write(mf, name, VisMF::OneFilePerCPU, true);

        std::cout << "done" << std::endl;
    }
}
