#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <Utility.H>

int
main (int argc, char** argv)
{
    if (argc == 2)
    {
        if (!BoxLib::UtilCreateDirectory(argv[1], 0755))
        {
            std::cout << "Utility::UtilCreateDirectory() failed!!!\n";
        }
    }
}
