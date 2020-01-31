#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>

#include <Array.H>
#include <ParmParse.H>

int
main (int argc, char** argv)
{
    BoxLib::Initialize(argc,argv);

    ParmParse pp;

    Array<int> arr;

    pp.queryarr("arr",arr);

    for (int i = 0; i < arr.size(); i++)
    {
        std::cout << arr[i] << std::endl;
    }

    BoxLib::Finalize();
}
