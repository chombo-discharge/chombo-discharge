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
#include <iomanip>

#include <Utility.H>

// this mersenne_ran_main() outputs first 1000 generated numbers
// compare against the output of mt19937int.out
int
main(int argc, char** argv)
{
  BoxLib::Initialize(argc,argv);
  BoxLib::mt19937 rr(4357UL);
    std::ios::fmtflags ofmtflags = std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout << std::setprecision(8);
    for ( int j=0; j<1000; j++ )
    {
        std::cout << std::setw(10) << rr.u_value() << ' ';
        if ( j%5==4 ) std::cout << std::endl;
    }
    std::cout << std::endl;
    BoxLib::Finalize();
}
