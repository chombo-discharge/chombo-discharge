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

#include <Thread.H>
#include <Profiler.H>
#include <ParallelDescriptor.H>

namespace BoxLib3
{
namespace testing
{
void profiler_main(int& argc, char**& argv);
}
}

extern "C"
void* test_profiler_routine(void*)
{
    CH_PROFILE("tp_routine");
    // BoxLib3::Thread::sleep(1.0);
    return 0;
}

namespace
{
void
thread_timing()
{
    CH_PROFILE("a_thread_timing()");
    FunctionThread ft(test_profiler_routine);
    FunctionThread ft1(test_profiler_routine);
    test_profiler_routine(0);
    ft.join();
    ft1.join();
}
}

int
main(int argc, char** argv)
{
    BoxLib::Initialize(argc, argv);
    const double million = 1000000.0;

    BoxLib::WallTimer wt;
    wt.start();
    wt.stop();
    std::cout << "Wall timer reports (us) = " << wt.accum_time()*million << std::endl;
    std::cout << "Wall timer tick = (us) " << wt.tick()*million << std::endl;
    thread_timing();

    BoxLib::Finalize();
}
