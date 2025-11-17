#include <CD_MultiIndex.H>
#include <CD_Driver.H>

using namespace ChomboDischarge;

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  constexpr int a_order = 3;

  MultiIndex mi(a_order);

  if (procID() == 0) {
    for (MultiIndex mi(a_order); mi.ok(); ++mi) {
      std::cout << mi.getCurrentIndex() << "\t" << mi.factorial() << std::endl;
    }
  }

  ChomboDischarge::finalize();
}
