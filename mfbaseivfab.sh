#Inclusion guards
for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py \)`; do
    sed -i 's/\#include \"units.H\"/\#include <CD_Units.H>/g' $i
    sed -i 's/\#include <units.H>/\#include <CD_Units.H>/g' $i
done

for i in `find . -type f \( -iname \*.H -o -iname \*.cpp -o -iname \*.py -o -iname \*.inputs -o -iname \*.options -o -iname *GNUmakefile \)`; do
    sed -i 's/namespace units/namespace Units/g' $i
    sed -i 's/units::/Units::/g' $i
    
    sed -i 's/const static Real s_eps0/constexpr Real eps0/g' $i
    sed -i 's/Units::s_eps0/Units::eps0/g' $i

    sed -i 's/const static Real s_Qe/constexpr Real Qe/g' $i
    sed -i 's/Units::s_Qe/Units::Qe/g' $i

    sed -i 's/const static Real s_c0/constexpr Real c0/g' $i
    sed -i 's/Units::s_c0/Units::c/g' $i

    sed -i 's/const static Real s_kb/constexpr Real kb/g' $i
    sed -i 's/Units::s_kb/Units::kb/g' $i

    sed -i 's/const static Real s_Td/constexpr Real Td/g' $i
    sed -i 's/Units::s_Td/Units::Td/g' $i

    sed -i 's/const static Real s_Na/constexpr Real Na/g' $i
    sed -i 's/Units::s_Na/Units::Na/g' $i

    sed -i 's/const static Real s_R/constexpr Real R/g' $i
    sed -i 's/Units::s_R/Units::R/g' $i

    sed -i 's/const static Real s_R/constexpr Real R/g' $i
    sed -i 's/Units::s_R/Units::R/g' $i

    sed -i 's/const static Real s_atm2pascal/constexpr Real atm2pascal/g' $i
    sed -i 's/Units::s_atm2pascal/Units::atm2pascal/g' $i

    sed -i 's/const static Real s_pi/constexpr Real pi/g' $i
    sed -i 's/Units::s_pi/Units::pi/g' $i

    sed -i 's/const static Real s_me/constexpr Real me/g' $i
    sed -i 's/Units::s_me/Units::me/g' $i

    
    sed -i 's/const static Real s_eV   = s_Qe/constexpr Real eV = 1.6021766208E-19/g' $i
    sed -i 's/Units::s_eV/Units::eV/g' $i
done

# Move files
mv Source/Utilities/units.H      Source/Utilities/CD_Units.H
