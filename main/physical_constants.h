#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

namespace phys
{
    constexpr double G = 6.67408e-11;
    constexpr double newton_limit = 1.0e8;
    constexpr double barnes_hut_theta = 1.5;
    constexpr double softening_length = 1.0e27; //approx 60 AU squared
}
#endif // PHYSICAL_CONSTANTS_H
