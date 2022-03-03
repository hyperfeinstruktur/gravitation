#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

namespace phys
{
    constexpr double G = 6.67430e-11;
    constexpr double newton_limit = 1.0e8;
    constexpr double barnes_hut_theta = 1.5;
    constexpr double softening_length = 0.0; //1.5e32; // ATTENTION: softening length SQUARED, 1pc = 3.0857e16 m
}
#endif // PHYSICAL_CONSTANTS_H
