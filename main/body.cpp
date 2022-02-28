#include "body.h"

body::body(double const& x, double const& y, double const& z, \
     double const& vx, double const& vy, double const& vz, \
     double const& m)
    : position({x,y,z}) , velocity({vx,vy,vz}) , mass(m)
{}

valarray<double> body::getpos() const
{
    return position;
}

valarray<double> body::getvel() const
{
    return velocity;
}

double body::getmass() const
{
    return mass;
}


void body::setpos(valarray<double> const& pos)
{
    position = pos;
}

void body::setvel(valarray<double> const& vel)
{
    velocity = vel;
}

void body::setmass(double const& m)
{
    mass = m;
}

void body::display() const
{
    cout << "Position: ( " << position[0] << " , " << position[1] << " , " << position[2] << " )" << endl;
    cout << "Velocity: ( " << velocity[0] << " , " << velocity[1] << " , " << velocity[2] << " )" << endl;
    cout << "Mass:       " << mass << " kg " << endl << endl;
}
