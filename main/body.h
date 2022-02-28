#ifndef BODY_H
#define BODY_H
#pragma once
#include <iostream>
#include <valarray>

using namespace std;

// A "body" is any object with a location in space, a velocity and a mass. This may be an
// actual fundamental body or an abstract body defined by the center of mass of a set of
// fundamental bodies. This is then used in the tree structure.
class body
{
public:
    // Constructs the body with initial position, velocity and mass
    body(double const&, double const&, double const&, \
         double const&, double const&, double const&, \
         double const&);

    // Getters for the engine class
    valarray<double> getpos() const;

    valarray<double> getvel() const;

    double getmass() const;

    // Setters for engine class
    void setpos(valarray<double> const&);

    void setvel(valarray<double> const&);

    void setmass(double const&);

    // For testing purposes (FIXME)
    void display() const;

private:

    valarray<double> position;

    valarray<double> velocity;

    double mass;
};

#endif // BODY_H
