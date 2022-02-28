#pragma once
#include<iostream>
#include<valarray>
#include<vector>
#include <ctime>
#include<fstream>
#include "ConfigFile.h"
#include "body.h"
#include "tree.h"
#pragma once
using namespace std;

// This engine class is the heart of the simulation. It loads the parameters, initializes the bodies,
// constructs the initial tree and then executes the simulation loop (therby calling the specific
// functions of the other classes.
class engine
{
public:

    // Constructor
    engine(int argc, char* argv[]);

    ~engine();

    // Core function, called from main
    void run();

private:

    // Writes to output file
    void printOut(bool force);

    // Performs one timestep
    void leapfrog_step();

    // Runs the simulation (called by run())
    void run_leapfrog();

    // Methods to convert between body and position-velocity representation
    // e.g. position vector and velocity vector of the form r=(r1_1,r1_2,r1_3,r2_1,...)
    // TODO: bypass this conversion step as it takes up unnecessary memory

    // Sets the bodies parameters from pos/vel vectors
    void setbodies_pos(valarray<double> const&);
    void setbodies_vel(valarray<double> const&);

    // Creates pos/vel vectors from bodies parameters
    valarray<double> getbodies_pos() const;
    valarray<double> getbodies_vel() const;



private:
    // Fundamental Bodies in the simulation
    vector<body> bodies;

    // Spatial boundaries of the simulation, required to construct the tree
    array<array<double,2>,3> boundaries;

    // Tree
    tree octree;

    // Numerical parameters
    unsigned int nb_obj;  // Number of bodies
    double tF;            // End time
    int sampling;         // output every sampling-th step
    int last;             // For writing output
    double dt;            // Timestep for fixed alg
    double t;             // Current time
    ofstream *outputFile; // Output File
};

