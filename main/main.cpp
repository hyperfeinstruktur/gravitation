#include <iostream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <valarray>
#include <iomanip>

#include "tree.h"
#include "body.h"
#include "engine.h"

using namespace std;

int main(int argc, char* argv[])
{
    // ===== Main Section ===== //
    try
    {
        engine engine(argc, argv);
        engine.run();
    }


    // ===== Exception Handling ===== //
    catch (int k)
    {
        cerr << "Fatal Error: ";
        if      ( k == 1) { cerr << "Problematic Bounds Provided to a node." << endl; }
        else if ( k == 2) { cerr << "Attempt to work with non 3D position." << endl; }
        else if ( k == 3) { cerr << "Attempt to compute force with massless body." << endl;}
        else if ( k == 4) { cerr << "Massless body added to node." << endl;}
        else if ( k == 5) { cerr << "Position of a body added to a node is outside its boundaries." << endl;}
        else if ( k == 6) { cerr << "Simulation canceled by user." << endl;}
        else              { cerr << "An unknown error occured." << endl; }

        cerr << "Terminating Program..." << endl;
        return k;
    }

  return 0;
}
