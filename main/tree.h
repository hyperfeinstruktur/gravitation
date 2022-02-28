#ifndef TREE_H
#define TREE_H
#pragma once
#include <vector>
#include <cstdlib>
#include <sstream>
#include <thread>
#include <mutex>

#include "tree_node.h"
#include "body.h"

class tree
{
public:
    // root boundaries and gravitational constant to be used
    tree(array<array<double,2>,3> const&);

    // Builds the tree from the positions of the particles provided in argument
    void build(vector<body> const&);

    bool inside(valarray<double> const&) const;

    valarray<double> force_on_body(body const&) const;

    valarray<double> force_on_bodies(vector<body> const&, unsigned int const&) const;

    void display() const;


private:
    // Root of the tree (node at level 0)
    tree_node root;

    array<array<double,2>,3> boundaries;

    valarray<double> bodies_attraction(body const&, body const&) const;

    valarray<double> law_of_gravity(valarray<double> const&, valarray<double> const&, double const&, double const&) const;

    void force_thread(vector<body> const&,valarray<double>&,size_t,size_t) const;
};

#endif // TREE_H
