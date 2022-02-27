#ifndef TREE_NODE_H
#define TREE_NODE_H
#pragma once
#include <array>
#include <memory>
#include "body.h"
#include "physical_constants.h"

using namespace std;

class tree_node
{
public:
    // Constructs new node with specified bounds
    tree_node(array<array<double,2>,3> const&);

    // Checks if a point is inside the node
    bool inside(valarray<double> const&) const; // check if position is inside cube represented by the node

    // Get info if external node
    bool external() const { return ext;}

    // divides the node into 8 children
    void divide();

    // Displays this node and subnodes
    void display(int const&) const;

    // For testing
    void divide_child(int k) {children[k]->divide();}

    // Adds a body to the node
    bool add_body(body const&); // true if body has been placed in an empty external node

    // To calculate the G force
    valarray<double> force_on_body(body const&) const;

    // Clears the node (i.e. deletes subnodes)
    void clear();

private:
    // boundaries(i,j) : i = x,y,z (0,1,2) j = inf,sup (0,1)
    array<array<double,2>,3> boundaries;

    // Width
    double s;
    // true: external node , false: internal node
    bool ext;

    // Does the node contain a body (true=yes)
    bool empty;
    // "body" contained in the node (either an actual particle or the center of mass)
    body center_of_mass;

    // children of the node. Null if external=true
    array<unique_ptr<tree_node>,8> children;

    // Displays bounds of this node only (used by display() recursively)
    void display_bounds(int const&) const;

    void update_com(body const&);

};

#endif // TREE_NODE_H
