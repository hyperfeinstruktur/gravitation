#include "tree_node.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <memory>

// Utilitaries, defined below
double           norm_3d(valarray<double> const& v);
double           grav_pot(valarray<double> const&, valarray<double> const&, double const&, double const&);
valarray<double> grav_force(valarray<double> const&, valarray<double> const&, double const&, double const&);
valarray<double> grav_force_pot(valarray<double> const&, valarray<double> const&, double const&, double const&);



tree_node::tree_node(array<array<double,2>,3> const& bounds)
    : boundaries(bounds) , ext(true) , empty(true) ,
      center_of_mass(0,0,0,0,0,0,0) // Massless by default
{
    for (auto& ptr : this->children)
    {
        ptr = nullptr;
    }
    //if (!check_bounds(bounds))
    //{
    //    throw 1;
    //}
    s = abs(boundaries[0][1]-boundaries[0][0]);
}

void tree_node::clear()
{
    for (auto& sub : children)
    {
        sub.reset();
    }
}
bool tree_node::inside(valarray<double> const& v) const
{
    if (v.size() > 3) { throw 2;}
    for (size_t j(0) ; j<3 ; ++j)
    {
        if ( (v[j] < boundaries[j][0] ) || ( v[j] > boundaries[j][1] ) )
        {
            //cout << v[j] <<" < "<< boundaries[j][0] << endl << "or" << endl;
            //cout << v[j] <<" > "<< boundaries[j][1] << endl;
            return false;
        }
    }
    return true;
}

void tree_node::divide()
{
    // Node marked as internal
    ext = false;

    // "Halving" the cube in each direction
    double mesh = (boundaries[0][1] - boundaries[0][0])/2.0;
    array<array<array<double,2>,3>,8> new_bounds;
    double x_midplane = boundaries[0][0] + mesh;
    double y_midplane = boundaries[1][0] + mesh;
    double z_midplane = boundaries[2][0] + mesh;

    // Bounds of the 8 subcubes
    new_bounds[0] = { { {boundaries[0][0],x_midplane} , {boundaries[1][0],y_midplane} , {boundaries[2][0],z_midplane} } };
    new_bounds[1] = { { {boundaries[0][0],x_midplane} , {boundaries[1][0],y_midplane} , {z_midplane,boundaries[2][1]} } };
    new_bounds[2] = { { {boundaries[0][0],x_midplane} , {y_midplane,boundaries[1][1]} , {boundaries[2][0],z_midplane} } };
    new_bounds[3] = { { {boundaries[0][0],x_midplane} , {y_midplane,boundaries[1][1]} , {z_midplane,boundaries[2][1]} } };
    new_bounds[4] = { { {x_midplane,boundaries[0][1]} , {boundaries[1][0],y_midplane} , {boundaries[2][0],z_midplane} } };
    new_bounds[5] = { { {x_midplane,boundaries[0][1]} , {boundaries[1][0],y_midplane} , {z_midplane,boundaries[2][1]} } };
    new_bounds[6] = { { {x_midplane,boundaries[0][1]} , {y_midplane,boundaries[1][1]} , {boundaries[2][0],z_midplane} } };
    new_bounds[7] = { { {x_midplane,boundaries[0][1]} , {y_midplane,boundaries[1][1]} , {z_midplane,boundaries[2][1]} } };

    // Creating subnodes and linking pointers of node to them
    for (size_t h(0) ; h<8 ; ++h)
    {
        children[h] = unique_ptr<tree_node>(new tree_node(new_bounds[h]));
    }
}

bool tree_node::add_body(body const& B)
{
    //if (!B.massive()) { throw 4; }
    //if (!inside(B.getpos())) { throw 5; } // body must be inside bounds of node

    if (inside(B.getpos())) // Only add body if it is inside node
    {
        if (empty)
        {
            empty = false;
            center_of_mass = B;
            return true;
        }
        else if (not ext)
        {
            update_com(B);
            bool placed = false;
            for (auto const& sub : children)
            {
                if (not placed) {placed = sub->add_body(B);}
            }
        }
        else
        {
            this->divide();
            bool p1 = false; bool p2 = false;
            for (auto const& sub : children)
            {
                if (not p1) {p1 = sub->add_body(this->center_of_mass);}
                if (not p2) {p2 = sub->add_body(B);}
            }
            update_com(B);
        }
    }
    return false;
}

void tree_node::update_com(const body & B)
{
    valarray<double> mr = center_of_mass.getmass()*center_of_mass.getpos(); // MR = sum(m_i*r_i)
    mr += B.getmass()*B.getpos();                                           // Adds m_n+1*r_n+1
    center_of_mass.setmass(center_of_mass.getmass() + B.getmass());         // Increases mass M
    center_of_mass.setpos( ( 1.0/center_of_mass.getmass() ) * mr);          // Sets R = 1/M*sum(m_i*r_i)
}


// valarray<double> tree_node::force_on_body(body const& B) const
// {
//     if (not empty)
//     {
//         double d = norm_3d(B.getpos() - center_of_mass.getpos());
//         if ((d < phys::newton_limit)) { return {0,0,0};} // TODO: remove newton_limit (replaced by softening)

//         if (ext || (s/d < phys::barnes_hut_theta))
//         {
//             //return bodies_attraction(B,center_of_mass);
//             return -grav_force(B.getpos(), center_of_mass.getpos(), B.getmass(), center_of_mass.getmass());
//         }
//         else
//         {
//             valarray<double> F = {0,0,0};
//             for (auto const& sub : children)
//             {
//                 F+=sub->force_on_body(B);
//             }
//             return F;
//         }
//     }
//     else return {0,0,0};
// }
valarray<double> tree_node::force_on_body(body const& B) const
{
    double d = norm_3d(B.getpos() - center_of_mass.getpos());
    if ((d < phys::newton_limit)) { return {0,0,0};} // TODO: remove newton_limit (replaced by softening)

    if (ext || (s/d < phys::barnes_hut_theta))
    {
        //return bodies_attraction(B,center_of_mass);
        return -grav_force(B.getpos(), center_of_mass.getpos(), B.getmass(), center_of_mass.getmass());
    }
    else
    {
        valarray<double> F = {0,0,0};
        for (auto const& sub : children)
        {
            if (not sub->empty) {F+=sub->force_on_body(B);}
        }
        return F;
    }

}

// ==== UTILITARIES === //

// 3-norm
double norm_3d(valarray<double> const& v)
{
    return sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// Gravitational Potential
double grav_pot(valarray<double> const& i , valarray<double> const& j, double const& mj)
{
    if ( mj>1e-10 )
    {
        return -phys::G * mj / norm_3d(i-j);
    }
    else
    {
        throw 3;
    }
}

// Gravitational Force
valarray<double> grav_force(valarray<double> const& i , valarray<double> const& j, double const& mi, double const& mj)
{
    if ( mi>1e-10 && mj>1e-10 )
    {
        return phys::G*mi*mj*(i-j) / pow( pow(norm_3d(i-j),2) + phys::softening_length , 1.5);
        //double d = norm_3d(i-j);
        //return phys::G*mi*mj*(i-j) / pow( d*d + phys::softening_length , 1.5);
        //valarray <double> F = (-phys::G*mi*mj) * pow( (1/norm_3d(i-j)) , 3 ) * (i-j);
        //double r_temp = norm_3d(i-j);
        //valarray <double> F = (-phys::G*mi*mj) * pow((pow(r_temp,2) + phys::softening_length),-1.5) * (i-j)/r_temp;
        //return -F;
    }
    else
    {
        throw 3;
    }
};

// Gravitational Force + Potential (returns size 4 vector: (Fx,Fy,Fz,phi)) TODO: implement softening
valarray<double> grav_force_pot(valarray<double> const& i , valarray<double> const& j, double const& mi, double const& mj)
{
    if ( mi>1e-10 && mj>1e-10 )
    {
        valarray<double> temp(4);
        double d = norm_3d(i-j);
        temp[3] = -phys::G * mj / d ;
        temp[slice(0,3,1)] = -phys::G*mi*mj*(i-j) / pow(d,3);
        return temp;
    }
    else
    {
        throw 3;
    }
};



//////////// Old Testing methods
void tree_node::display_bounds(int const& lvl) const
{
    string spaces;
    for (int j(0) ; j<lvl ; ++j) {spaces+="====";}
    cout << spaces << "x in [" << boundaries[0][0] << " , " << boundaries[0][1] << "]" << endl;
    cout << spaces << "y in [" << boundaries[1][0] << " , " << boundaries[1][1] << "]" << endl;
    cout << spaces << "z in [" << boundaries[2][0] << " , " << boundaries[2][1] << "]" << endl;
}


void tree_node::display(int const& lvl) const
{
    string spaces;
    for (int j(0) ; j<lvl ; ++j) {spaces+="====";}
    cout << spaces << "Bounds of Node: " << endl;
    display_bounds(lvl); cout << endl;
    cout << spaces << "Center of Mass: " << endl;
    cout << spaces << center_of_mass.getmass() << " kg, located at [" << center_of_mass.getpos()[0] << " , " << center_of_mass.getpos()[1] << " , " << center_of_mass.getpos()[2] << "]" << endl << endl;
    if (not external())
    {
        cout << spaces << "Bounds of children:" << endl;
        for (auto const& node : children)
        {
            node->display(lvl+1);
            cout << endl << endl;
        }
    }
    else
    {
       cout << spaces << "This node is external" << endl;
    }
}