#include "tree.h"

mutex force_lock;

tree::tree(array<array<double,2>,3> const& bounds)
    :root(bounds) , boundaries(bounds)
{}

bool tree::inside(valarray<double> const& v) const
{
    return root.inside(v);
}

void tree::build(vector<body> const& bodies)
{
    root.clear();
    root = tree_node(boundaries);
    for (auto const& B : bodies)
    {
        if (root.inside(B.getpos())) {root.add_body(B);}
    }
}

void tree::display() const
{
    root.display(0);
}


valarray<double> tree::force_on_body(body const& B) const
{
    return root.force_on_body(B) / B.getmass();
}

void tree::force_thread(const vector<body> & bodies,valarray<double> & forces, size_t start, size_t end) const
{
/*
    unsigned int size = 3*(end-start+1);
    valarray<double> temp = valarray<double>(size);
    size_t k(0);
    for (size_t j(start) ; j<=end ; ++j)
    {
        temp[slice(k,3,1)] = force_on_body(bodies[j]);
        k+=3;
    }
    //force_lock.lock();
    forces[slice(3*start,size,1)] = temp;
    //force_lock.unlock();
    //cout << "A thread has finished." << endl;
*/

    for (size_t j(start) ; j<=end ; ++j)
    {
        //force_lock.lock();
        forces[slice(3*j,3,1)] = force_on_body(bodies[j]);
        //force_lock.unlock();
    }

}

valarray<double> tree::force_on_bodies(const vector<body> & bodies,unsigned int const& nbt) const
{
    valarray<double> forces(3*bodies.size());
    /*

    unsigned int nb_threads = nbt;
    vector<thread> threadvect;
    size_t spread = bodies.size() / nb_threads;
     unsigned int mod = bodies.size()%nb_threads;
    size_t newend = spread -1;
    size_t start = 0;
    for (unsigned int k = 0 ; k < nb_threads ; ++k)
    {
        if (k == nb_threads -1) {newend+=mod;}
        threadvect.emplace_back(&tree::force_thread,this,ref(bodies),ref(forces),start,newend);
        start += spread;
        newend += spread;
    }

    for (auto& t : threadvect)
    {
        t.join();
    }
    */
    for (size_t k = 0;k<bodies.size();++k)
    {
        forces[slice(3*k,3,1)] = force_on_body(bodies[k]);
    }
    return forces;

}
