#include "engine.h"

// Conversion Factors
double pc_in_meters = 3.0857e16;
double solarmass_in_kg = 1.98847e30;
double km_to_m = 1000.0;

void disp(valarray<double> v);

void engine::printOut(bool force){
    if((!force && last>=sampling) || (force && last!=1))
    {
        // TODO: Bypass pos vector (iterate over bodies directly
      valarray<double> pos = getbodies_pos();
      //valarray<double> vel = getbodies_vel(); if necessary to print velocity
      for (size_t i(0) ; i<3*nb_obj ; ++i)
      {
        *outputFile << pos[i]/pc_in_meters << " ";
      }
      //for (size_t i(0) ; i<3*nb_obj ; ++i)
      //{
      //  *outputFile << vel[i] << " ";
      //}
            *outputFile << t << " " << endl;
      last=1;
    }
    else
    {
      last++;
    }
};



engine::engine(int argc, char* argv[])
    : octree({{{0,0},{0,0},{0,0}}},1)
{

    // ====== Global Parameter Initialization //

    // Input File
    string inputPath("config.in"); // Default configuration
    if(argc>1)                            // Custom input
    {
      inputPath = argv[1];
    }
    ConfigFile configFile(inputPath);
    for(int i(2); i<argc; ++i)
    {
        configFile.process(argv[i]);
    }

    // Output File
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);

    // Display
    cout << "Loaded configuration file   : '" << inputPath << "'." << endl;
    cout << "Writing to output file      : '";
    cout << configFile.get<string>("output").c_str() << "'." << endl;

    // Numerical Parameters
    tF     =                configFile.get<double>("tFin");
    dt =                 tF/configFile.get<double>("nsteps");
    sampling 		 =      configFile.get<double>("sampling");
    self_gravity = configFile.get<bool>("self_gravity");
     
    
    // Tree Init
    double inf = pc_in_meters*configFile.get<double>("inf");
    double sup = pc_in_meters*configFile.get<double>("sup");
    octree = tree({{ {inf,sup} , {inf,sup} , {inf,sup} }},self_gravity); // Cubic Domain
    nb_obj = configFile.get<double>("nb_obj");

    // Bodies Initialization
    for (size_t j(1) ; j<=nb_obj ; ++j)
    {
        body to_add = body( pc_in_meters*configFile.get<double>("y" + to_string(j) + "_" + "1") , \
                            pc_in_meters*configFile.get<double>("y" + to_string(j) + "_" + "2"), \
                            pc_in_meters*configFile.get<double>("y" + to_string(j) + "_" + "3"), \
                            km_to_m * configFile.get<double>("v" + to_string(j) + "_" + "1") , \
                            km_to_m * configFile.get<double>("v" + to_string(j) + "_" + "2"), \
                            km_to_m * configFile.get<double>("v" + to_string(j) + "_" + "3"), \
                            solarmass_in_kg*configFile.get<double>("m" + to_string(j))           );
        if (octree.inside(to_add.getpos()))
        {
            bodies.push_back(to_add);
        }
        else
        {
            cout << "Tree bounds: " << inf << " , " << sup << endl;
            disp(to_add.getpos());
            throw 5;
        }
    }

    unsigned int nb = configFile.get<double>("nsteps");
    cout << "Construction Done. Initialized Simulation with the following parameters: " << endl << endl;
    cout << "Number of Bodies:                           " << nb_obj << endl;
    cout << "Self Gravity (1: on, 0: off):               " << self_gravity << endl;
    cout << "Simulation endtime (millions of years):     " << tF/(86400 * 1e6) << endl;
    cout << "Number of steps:                            " << nb << endl;
    cout << "Timestep (seconds):                         " << dt << endl;
    cout << "Sampling of output file:                    " << sampling << endl;
    cout << "Theta Parameter of Barnes-Hut Algorithm:    " << phys::barnes_hut_theta << endl;
    cout << "Softening length:                           " << sqrt(phys::softening_length) / pc_in_meters << " pc" << endl;
    cout << "With these parameters, the output file size will be approximately " << (3*nb_obj*(nb/sampling)*20)/1e6 << " Mb. " << endl;
    //cout << "Proceed? (y: yes / any other key: no)";
    //char ans;
    //cin >> ans;
    //if (ans != 'y') {throw 6;}
    cout << endl << "Starting Simulation..." << endl << endl;
};






engine::~engine(){
    cout << "Destructor Called" << endl;
    outputFile->close();
    delete outputFile;
};

void disp(valarray<double> v)
{
    cout << "( ";
    for (size_t k(0) ; k<v.size(); ++k)
    {
      cout << v[k] << " , ";
    } cout << " )" << endl;
}

void engine::leapfrog_step()
{
    if (self_gravity)
    {
            valarray<double> xi = getbodies_pos();
            valarray<double> vi = getbodies_vel();
            octree.build(bodies);
            valarray<double> ai = octree.force_on_bodies(bodies,1);
            setbodies_pos(xi + dt*vi + 0.5*ai*dt*dt);
            octree.build(bodies);
            valarray<double> aip1 = octree.force_on_bodies(bodies,1);
            setbodies_vel(vi + 0.5*(ai+aip1)*dt);
    }
    else
    {
            valarray<double> xi = getbodies_pos();
            valarray<double> vi = getbodies_vel();
            valarray<double> ai = octree.force_on_bodies(bodies,1);
            setbodies_pos(xi + dt*vi + 0.5*ai*dt*dt);
            valarray<double> aip1 = octree.force_on_bodies(bodies,1);
            setbodies_vel(vi + 0.5*(ai+aip1)*dt);
    }
}

void engine::run_leapfrog()
{
    double prop = 10;
    clock_t begin = clock();
    float elapsed;
    while( t<(tF-0.5*dt) )
    {
        leapfrog_step();
        t += dt;
        printOut(false);
        elapsed = double(clock() - begin)/CLOCKS_PER_SEC;
        if ((int)(t*100/tF) == prop)
        {
            cout << "Simulation " << (t/tF)*100 << "% complete. Elapsed time: " << elapsed << " seconds, (" << elapsed/60.0 << " minutes). " << endl;
            //cout << "==> Remaining time is approximately "<< 1/(t/tF)*elapsed << "seconds, (" << (1-t/tF)*elapsed/60.0 << " minutes)." << endl;
            if (prop == 10) { cout << "==> The simulation will take approximately " << 10*elapsed << " seconds, (" << 10*elapsed/60.0 << " minutes). " << endl; }
            //if (prop == 100) {cout << endl << "Simulation fihisned. Total computing time: " << elapsed << "seconds, (" << elapsed/60.0 << " minutes). " << endl;}
            prop+=10;
        }
    }
}

void engine::setbodies_pos(valarray<double> const& v)
{
    for (size_t j(0) ; j<bodies.size() ; ++j)
    {
        bodies[j].setpos(v[slice(3*j,3,1)]);
    }
}

void engine::setbodies_vel(valarray<double> const& v)
{
    for (size_t j(0) ; j<bodies.size() ; ++j)
    {
        bodies[j].setvel(v[slice(3*j,3,1)]);
    }
}

valarray<double> engine::getbodies_pos() const
{
    valarray<double> temp = valarray<double>(3*nb_obj);
    for (size_t j(0) ; j<bodies.size() ; ++j)
    {
        temp[slice(3*j,3,1)] = bodies[j].getpos();
    }
    return temp;
}

valarray<double> engine::getbodies_vel() const
{
    valarray<double> temp = valarray<double>(3*nb_obj);
    for (size_t j(0) ; j<bodies.size() ; ++j)
    {
        temp[slice(3*j,3,1)] = bodies[j].getvel();
    }
    return temp;
}

void engine::run(){
    last = 0;
    t = 0;
    printOut(true);
    run_leapfrog();
    printOut(true);
};
