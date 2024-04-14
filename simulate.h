#ifndef MDS
#define MDS

#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>

// may change later
using namespace std;
typedef double(*U_pot)(vector<double>);


class particle
{
public:
    particle(vector<double>, vector<double>, double);
    particle();
    particle(double);


    vector<double> get_pos();
    vector<double> get_vel();
    double get_mass();

    void update_pos(vector<double>);
    void update_vel(vector<double>);
    void set_mass(double);

private:
    vector<double> pos, vel;
    double mass;
};

class box
{
public:
    box(double,double,double);
    box(double);

    void setNumParticles(double, double);
    void setPositions_random();
    void setTemperature(double);  
    void setTime(double);
    void setPotential(U_pot);  


    void showParts(); // strictly for debug
    

private:
    double l_x, l_y, l_z, vol;
    double t_max, t;
    U_pot potential;
    uint N;
    vector<particle> particles;
};

#endif