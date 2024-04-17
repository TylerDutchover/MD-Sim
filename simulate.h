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
    vector<double> get_acc();
    double get_mass();

    void update_pos(vector<double>);
    void update_vel(vector<double>);
    void update_acc(vector<double>);
    void set_mass(double);

private:
    vector<double> pos, vel, acc;
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
    void setAcc();


    // specific
    vector<double> calc_force(vector<double>, vector<double>);
      

    void run_sim();
    void update_pos(double);


    vector<double> sum(vector<double>, vector<double>);
    vector<double> scale(vector<double>, double);
    double dot_product(vector<double>, vector<double>);
    vector<double> reflect_v(vector<double>, vector<double>);


    void showParts();
    

private:
    double l_x, l_y, l_z, vol;
    double t_max, t;
    uint N;
    vector<particle> particles;
};

#endif