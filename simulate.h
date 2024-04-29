#ifndef MDS
#define MDS

#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <map>

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
    void setInfluenceZone(double, double);
    void setPositions_random();
    void setVelocity(double);  
    void setTime(double, double);
    void setAcc();
    void build_ring_for_all();
    void build_ring(uint);

    vector<double> calc_force(vector<double>, vector<double>);
      
    void run_sim();
    void update_pos(double);

    vector<double> sum(vector<double>, vector<double>);
    vector<double> scale(vector<double>, double);
    double dot_product(vector<double>, vector<double>);
    double sep_magnitude(vector<double>, vector<double>);
    vector<double> reflect_v(vector<double>, vector<double>);
    

private:
    double l_x, l_y, l_z, vol, M, t_max, t_step, r_min, r_max;
    uint N;
    bool recalculate_ring;
    map<uint, vector<particle>> influence_ring;
    vector<particle> particles;
};

#endif