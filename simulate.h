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
    void setVelocity();  
    void setTime(double, double);
    void calc_acceleration();
    void build_ring_for_all();
    void build_ring(uint);
    void show_ring();
    void run_sim();
    void update_pos(double);
    void init_walls();
    void special_setup();
    void set_temp(double t);

    double  calc_temp();
    double  calc_total_momentum();
    double  calc_pressure(double, double);
    bool    in_equillibrium();
    vector<double> calc_force(uint, uint);
      
    

    vector<double> sum(vector<double>, vector<double>);
    vector<double> scale(vector<double>, double);
    double dot_product(vector<double>, vector<double>);
    double sep_magnitude(vector<double>, vector<double>);
    vector<double> reflect_v(vector<double>, vector<double>);
    

private:
    double l_x, l_y, l_z, vol, M, t_max, t_step, r_min, r_max, p_offset_avg, p_offset;
    uint N;
    bool recalculate_ring, equillibrium;
    map<uint, uint> walls;
    map<uint, vector<particle*>> influence_ring;
    vector<particle> particles;
};

#endif