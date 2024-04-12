#include "simulate.h"


// def in the main file...i guess it okay
extern double get_rand();


// for rectangular box
box::box(double x, double y, double z) { 
    l_x =   x;
    l_y =   y;
    l_z =   z; 
    vol = l_x * l_y * l_z;
    N = 0;  
}

// for a square box
box::box(double x) {
    l_x =   x;
    l_y =   x;
    l_z =   x;
    vol = l_x * l_y * l_z;
    N = 0;
}

// create all particles in the particles vector
void box::setNumParticles(double n) {
    N = n;
    particle myPar;
    for(uint i = 0; i < N; i++) {
        particles.push_back(myPar);
    }
}

// randomly distribute the atoms inside the box
void box::setPositions_random() {
    double x,y,z;
    for(uint i = 0; i < N; i++) {
        x = get_rand() * l_x;
        y = get_rand() * l_y;
        z = get_rand() * l_z;
        particles.at(i).update_pos({x,y,z});
    }
}

// max time for sim to run from 0 -> t_max
void box::setTime(double tMax) {
    t_max = tMax;
}

// sets the function pointer for potential
void box::setPotential(U_pot pot) {
    potential = pot;
}

void box::showParts() {
    for(uint i = 0; i < N; i++) {
        for(auto i : particles.at(i).get_pos())
            cout << setw(10) << setprecision(10) << left << i << setw(1) << ' ';
        cout << endl;
    }
}



//---------------------------------------------------------------//
particle::particle(vector<double> r, vector<double> v, double M) {
    pos     = r;
    vel     = v;
    mass    = M;
}

particle::particle(double M) {
    mass    = M;
    pos     = {0,0,0};
    vel     = {0,0,0};
}
    
particle::particle() {
    mass    = 0;
    pos     = {0,0,0};
    vel     = {0,0,0};
}

vector<double> particle::get_pos() {
    return pos;
}

vector<double> particle::get_vel() {
    return vel;
}

void particle::update_pos(vector<double> r) {
    pos = r;
}

void particle::update_vel(vector<double> V) {
    vel = V;
}

void particle::set_mass(double M) {
    mass = M;
}