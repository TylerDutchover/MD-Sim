#include "simulate.h"
#include "const.h"

// def in the main file...i guess it okay
extern double get_rand();


// for rectangular box
box::box(double x, double y, double z) { 
    l_x =   x;
    l_y =   y;
    l_z =   z; 
    vol = l_x * l_y * l_z;
    N = 0; 
    t_max = 0; 
}

// for a square box
box::box(double x) {
    l_x =   x;
    l_y =   x;
    l_z =   x;
    vol = l_x * l_y * l_z;
    N = 0;
    t_max = 0;
}

// create all particles in the particles vector
void box::setNumParticles(double n, double mass) {
    N = n;
    particle myPar(mass);
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

// use the rms speed for init ?
void box::setTemperature(double T) {
    double v_rms;
    v_rms = (3*K_B*T) / particles.at(0).get_mass();
    v_rms = sqrt(v_rms);
    v_rms = v_rms / sqrt(3);
    for(uint i = 0; i < N; i++) {
        particles.at(i).update_vel({v_rms, v_rms, v_rms});
    }

}

// max time for sim to run from 0 -> t_max
void box::setTime(double tMax) {
    t_max = tMax;
}

void box::setAcc() {

    vector<double> total_force, force, r1, r2, accl;
    double M = particles.at(0).get_mass();

    for(uint i = 0; i < N; i++) 
    {
        total_force = {0,0,0};                              // total force on a particle due to the other N-1 particles
        r1 = particles.at(i).get_pos();                               // grab positional vector of particle

        for(uint j = 0; j < N; j++) 
        {
            if(j != i) 
            {
                r2 = particles.at(j).get_pos();                 // grab second postitional vector
                force = calc_force(r1, r2);                 // calculate the force between the two
                for(uint q = 0; q < 3; q++)
                    total_force.at(q) += force.at(q);       // sum total force
            }
        }

        accl = {total_force.at(0) / M, total_force.at(1) / M, total_force.at(2) / M};
        particles.at(i).update_acc(accl);
    }
}


// calculates the force between two particles from the lj potential
vector<double> box::calc_force(vector<double> r1, vector<double> r2) {
    vector<double> f_LJ;
    double  r_x  =  r2.at(0) - r1.at(0),
            r_y  =  r2.at(1) - r1.at(1),
            r_z  =  r2.at(2) - r1.at(2),
            r,                              // magnitude of seperation vector        
            A,                              // sigma    / r
            B,                              // epsilon  / r
            LJ_F_scale;                     // prefactor for LJ force
        
    r = sqrt(r_x*r_x + r_y*r_y + r_z*r_z); 
    f_LJ = {(r_x / r), (r_y / r), (r_z / r) };

    A =     SGMA    / r;
    B =     EPS     / r;  
    LJ_F_scale = 24 * B * ( (2*pow(A, 12)) - pow(A,6) );

    f_LJ = { (LJ_F_scale * f_LJ.at(0)), (LJ_F_scale * f_LJ.at(1)), (LJ_F_scale * f_LJ.at(2)) };
    return f_LJ;
}

void box::update_pos(double time_step) {
    vector<double> r_tmp, v_tmp, a1, a0;

    for(uint i = 0; i < N; i++) {
        r_tmp = particles.at(i).get_pos();
        v_tmp = particles.at(i).get_vel();
        a1 = particles.at(i).get_acc();
        a0 = a1;

        if(i == 0) {
            for(uint q = 0; q < 3; q++)
                cout << r_tmp.at(q) << ',';
            cout << '\n';
        }

        // update postitional vectors before anything else
        for(uint q = 0; q < 3; q++)
            r_tmp.at(q) = r_tmp.at(q) + (v_tmp.at(q)*time_step) + (0.5*a1.at(q) * pow(time_step,2));
        
        // recalculate the force with the updated positions, then update velocity vector
        setAcc();
        for(uint q = 0; q < 3; q++)
            v_tmp.at(q)  = v_tmp.at(q) + (0.5* (a1.at(q) + a0.at(q)) * time_step);

        particles.at(i).update_pos(r_tmp);
        particles.at(i).update_vel(v_tmp);
    }
    return;
}

void box::run_sim() {
    double t = 0, t_step = 1e-13;

    while(t < t_max) {
        update_pos(t_step);
        t += t_step;
    }
    ;
}


void box::showParts() {
    double mVel = 0, pVel;
    for(uint i = 0; i < N; i++) {
        mVel = 0;
        for(uint j = 0; j < 3; j++) {
            pVel = particles.at(i).get_acc().at(j);
            cout << setw(10) << setprecision(10) << left << pVel << ' ';
        }
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

vector<double> particle::get_acc() {
    return acc;
}

double particle::get_mass() {
    return mass;
}

void particle::update_pos(vector<double> r) {
    pos = r;
}

void particle::update_vel(vector<double> V) {
    vel = V;
}

void particle::update_acc(vector<double> a) {
    acc = a;
}

void particle::set_mass(double M) {
    mass = M;
}
