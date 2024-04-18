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
void box::setTime(double tMax, double tstep) {
    t_max = tMax;
    t_step = tstep;
}

void box::setAcc() {

    vector<double> total_force, r1, r2, accl;
    double M = particles.at(0).get_mass();

    for(uint i = 0; i < N; i++) 
    {
        total_force = {0,0,0};                                          // total force on a particle due to the other N-1 particles
        r1 = particles.at(i).get_pos();                                 // grab positional vector of particle

        for(uint j = 0; j < N; j++) 
        {
            if(j != i) 
            {
                r2 = particles.at(j).get_pos();                         // grab second postitional vector
                total_force = sum(total_force, calc_force(r1,r2));
            }
        }

        accl = scale(total_force, (1/M));
        particles.at(i).update_acc(accl);
    }
}

// calculates the force between two particles from the lj potential
vector<double> box::calc_force(vector<double> r1, vector<double> r2) {
    vector<double> r_hat;
    double  r,                              // magnitude of seperation vector        
            A,                              // sigma    / r
            B,                              // epsilon  / r
            LJ_F_scale;                     // prefactor for LJ force
    

    r_hat = sum(r2, scale(r1, -1));
    r = sqrt( dot_product(r_hat, r_hat) );
    r_hat = scale(r_hat, 1/r);

    A =     SGMA    / r;
    B =     EPS     / r;  
    LJ_F_scale = 24 * B * ( (2*pow(A, 12)) - pow(A,6) );

    return scale(r_hat, LJ_F_scale);
}

void box::update_pos(double time_step) {
    vector<double> r_tmp, v_tmp, a1, a0;
    double x;

    for(uint i = 0; i < N; i++) {
        r_tmp   = particles.at(i).get_pos();
        v_tmp   = particles.at(i).get_vel();
        a1      = particles.at(i).get_acc();
        a0      = a1;
        
        if(i == 0) {
            //for(uint q = 0; q < 3; q++) 
            for(uint q = 0; q < 3; q++)
            {
                if(q != 2)
                    cout << r_tmp.at(q) << ',';
                else
                    cout << r_tmp.at(q) << '\n'; 
            }
            //for(auto a : v_tmp)
            //    cout << setw(15) << setprecision(10) <<  right << a << setw(2) << ' ';
            cout << endl;
        }
    


        // update postitional vectors before anything else
        for(uint q = 0; q < 3; q++) {
            x = r_tmp.at(q);
            x = x + (v_tmp.at(q)*time_step) + (0.5*a1.at(q) * pow(time_step,2));
            switch(q) {
                case 0:
                    if(x < 0) 
                        v_tmp = reflect_v(v_tmp, {1,0,0});
                    else if( x > l_x)
                        v_tmp = reflect_v(v_tmp, {-1,0,0});
                break;
                case 1:
                    if(x < 0) 
                        v_tmp = reflect_v(v_tmp, {0,1,0});
                    else if( x > l_y)
                        v_tmp = reflect_v(v_tmp, {0,-1,0});
                break;
                case 2:
                    if(x < 0) 
                        v_tmp = reflect_v(v_tmp, {0,0,1});
                    else if( x > l_x)
                        v_tmp = reflect_v(v_tmp, {0,0,-1});
                break;
            }
            r_tmp.at(q) = r_tmp.at(q) + (v_tmp.at(q)*time_step) + (0.5*a1.at(q) * pow(time_step,2));
        }
        setAcc();                                   // recalculate the forces


        a1      =   particles.at(i).get_acc();      // grab the updated acceleration

        for(uint q = 0; q < 3; q++)
            v_tmp.at(q)  = v_tmp.at(q) + (0.5* (a1.at(q) + a0.at(q)) * time_step);

        // finally update particle position and velocity
        particles.at(i).update_pos(r_tmp);
        particles.at(i).update_vel(v_tmp);
    }
}

void box::run_sim() {
    double t = 0;
    while(t < t_max) {
        update_pos(t_step);
        t += t_step;
    }
    ;
}

vector<double> box::sum(vector<double> v1, vector<double> v2) {
    uint n = v1.size();
    for(uint q = 0; q < n; q++)
        v1.at(q) = v1.at(q) + v2.at(q);
    return v1;
}

vector<double> box::scale(vector<double> v, double a) {
    uint n = v.size();
    for(uint q = 0; q < n; q++) 
        v.at(q) = a * v.at(q);
    return v;
}

double box::dot_product(vector<double>v1, vector<double> v2) {
    uint n = v1.size();
    double product = 0;
    for(uint q = 0; q < n; q++)
        product += v1.at(q)*v2.at(q);
    return product;
}

vector<double> box::reflect_v(vector<double> v, vector<double> n) {
    vector<double> m = scale(n, -2 * dot_product(v,n));
    v = sum(v, m);
    return v;
}

void box::showParts() {
    return;
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