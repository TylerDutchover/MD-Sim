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
    recalculate_ring = false; 
}

// for a square box
box::box(double x) {
    l_x =   x;
    l_y =   x;
    l_z =   x;
    vol = l_x * l_y * l_z;
    N = 0;
    t_max = 0;
    recalculate_ring = false;
}

// create all particles in the particles vector
void box::setNumParticles(double n, double mass) {
    N = n;
    particle myPar(mass);
    M = mass;
    vector<particle> nada;
    for(uint i = 0; i < N; i++) {
        particles.push_back(myPar);
        influence_ring.insert(pair<uint, vector<particle>>(i, nada));
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
void box::setVelocity(double T) {
    double v_rms,x,y,z;
    v_rms = (3*K_B*T) / particles.at(0).get_mass();
    v_rms = sqrt(v_rms);                             // in angstrom / sec

    
    for(uint i = 0; i < N; i++) {
        x = get_rand() * v_rms;
        y = get_rand() * v_rms;
        z = get_rand() * v_rms;
        particles.at(i).update_vel({x,y,z});
    }
}

// max time for sim to run from 0 -> t_max
void box::setTime(double tMax, double tstep) {
    t_max = tMax;
    t_step = tstep;
}

void box::setAcc() {
    vector<double> total_force, r1, r2;
    vector<particle> ring_particles;        // particles that influence within rmin < r < rmax
    double n;

    for(uint i = 0; i < N; i++) 
    {
        total_force = {0,0,0};                                          
        r1 = particles.at(i).get_pos(); 
        ring_particles = influence_ring.at(i);
        n = ring_particles.size();

        for(uint j = 0; j < n; j++) 
        {
            r2 = ring_particles.at(j).get_pos();             
            total_force = sum(total_force, calc_force(r1,r2));
        }

        if(recalculate_ring) {
            build_ring(i);
            recalculate_ring = false;
        }

        particles.at(i).update_acc(scale(total_force, (1/M)));
    }
}

// O(N^2)
void box::build_ring_for_all() {
    double r;
    vector<double> r1, r2;

    for(uint i = 0; i < N; i++) {
        r1 = particles.at(i).get_pos();
        influence_ring.at(i).clear();
        for(uint j = 0; j < N; j++) {

            if(i ^ j){
                r2 = particles.at(j).get_pos();
                r2 = sum(r2, scale(r1, -1));            // use r2 again to avoid long function calls or another variable creation
                r = sqrt(dot_product(r2,r2));

                if((r < r_max) && (r > r_min))
                    influence_ring.at(i).push_back(particles.at(j));
            }
        }
    }
}

// O(N);
void box::build_ring(uint q) {

    vector<double> r2, r1 = particles.at(q).get_pos();
    particle p;
    double r;
    influence_ring.at(q).clear();

    for(uint i = 0; i < N; i++) {

        if( i ^ q ) {
            p = particles.at(i);
            r2 = p.get_pos();
            r2 = sum(r2, scale(r1, -1));
            r = sqrt(dot_product(r2,r2));
            if((r < r_max) && (r > r_min))
                    influence_ring.at(q).push_back(p);
        }
    }
}

void box::setInfluenceZone(double rmin, double rmax) {
    r_min = rmin;
    r_max = rmax;
    build_ring_for_all();
}

// calculates the force between two particles from the lj potential
vector<double> box::calc_force(vector<double> r1, vector<double> r2) {
    
    vector<double> r_hat;
    double  r, A, B, LJ_F_scale;                     

    r_hat = sum(r2, scale(r1, -1));
    r = sqrt( dot_product(r_hat, r_hat) );


    // if the radius is outside or inside the zone we won't consider
    if((r > r_max) || (r < r_max))
        recalculate_ring = true;


    r_hat = scale(r_hat, 1/r);

    A =     SGMA    / r;
    B =     EPS     / r;  
    LJ_F_scale = 24 * B * ( (2*pow(A, 12)) - pow(A,6) );

    return scale(r_hat, LJ_F_scale);
}

void box::update_pos(double time_step) {
    vector<double> r,v,a;
    double x, del_x;

    for(uint i = 0; i < N; i++) {
        r   =   particles.at(i).get_pos();
        v   =   particles.at(i).get_vel();
        a   =   particles.at(i).get_acc();
        
        if(i == 0) {
            //for(uint q = 0; q < 3; q++) 
            for(uint q = 0; q < 3; q++)
            {
                if(q != 2)
                    cout << r.at(q) << ',';
                else
                    cout << r.at(q) << '\n'; 
            }
        }

        // update postitional vectors before anything else
        for(uint q = 0; q < 3; q++) {
            x = r.at(q);
            del_x = (v.at(q)*time_step) + (0.5 * a.at(q) * time_step*time_step);
            x = x + del_x;
            switch(q) {
                case 0:
                    if(x < 0) 
                        v = reflect_v(v, {1,0,0});
                    else if( x > l_x)
                        v = reflect_v(v, {-1,0,0});
                break;
                case 1:
                    if(x < 0) 
                        v = reflect_v(v, {0,1,0});
                    else if( x > l_y)
                        v = reflect_v(v, {0,-1,0});
                break;
                case 2:
                    if(x < 0) 
                        v = reflect_v(v, {0,0,1});
                    else if( x > l_x)
                        v = reflect_v(v, {0,0,-1});
                break;
            }
            r.at(q) = r.at(q) + del_x;
        }

        setAcc();                                   // recalculate the forces
        a = particles.at(i).get_acc();

        for(uint q = 0; q < 3; q++)
            v.at(q)  = v.at(q) + (0.5*a.at(q)* time_step);

        particles.at(i).update_pos(r);
        particles.at(i).update_vel(v);
    }
}

void box::run_sim() {
    double t = 0;
    while(t < t_max) {
        update_pos(t_step);
        t += t_step;
    }
}

vector<double> box::sum(vector<double> v1, vector<double> v2) {
    for(uint q = 0; q < 3; q++)
        v1.at(q) = v1.at(q) + v2.at(q);
    return v1;
}

vector<double> box::scale(vector<double> v, double a) {
    for(uint q = 0; q < 3; q++) 
        v.at(q) = a * v.at(q);
    return v;
}

double box::dot_product(vector<double>v1, vector<double> v2) {
    double product = 0;
    for(uint q = 0; q < 3; q++)
        product += v1.at(q)*v2.at(q);
    return product;
}

vector<double> box::reflect_v(vector<double> v, vector<double> n) {
    vector<double> m = scale(n, -2 * dot_product(v,n));
    return sum(v, m);
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
    acc     = {0,0,0};
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