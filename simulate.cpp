#include "simulate.h"
#include "const.h"

// def in the main file...i guess it okay
extern double get_rand();

int sign() {
    int n = rand() % 2;
    return (n == 0) ? 1 : -1;
}



// for rectangular box
box::box(double x, double y, double z) { 
    l_x =   x;
    l_y =   y;
    l_z =   z; 
    vol = l_x * l_y * l_z;
    N = 0; 
    t_max = 0;
    recalculate_ring = false; 
    p_offset = 0;
    equillibrium = false;
    init_walls();
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
    p_offset = 0;
    equillibrium = false;
    init_walls();
}

// create all particles in the particles vector
void box::setNumParticles(double n, double mass) {
    N = n;
    particle myPar(mass);
    M = mass;
    vector<particle*> nada;
    for(uint i = 0; i < N; i++) {
        particles.push_back(myPar);
        influence_ring.insert(pair<uint, vector<particle*>>(i, nada));
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
void box::setVelocity() {
    double x,y,z;
    
    for(uint i = 0; i < N; i++) {
        x = get_rand() * 10e10 * sign();
        y = get_rand() * 10e10 * sign();
        z = get_rand() * 10e10 * sign();
        particles.at(i).update_vel({x,y,z});
    }

}

// max time for sim to run from 0 -> t_max
void box::setTime(double tMax, double tstep) {
    t_max = tMax;
    t_step = tstep;
}

void box::set_temp(double T) {
    double t_cur = calc_temp(), t_scale;
    t_scale = sqrt(T/ t_cur);
    for(uint q = 0; q < N; q++) {
        particles.at(q).update_vel( scale(particles.at(q).get_vel(), t_scale) );
    }
}

void box::init_walls() {
    for(uint q = 1; q <= 6; q++) 
        walls.insert( pair<uint,uint>(q,0) );
}

void box::special_setup() {
    particles.at(0).update_pos({l_x/2,0,l_x/2});
    particles.at(1).update_pos({{l_x/2,l_x/2,l_x/2}});
    particles.at(0).update_vel( {0,0.5*10e10,0} );
    particles.at(1).update_vel( {0,-1*0.5*10e10,0} );
}

void box::calc_acceleration() {
    vector<double> total_force, r1, r2;
    vector<particle*> ring_particles;        // particles that influence within rmin < r < rmax
    uint n;

    for(uint i = 0; i < N; i++) 
    {
        total_force = {0,0,0};                                          
        r1 = particles.at(i).get_pos(); 
        ring_particles = influence_ring.at(i);
        n = ring_particles.size();

        for(uint j = 0; j < n; j++) 
        {
            r2 = ring_particles.at(j)->get_pos(); 
            total_force = sum( total_force, calc_force(i,j) );
        }

        if(recalculate_ring) {
            //build_ring(i);
            build_ring_for_all();
            recalculate_ring = false;
        }

        particles.at(i).update_acc(scale(total_force, (1/M)));
    }
}

double box::calc_temp() {
    double T = 0, A = M / (K_B * 3*N);
    vector<double> v;
    for(uint q = 0; q < N; q++) {
        v = particles.at(q).get_vel();
        T += A * dot_product(v,v);
    }
    return T;
}

double box::calc_pressure(double time, double temp) {

    double ideal_part, pair_part, pressure;

    ideal_part = N * K * temp;
    pair_part = p_offset / 3;
    pair_part = pair_part * 10e-20;

    pressure = (ideal_part + pair_part) / (vol*10e-30);
    return pressure * 10e-5;
}

double box::calc_total_momentum() { 

    vector<double> p_tot = {0,0,0};


    for(uint q = 0; q < N; q++) {
        p_tot = sum(p_tot, scale( particles.at(q).get_vel(), M) );
    }
    return 10e10*sqrt(dot_product(p_tot,p_tot));
}

// O(N^2)
void box::build_ring_for_all() {
    double r;
    vector<double> r1, r2;

    for(uint i = 0; i < N; i++) {
        r1 = particles.at(i).get_pos();
        influence_ring.at(i).clear();
        for(uint j = 0; j < N; j++) {

            if(j != i) {
                r2 = particles.at(j).get_pos();
                r2 = sum(r2, scale(r1, -1));            // use r2 again to avoid long function calls or another variable creation
                r = sqrt(dot_product(r2,r2));



                if((r < r_max) && (r > r_min))
                    influence_ring.at(i).push_back(&particles.at(j));
            }
        }
    }
}

// O(N);
void box::build_ring(uint q) {

    vector<double> r2, r1 = particles.at(q).get_pos();
    double r;
    influence_ring.at(q).clear();

    for(uint i = 0; i < N; i++) {

        if( i != q ) {
            r2 = particles.at(i).get_pos();
            r2 = sum(r2, scale(r1, -1));
            r = sqrt(dot_product(r2,r2));
            if((r < r_max) && (r > r_min))
                    influence_ring.at(q).push_back(&particles.at(i));
        }
    }
}

void box::setInfluenceZone(double rmin, double rmax) {
    r_min = rmin;
    r_max = rmax;
    build_ring_for_all();
}


void box::show_ring() {
    uint n;
    for(uint q = 0; q < N; q++) {
        n = influence_ring.at(q).size();
        cout << setw(10) << left << q << setw(10) << left << n << endl;
    }
}

// calculates the force between two particles from the lj potential
vector<double> box::calc_force(uint base, uint neighbor) {
    vector<double> r_hat, r1, r2, v1, v2, q, force;
    double  r, A, B, LJ_F_scale; 

    // retarted!!!!!!
    r1      = particles.at(base).get_pos();
    r2      = influence_ring.at(base).at(neighbor)->get_pos(); 

    v1      = particles.at(base).get_vel();
    v2      = influence_ring.at(base).at(neighbor)->get_vel();                 

    r_hat   = sum(r1, scale(r2, -1));
    r       = sqrt( dot_product(r_hat, r_hat) );
    r_hat   = scale(r_hat, 1/r);
 


    // elastic collistion between particles
    if(r == r_min) {
        q = sum(v1, scale(v2, -1));
        q = scale(r_hat, dot_product(r_hat, q));
        q = scale(q, -2*M);


        influence_ring.at(base).at(neighbor)->update_vel(  sum(v2, scale(q, (-1/M)))   );
        particles.at(base).update_vel(     sum(v1, scale(q, (1/M)))        );
    }

    // if the radius is outside or inside the zone we won't consider
    if(r > r_max) {
        recalculate_ring = true;
    }

    A =     SGMA    / r;
    B =     EPS     / r;  
    LJ_F_scale = 24 * B * ( (2*pow(A, 12)) - pow(A,6) );

    force = scale(r_hat, LJ_F_scale);

    //pressure calculation
    p_offset += dot_product(r_hat, force);


    return scale(r_hat, LJ_F_scale);
}

void box::update_pos(double time_step) {
    vector<double> r,v,a;
    double x, acc;
    p_offset = 0;

    for(uint i = 0; i < N; i++) {
        r   =   particles.at(i).get_pos();
        v   =   particles.at(i).get_vel();
        a   =   particles.at(i).get_acc();
        
        // update postitional vectors before anything else
        for(uint q = 0; q < 3; q++) {
            x = r.at(q);
            acc = a.at(q);
            x = x + (v.at(q)*time_step) + (0.5 * acc * time_step*time_step);
            switch(q) {
                // WALLS PARALLEL TO X-AXIS
                case 0:
                    if(x < 0) {
                        v = reflect_v(v, {1,0,0});
                        walls.at(1) += 1;
                    }
                    else if( x > l_x ) {
                        v = reflect_v(v, {-1,0,0});
                        walls.at(2) += 1;
                    }
                break;

                // WALLS PARALLEL TO Y-AXIS
                case 1:
                    if(x < 0) {
                        v = reflect_v(v, {0,1,0});
                        walls.at(3) += 1;
                    }
                    else if( x > l_y ) {
                        v = reflect_v(v, {0,-1,0});
                        walls.at(4) += 1;
                    }
                break;

                // WALLS PARALLEL TO Z-AXIS
                case 2:
                    if(x < 0) {
                        v = reflect_v(v, {0,0,1});
                        walls.at(5) += 1;
                    }
                    else if( x > l_z) {
                        v = reflect_v(v, {0,0,-1});
                        walls.at(6) += 1;
                    }
                break;
            }
            r.at(q) = r.at(q) + (v.at(q)*time_step) + (0.5 * acc * time_step*time_step);
        }

        calc_acceleration();
        a = particles.at(i).get_acc();

        for(uint q = 0; q < 3; q++)
            v.at(q)  = v.at(q) + (0.5*a.at(q)* time_step);

        particles.at(i).update_pos(r);
        particles.at(i).update_vel(v);
    }

    p_offset_avg = p_offset / N;

}

void box::run_sim() {
    double t = 0, T;
    while(t < t_max) {
        update_pos(t_step);
        p_offset = p_offset_avg;
        T = calc_temp();
        //cout << setw(5) << left << "T = " << setw(15) << setprecision(10) << left << T << setw(5) << left << " P = " << setw(15) << setprecision(10) << left << calc_pressure(t, T) << endl;
        //     << setw(8) << left << " P_tot = " << setw(15) << setprecision(10) << left << calc_total_momentum();
        //cout << setprecision(10) << T << endl;//',' << setprecision(10) << calc_pressure(t, T) << endl;
        cout << setprecision(10) << T << ',' << setprecision(10) << calc_pressure(t,T) << endl;
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