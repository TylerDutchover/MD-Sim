#include "simulate.h"

double get_rand() {
    return (static_cast<double>(rand())) / RAND_MAX;
}


int main()
{
    double gas_mass = 6.6e-26;
    srand(time(NULL));
    box testBox(10000); // 50 angstroms

    
    testBox.setNumParticles(8000, gas_mass);     // requires mass
    testBox.setPositions_random();              // random postitions
    testBox.setVelocity(300);                   // set velocity st |v| <= v_rms
    testBox.setTime(1e-10, 1e-15);               // max time to run sim, and time step
    testBox.setInfluenceZone(5, 10);
    testBox.setAcc();
    cout << "x,y,z\n";
    testBox.run_sim();

    return 0;
}