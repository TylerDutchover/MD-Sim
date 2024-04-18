#include "simulate.h"

double get_rand() {
    return (static_cast<double>(rand())) / RAND_MAX;
}


int main()
{
    double gas_mass = 6.6e-26;
    srand(time(NULL));
    box testBox(1e-8);

    
    testBox.setNumParticles(100, gas_mass);   // requires mass
    testBox.setPositions_random();          // random postitions
    testBox.updateSizeBox(1e-5);
    testBox.setTemperature(300);            // set velocity st |v| = v_rms
    testBox.setTime(1e-9, 1e-15);                 // max time to run sim
    testBox.setAcc();
    cout << "x,y,z\n";
    testBox.run_sim();

    return 0;
}