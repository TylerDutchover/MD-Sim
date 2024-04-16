#include "simulate.h"

double get_rand() {
    return (static_cast<double>(rand())) / RAND_MAX;
}


int main()
{
    double gas_mass = 6.6e-26;
    srand(time(NULL));
    box testBox(50);

    testBox.setNumParticles(5, gas_mass);   // requires mass
    testBox.setPositions_random();          // random postitions
    testBox.setTemperature(300);            // set velocity st |v| = v_rms
    testBox.setTime(1e-5);                 // max time to run sim
    testBox.setAcc();                       // set acc vector according to LJ force

    testBox.run_sim();
    //testBox.showParts();

    return 0;
}