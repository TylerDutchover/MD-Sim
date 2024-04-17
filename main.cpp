#include "simulate.h"

double get_rand() {
    return (static_cast<double>(rand())) / RAND_MAX;
}


int main()
{
    double gas_mass = 6.6e-26;
    srand(time(NULL));
    box testBox(1e-4);

    cout << "creating particles...";
    testBox.setNumParticles(5, gas_mass);   // requires mass
    cout << "done.\n";

    cout << "setting up positions...";
    testBox.setPositions_random();          // random postitions
    cout << "done.\n";

    cout << "initializing velocities...";
    testBox.setTemperature(300);            // set velocity st |v| = v_rms
    cout << "done.\n";

    testBox.setTime(1e-6);                 // max time to run sim

    cout << "calculating forces...";
    testBox.setAcc();                       // set acc vector according to LJ force
    cout << "done.\n";

    cout << "Starting simulations.\n\n";
    testBox.run_sim();
    //testBox.showParts();

    return 0;
}