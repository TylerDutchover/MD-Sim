#include "simulate.h"

double get_rand() {
    return (static_cast<double>(rand())) / RAND_MAX;
}


int main()
{
    //cout << "size = ";
    //cin >> size;

    //cout << "\nNumOfParticles = ";
    //cin >> numofParts;


    double gas_mass = 6.6e-26;
    srand(time(NULL));
    box testBox(200);

    
    
    testBox.setNumParticles(100, gas_mass);        // requires mass
    testBox.setPositions_random();                  // random postitions
    testBox.setVelocity();                          // set velocity randomly
    testBox.set_temp(400);
    testBox.setTime(1e-10, 1e-15);                  // max time to run sim, and time step
    testBox.setInfluenceZone(5, 25);


    cout << "done.\nCalculating forces...\n";
    testBox.calc_acceleration();
    cout << "done. Starting\n\n";

    //testBox.show_ring();
    //getchar();


    //cout << "x,y,z\n";
    testBox.run_sim();

    return 0;
}