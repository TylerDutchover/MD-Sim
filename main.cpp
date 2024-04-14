#include "simulate.h"

double get_rand() {
    return (static_cast<double>(rand())) / RAND_MAX;
}


int main()
{
    double gas_mass = 6.6e-26;
    srand(time(NULL));
    box testBox(50);

    testBox.setNumParticles(5, gas_mass);
    testBox.setPositions_random();
    testBox.setTemperature(300);

    testBox.showParts();

    return 0;
}