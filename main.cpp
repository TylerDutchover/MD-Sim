#include "simulate.h"

double get_rand() {
    return (static_cast<double>(rand())) / RAND_MAX;
}


int main()
{
    srand(time(NULL));
    box testBox(50);

    testBox.setNumParticles(5);
    testBox.setPositions_random();

    testBox.showParts();

    return 0;
}