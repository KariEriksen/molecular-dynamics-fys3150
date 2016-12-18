#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"
#include <iostream>
using namespace std;

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();

}


void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention


    for (int i=0; i < m_atoms.size(); i++){

        Atom* atom = m_atoms.at(i);
        vec3 position = atom->position;
        vec3 initPosition = atom->initialPosition;
        //vec3 actPosition = atom->position;
        //vec3 counter;
        double x = position.x();
        double y = position.y();
        double z = position.z();

        double xi = initPosition.x();
        double yi = initPosition.y();
        double zi = initPosition.z();

        double xSize = m_systemSize.x();
        double ySize = m_systemSize.y();
        double zSize = m_systemSize.z();

        if (x < 0) { x = x + xSize, xi += xSize; }
        if (y < 0) { y = y + ySize, yi += ySize; }
        if (z < 0) { z = z + zSize, zi += zSize; }

        if (x >= xSize) { x = x - xSize, xi -= xSize; }
        if (y >= ySize) { y = y - ySize, yi -= ySize; }
        if (z >= zSize) { z = z - zSize, zi -= zSize; }

        atom->position.set(x,y,z);
        atom->initialPosition.set(xi, yi, zi);
    }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.

    vec3 momentum(0,0,0);

    for (int i=0; i < m_atoms.size(); i++){
        Atom* atom = m_atoms.at(i);
        vec3 velocity = atom->velocity;
        momentum += atom->mass()*velocity;
    }

    momentum /= m_atoms.size();
    for (int i=0; i < m_atoms.size(); i++){
        m_atoms.at(i)->velocity -= momentum / m_atoms.at(i)->mass();
    }

    vec3 momentumTotal(0,0,0);

    for (int i=0; i < m_atoms.size(); i++){
        Atom* atom = m_atoms.at(i);
        vec3 velocity = atom->velocity;
        momentumTotal += atom->mass()*velocity;
    }
    //std::cout << momentumTotal << std::endl;
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).

    /*for(int i=0; i<100; i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        double x = Random::nextDouble(0, 10); // random number in the interval [0,10]
        double y = Random::nextDouble(0, 10);
        double z = Random::nextDouble(0, 10);
        atom->position.set(x,y,z);
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);
    }
    */
    int Nx, Ny, Nz;
    Nx = numberOfUnitCellsEachDimension;
    Ny = numberOfUnitCellsEachDimension;
    Nz = numberOfUnitCellsEachDimension;
    double b = latticeConstant;
    double b2 = b/2.;
    Random::randomSeed();

    for (int i=0; i < Nx; i++){
        for (int j=0; j < Ny; j++){
            for(int k=0; k < Nz; k++){

                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));

                atom1->position.set(i*b,      j*b,      k*b);
                atom2->position.set(i*b + b2, j*b + b2, k*b);
                atom3->position.set(i*b,      j*b + b2, k*b + b2);
                atom4->position.set(i*b + b2, j*b,      k*b + b2);

                atom1->initialPosition = atom1->position;
                atom2->initialPosition = atom2->position;
                atom3->initialPosition = atom3->position;
                atom4->initialPosition = atom4->position;

                double boltzmannConstant = 1.0; // In these units, the boltzmann constant equals 1
                double standardDeviation = sqrt(boltzmannConstant*temperature/atom1->mass());

                atom1->velocity.set(Random::nextGaussian(0, standardDeviation), Random::nextGaussian(0, standardDeviation), Random::nextGaussian(0, standardDeviation));
                atom2->velocity.set(Random::nextGaussian(0, standardDeviation), Random::nextGaussian(0, standardDeviation), Random::nextGaussian(0, standardDeviation));
                atom3->velocity.set(Random::nextGaussian(0, standardDeviation), Random::nextGaussian(0, standardDeviation), Random::nextGaussian(0, standardDeviation));
                atom4->velocity.set(Random::nextGaussian(0, standardDeviation), Random::nextGaussian(0, standardDeviation), Random::nextGaussian(0, standardDeviation));

                m_atoms.push_back(atom1);
                m_atoms.push_back(atom2);
                m_atoms.push_back(atom3);
                m_atoms.push_back(atom4);

            }
        }
    }
    setSystemSize(vec3(Nx*b, Ny*b, Nz*b)); // Remember to set the correct system size!
}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
