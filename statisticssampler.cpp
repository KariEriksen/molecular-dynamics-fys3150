#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include <iostream>

using std::ofstream; using std::cout; using std::endl;

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already
    if(!m_file.is_open()) {
        m_file.open("statistics.txt", ofstream::out);
        // If it's still not open, something bad happened...
        if(!m_file.is_open()) {
            cout << "Error, could not open statistics.txt" << endl;
            exit(1);
        }
    }

    //m_kineticEnergy;
    //for(int i=0; i < system.systemSize(); i++) {
    //m_file <<  m_temperature << " " << m_kineticEnergy << " " << m_potentialEnergy << " " << m_kineticEnergy+m_potentialEnergy << " " << m_diffusion << " ";
    //}
    m_file << system.steps() << " " <<  system.time() << " " << m_temperature << " " << m_kineticEnergy << " " << m_potentialEnergy << " " <<  m_diffusion << " " << m_r_msd << endl;
}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleDiffusion(system);
    saveToFile(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_potentialEnergy = system.potential().potentialEnergy();
}

void StatisticsSampler::sampleTemperature(System &system)
{
    // Hint: reuse the kinetic energy that we already calculated
    m_temperature = (2./3)*m_kineticEnergy/(system.atoms().size());
}

void StatisticsSampler::sampleDensity(System &system)
{
    m_density = system.atoms().size()/(system.systemSize().x()*system.systemSize().y()*system.systemSize().z());
}

void StatisticsSampler::sampleDiffusion(System &system)
{
    vec3 ri_vec;
    double ri;
    for(Atom *atom : system.atoms()) {
        ri_vec = (atom->position - atom->initialPosition);
        ri += ri_vec.lengthSquared();
    }
    m_r_msd = ri/system.atoms().size();
    m_diffusion = m_r_msd/(6*system.time());
}

