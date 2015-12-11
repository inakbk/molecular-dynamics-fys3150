#include "system.h"
#include "statisticssampler.h"
#include "potentials/potential.h"

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::sample(System &system)
{
    // Measuring different kinds of statistical properties, saving to file is done in IO class.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleDiffusionConstant(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_potentialEnergy = system.potential()->potentialEnergy();
}

void StatisticsSampler::sampleTemperature(System &system)
{
    //instantaneous temperature:
    m_temperature = 2*m_kineticEnergy/(3*system.atoms().size());
}

void StatisticsSampler::sampleDensity(System &system)
{
    //the (mass)density is constant 4m/b^3, this is particle desity:
    m_density = system.atoms().size() / ( system.systemSize().x()*system.systemSize().y()*system.systemSize().z() );
}

void StatisticsSampler::sampleDiffusionConstant(System &system)
{
    double sumSquaredR = 0;
    for(Atom *atom : system.atoms()) {
        double x = atom->position[0]; // x(t)
        double y = atom->position[1];
        double z = atom->position[2];

        double x0 = atom->initialPosition[0]; // x(0)
        double y0 = atom->initialPosition[1];
        double z0 = atom->initialPosition[2];

        double xi2 = (x - x0)*(x - x0);
        double yi2 = (y - y0)*(y - y0);
        double zi2 = (z - z0)*(z - z0);

        sumSquaredR += xi2 + yi2 - zi2;
    }
    m_meanSquareDisplacement = sumSquaredR/system.atoms().size();
    m_diffusionConstant = m_meanSquareDisplacement/( 6*system.time() );
}




