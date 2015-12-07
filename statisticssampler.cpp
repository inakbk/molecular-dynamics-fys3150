#include "system.h"
#include "statisticssampler.h"
#include "potentials/potential.h"

StatisticsSampler::StatisticsSampler()
{

}

//void StatisticsSampler::saveToFile(System &system)
//{
//    // Save the statistical properties for each timestep for plotting etc.
//}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    //saveToFile(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }

//    for(int i=0; i<system.atoms().size(); i++) {
//        Atom *atom = system.atoms()[i];

//    }
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
    //the (mass)density is constant 4m/b^3, but this is particle desity:
    m_density = system.atoms().size() / ( system.systemSize().x()*system.systemSize().y()*system.systemSize().z() );
}

void StatisticsSampler::sampleDiffusionConstant(System &system)
{
    for(Atom *atom : system.atoms()) {
        //sumSquaredR +=
    }

    //m_diffusionConstant = sumSquaredR/( 6*system.atoms().size()*system.time() ); //add here
}




