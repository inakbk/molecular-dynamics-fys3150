#include "lennardjones.h"
#include <cmath>

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0;
    double sigma6 = pow(m_sigma, 6);
    double rCut = 2.5*m_sigma;
    double energyAtRcut = 4*m_epsilon*sigma6*pow(rCut,-6)*(sigma6*pow(rCut,-6) - 1);

    for(int i=0; i< system->atoms().size(); i++){
        Atom *atom_i = system->atoms()[i];
        for(int j=i+1; j< system->atoms().size(); j++){
            Atom *atom_j = system->atoms()[j];

            double x_ij = atom_i->position[0] - atom_j->position[0];
            double y_ij = atom_i->position[1] - atom_j->position[1];
            double z_ij = atom_i->position[2] - atom_j->position[2];

            // minimum image convention:
            if (x_ij >   system->systemSize()[0] * 0.5) x_ij = x_ij - system->systemSize()[0];
            else if (x_ij <= -system->systemSize()[0] * 0.5) x_ij = x_ij + system->systemSize()[0];
            if (y_ij >   system->systemSize()[1] * 0.5) y_ij = y_ij - system->systemSize()[1];
            else if (y_ij <= -system->systemSize()[1] * 0.5) y_ij = y_ij + system->systemSize()[1];
            if (z_ij >   system->systemSize()[2] * 0.5) z_ij = z_ij - system->systemSize()[2];
            else if (z_ij <= -system->systemSize()[2] * 0.5) z_ij = z_ij + system->systemSize()[2];

            double r2 = x_ij*x_ij + y_ij*y_ij + z_ij*z_ij;
            if(r2 > rCut*rCut) continue;

            double oneOverR2 = 1.0/r2;
            double oneOverR6 = oneOverR2*oneOverR2*oneOverR2;

            double F = 24*m_epsilon*sigma6*oneOverR6*(2.0*sigma6*oneOverR6 - 1.0)*oneOverR2;

            // calculating final force
            atom_i->force[0] += F*x_ij;
            atom_i->force[1] += F*y_ij;
            atom_i->force[2] += F*z_ij;

            atom_j->force[0] -= F*x_ij;
            atom_j->force[1] -= F*y_ij;
            atom_j->force[2] -= F*z_ij;

            // calculating total potential energy in same loop (moved with the energy at rCut)
            m_potentialEnergy += 4*m_epsilon*sigma6*oneOverR6*(sigma6*oneOverR6 - 1) - energyAtRcut;
        }
    }
}
