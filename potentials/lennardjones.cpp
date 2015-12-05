#include "lennardjones.h"
#include <cmath>

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop

    for(int i=0; i< system->atoms().size(); i++){
        for(int j=i+1; j< system->atoms().size(); j++){
            Atom *atom_i = system->atoms()[i];
            Atom *atom_j = system->atoms()[j];

            double x_ij = atom_i->position[0] - atom_j->position[0];
            double y_ij = atom_i->position[1] - atom_j->position[1];
            double z_ij = atom_i->position[2] - atom_j->position[2];

            double r_ij = std::sqrt(x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);

            double r = sigma/r_ij;
            double r6 = r*r*r*r*r*r;
            double r8 = r6*r*r;
            double some_numbers = ( 24*m_epsilon/(m_sigma*m_sigma) )*r8*( 2*r6 - 1);

            atom_i->force[0] += some_numbers*x_ij;
            atom_i->force[1] += some_numbers*y_ij;
            atom_i->force[2] += some_numbers*z_ij;

            atom_j->force[0] -= some_numbers*x_ij;
            atom_j->force[1] -= some_numbers*y_ij;
            atom_j->force[2] -= some_numbers*z_ij;
        }
    }
}
