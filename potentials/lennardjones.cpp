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

    for(Atom *atom_i : m_atoms){
        for(Atom *atom_j : m_atoms){
            double x_ij = atom_i->position[0] - atom_j->position[0];
            double y_ij = atom_i->position[0] - atom_j->position[0];
            double z_ij = atom_i->position[0] - atom_j->position[0];

            double r_ij = std::sqrt(x_ij*x_ij + y_ij*y_ij + z_ij*z_ij);

            double r = sigma/r_ij;
            double r6 = r*r*r*r*r*r;
            double r8 = r6*r*r;
            double some_numbers = ( 24*m_epsilon/(m_sigma*m_sigma) )*r8*( 2*r6 - 1);

            m_force_x = some_numbers*x_ij;
            m_force_y = some_numbers*y_ij;
            m_force_y = some_numbers*z_ij;
        }

    }


}
