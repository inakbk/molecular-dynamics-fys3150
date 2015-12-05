#include "lennardjones.h"

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop

    x_ij = ;
    y_ij = ;
    z_ij = ;

    r_ij = ;

    r = sigma/r_ij;
    r6 = r*r*r*r*r*r;
    r8 = r6*r*r;
    some_numbers = ( 24*m_epsilon/(m_sigma*m_sigma) )*r8*( 2*r6 - 1);

    m_force_x = some_numbers*x_ij;
    m_force_y = some_numbers*y_ij;
    m_force_y = some_numbers*z_ij;

}
