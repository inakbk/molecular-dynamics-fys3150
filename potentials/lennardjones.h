#ifndef LENNARDJONES_H
#define LENNARDJONES_H
#include "potential.h"

class LennardJones : public Potential
{
private:
    double m_sigma = 1.0;
    double m_epsilon = 1.0;

    double m_force_x = 0;
    double m_force_y = 0;
    double m_force_z = 0;

public:
    LennardJones(double sigma, double epsilon);
    virtual void calculateForces(System *system);
};
#endif
