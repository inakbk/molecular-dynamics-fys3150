#include "velocityverlet.h"
#include "../system.h"
#include "../atom.h"

void VelocityVerlet::integrate(System *system, double dt)
{
    system->calculateForces();
    for(Atom *atom : system->atoms()) {
        vec3 velocity_half = atom->velocity + atom->force*0.5*dt / atom->mass();
        atom->position += velocity_half*dt;

    }
    system->calculateForces();
    for(Atom *atom : system->atoms()){
        atom->velocity = velocity_half + atom->force*0.5*dt / atom->mass();

    }

    system->applyPeriodicBoundaryConditions();
}
