#include "velocityverlet.h"
#include "../system.h"
#include "../atom.h"

void VelocityVerlet::integrate(System *system, double dt)
{
    system->calculateForces();
    for(Atom *atom : system->atoms()) {
        vec3 velocity_half = atom->velocity + atom->force*0.5*dt / atom->mass();
        atom->position += velocity_half*dt;
        vec3 force_next =
        atom->velocity = velocity_half + force_next*0.5*dt / atom->mass();

    }

    system->applyPeriodicBoundaryConditions();
}
