#include "velocityverlet.h"
#include "../system.h"
#include "../atom.h"

void VelocityVerlet::integrate(System *system, double dt)
{
    if(first_time) {
        system->calculateForces();
        first_time = false;
    }

    for(Atom *atom : system->atoms()) {
        atom->velocity = atom->velocity + atom->force*0.5*dt / atom->mass();
        atom->position += atom->velocity*dt;
    }
    system->applyPeriodicBoundaryConditions();
    // calculating forces again now that the position is updated:
    system->calculateForces();
    for(Atom *atom : system->atoms()){
        // calculating the new velocity with the new force:
        atom->velocity = atom->velocity + atom->force*0.5*dt / atom->mass();
    }
}
