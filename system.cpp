#include "system.h"
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()
{

}

System::~System()
{
    delete m_potential;
    delete m_integrator;
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention

    for(Atom *atom : m_atoms) {
        double x = atom->position[0];
        double y = atom->position[1];
        double z = atom->position[2];
        if(x<0) x += m_systemSize[0];
        if(y<0) y += m_systemSize[1];
        if(z<0) z += m_systemSize[2];
        if(x>m_systemSize[0]) x -= m_systemSize[0];
        if(y>m_systemSize[1]) y -= m_systemSize[1];
        if(z>m_systemSize[2]) z -= m_systemSize[2];
        atom->position[0] = x;
        atom->position[1] = y;
        atom->position[2] = z;

        //r_0 must also be corrected here
    }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.

    vec3 average_velocity = vec3(0,0,0);

    for(Atom *atom : m_atoms) {
        average_velocity[0] += atom->velocity[0];
        average_velocity[1] += atom->velocity[1];
        average_velocity[2] += atom->velocity[2];
    }

    average_velocity = average_velocity/atoms().size();
    //std::cout << "Average velocity before removing total momentum " << average_velocity << std::endl;

    // removing total momentum in the velocities (assuming all masses are equal)
    for(Atom *atom : m_atoms) {
        atom->velocity[0] = atom->velocity[0] - average_velocity[0];
        atom->velocity[1] = atom->velocity[1] - average_velocity[1];
        atom->velocity[2] = atom->velocity[2] - average_velocity[2];
    }

    //testing implementation
//    average_velocity = vec3(0,0,0);
//    for(Atom *atom : m_atoms) {
//        average_velocity[0] += atom->velocity[0];
//        average_velocity[1] += atom->velocity[1];
//        average_velocity[2] += atom->velocity[2];
//    }
//    std::cout << "Average velocity after removing total momentum " << average_velocity/atoms().size() << std::endl;
}

void System::resetForcesOnAllAtoms() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).

    for(int i=0; i<numberOfUnitCellsEachDimension-1; i++) {
        for(int j=0; j<numberOfUnitCellsEachDimension-1; j++){
            for(int k=0; k<numberOfUnitCellsEachDimension-1; k++){
                // one unit cell:
                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                double x = 0 + latticeConstant*i;
                double y = 0 + latticeConstant*j;
                double z = 0 + latticeConstant*k;
                atom1->position.set(x,y,z);
                atom1->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom1);

                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = latticeConstant*0.5 + latticeConstant*i;
                y = latticeConstant*0.5 + latticeConstant*j;
                z = 0 + latticeConstant*k;
                atom2->position.set(x,y,z);
                atom2->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom2);

                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = 0 + latticeConstant*i;
                y = latticeConstant*0.5 + latticeConstant*j;
                z = latticeConstant*0.5 + latticeConstant*k;
                atom3->position.set(x,y,z);
                atom3->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom3);

                Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = latticeConstant*0.5 + latticeConstant*i;
                y = 0 + latticeConstant*j;
                z = latticeConstant*0.5 + latticeConstant*k;
                atom4->position.set(x,y,z);
                atom4->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom4);
            }
        }
    }
    setSystemSize(vec3(numberOfUnitCellsEachDimension, numberOfUnitCellsEachDimension, numberOfUnitCellsEachDimension));
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
}

void System::step(double dt) {
    m_integrator->integrate(this, dt);
    m_steps++;
    m_time += dt;
}
