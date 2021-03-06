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

        double x0 = atom->initialPosition[0];
        double y0 = atom->initialPosition[1];
        double z0 = atom->initialPosition[2];

        // checking if particle has left the system, if so, placing on oposite side
        // of box and moving initial position for the diffusion constant
        if(x<0){
            x += m_systemSize[0];
            x0 += m_systemSize[0];
        }
        if(y<0){
            y += m_systemSize[1];
            y0 += m_systemSize[1];
        }
        if(z<0){
            z += m_systemSize[2];
            z0 += m_systemSize[2];
        }
        if(x>=m_systemSize[0]){
            x -= m_systemSize[0];
            x0 -= m_systemSize[0];
        }
        if(y>=m_systemSize[1]){
            y -= m_systemSize[1];
            y0 -= m_systemSize[1];
        }
        if(z>=m_systemSize[2]){
            z -= m_systemSize[2];
            z0 -= m_systemSize[2];
        }

        atom->position[0] = x;
        atom->position[1] = y;
        atom->position[2] = z;

        atom->initialPosition[0] = x0;
        atom->initialPosition[1] = y0;
        atom->initialPosition[2] = z0;
    }
}

void System::removeTotalMomentum() {
    // Finds the total momentum and removes momentum equally on each atom so the total momentum becomes zero.

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

    // to test implementation:
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

    for(int i=0; i<numberOfUnitCellsEachDimension; i++) {
        for(int j=0; j<numberOfUnitCellsEachDimension; j++){
            for(int k=0; k<numberOfUnitCellsEachDimension; k++){
                // one unit cell:
                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                double x = latticeConstant*i;
                double y = latticeConstant*j;
                double z = latticeConstant*k;
                atom1->position.set(x,y,z);
                atom1->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom1);

                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = latticeConstant*(0.5 + i);
                y = latticeConstant*(0.5 + j);
                z = latticeConstant*k;
                atom2->position.set(x,y,z);
                atom2->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom2);

                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = latticeConstant*i;
                y = latticeConstant*(0.5 + j);
                z = latticeConstant*(0.5 + k);
                atom3->position.set(x,y,z);
                atom3->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom3);

                Atom *atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = latticeConstant*(0.5 + i);
                y = latticeConstant*j;
                z = latticeConstant*(0.5 + k);
                atom4->position.set(x,y,z);
                atom4->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom4);

                if(firstTimeCreateFFCLattice){
                    atom1->initialPosition.set(x,y,z);
                    atom2->initialPosition.set(x,y,z);
                    atom3->initialPosition.set(x,y,z);
                    atom4->initialPosition.set(x,y,z);

                    firstTimeCreateFFCLattice = false;
                }
            }
        }
    }
    setSystemSize(latticeConstant*numberOfUnitCellsEachDimension*vec3(1, 1, 1));

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
