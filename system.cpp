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
    }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
}

void System::resetForcesOnAllAtoms() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).
    for(int i=0; i<100; i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        double x = Random::nextDouble()*10; // random number in the interval [0,10]
        double y = Random::nextDouble()*10;
        double z = Random::nextDouble()*10;
        atom->position.set(x,y,z);
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);
    }
    setSystemSize(vec3(10, 10, 10));
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
