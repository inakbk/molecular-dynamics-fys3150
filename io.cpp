#include "io.h"
#include "system.h"
#include "atom.h"
#include "unitconverter.h"
#include "statisticssampler.h"
#include <cstdlib>
using std::endl; using std::cout;

IO::IO()
{

}

IO::~IO() {
    close();
}

void IO::open(const char *filename) {
    if(file.is_open()) {
        std::cout << "<IO.cpp> Error, tried to open file " << filename << ", but some file is already open." << endl;
        exit(1);
    }

    file.open(filename);
}

void IO::close() {
    if(file.is_open()) {
        file.close();
    }
}

// This saves the current state to a file following the xyz-standard (see http://en.wikipedia.org/wiki/XYZ_file_format )
void IO::saveState(System *system)
{
    if(file.is_open()) {
        file << system->atoms().size() << endl;
        file << "The is an optional comment line that can be empty." << endl;
        for(Atom *atom : system->atoms()) {
            file << UnitConverter::lengthToAngstroms(atom->position.x()) << " " << UnitConverter::lengthToAngstroms(atom->position.y()) << " " << UnitConverter::lengthToAngstroms(atom->position.z()) << endl;
        }
    }
}

// This saves the statistics to a .txt file
void IO::saveStatistics(System *system, StatisticsSampler *statisticsSampler, double totalExecutionTime)
{
    if(file.is_open()) {
        if(first_time){
            file << "Timestep Time Temperature KineticEnergy PotentialEnergy TotalEnergy DiffusionConstant MeanSquareDisplacement" << endl;
            first_time = false;
        }
        //should convert units or not? add units in first line?
        file << system->steps() << "      " << system->time() << "      " << statisticsSampler->temperature() << "      " << statisticsSampler->kineticEnergy() << "      " << statisticsSampler->potentialEnergy() << "      " << statisticsSampler->totalEnergy() << "      " << statisticsSampler->diffusionConstant() << "      " << statisticsSampler->meanSquareDisplacement() << endl;
//      file << UnitConverter::lengthToAngstroms(atom->position.x()) << " " << UnitConverter::lengthToAngstroms(atom->position.y()) << " " << UnitConverter::lengthToAngstroms(atom->position.z()) << endl;

        if(totalExecutionTime){
            file << "Execution time: " << totalExecutionTime << endl;
        }
    }
}



