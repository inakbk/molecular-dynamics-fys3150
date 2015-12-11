#include "math/random.h"
#include "potentials/lennardjones.h"
#include "integrators/eulercromer.h"
#include "integrators/velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include "time.h"
#include <string>

using namespace std;

int main(int numberOfArguments, char **argumentList)
{
    clock_t executionTimeStart, executionTimeFinish;
    executionTimeStart = clock();

    int numberOfUnitCells = 5;
    double initialTemperature = UnitConverter::temperatureFromSI(2000.0); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms
    double dt = UnitConverter::timeFromSI(5e-15); // Measured in seconds
    int integratorNumber = 2; //initially set to Velocity Verlet, 1 is Euler-Cromer
    int numberOfTimesteps = 1000;

    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in MD units)
    if(numberOfArguments > 2) initialTemperature = atof(argumentList[2]);
    // If a third argument is provided, it is the lattice constant determining the density (measured in MD units)
    if(numberOfArguments > 3) latticeConstant = atof(argumentList[3]);
    // If a fourth argument is provided, it is the time step length dt (measured in MD units)
    if(numberOfArguments > 4) dt = atof(argumentList[4]);
    // If a fifth argument is provided, it is the integrator number
    if(numberOfArguments > 5) integratorNumber = atoi(argumentList[5]);
    // If a sixth argument is provided, it is the number of timesteps
    if(numberOfArguments > 6) numberOfTimesteps = atoi(argumentList[6]);

    /* // Possible to print the MD units used in the simulation
    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    */

    System system;

    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
    system.setPotential(new LennardJones(3.405, 1.0));

    if(integratorNumber == 1) system.setIntegrator(new EulerCromer());
    else if(integratorNumber == 2) system.setIntegrator(new VelocityVerlet());

    system.removeTotalMomentum();

    StatisticsSampler statisticsSampler;
    IO statisticsFile; // To write statistics to file for plotting
    IO movie; // To write the state to file to look at in Ovito

    if(numberOfArguments > 5){ // making unique filenames for plotting with python
        string filenameStatistics = "run_plot_python_output/statistics_file_NrOfCells" + to_string(numberOfUnitCells) + "_T" + to_string(int(initialTemperature*1000)) + "_b" + to_string(int(latticeConstant*1000)) + "_dt" + to_string(int(dt*10000)) + "_int" + to_string(integratorNumber) + ".txt";
        statisticsFile.open(filenameStatistics.c_str());
        string filenameMovie = "run_plot_python_output/movie_NrOfCells" + to_string(numberOfUnitCells) + "_T" + to_string(int(initialTemperature*1000)) + "_b" + to_string(int(latticeConstant*1000)) + "_dt" + to_string(int(dt*10000)) + "_int" + to_string(integratorNumber) + ".xyz";
        movie.open(filenameMovie.c_str());
    }
    else{
        statisticsFile.open("statistics_file.txt");
        movie.open("movie.xyz");
    }

    cout << "Timestep Time Temperature KineticEnergy PotentialEnergy TotalEnergy DiffusionConstant MeanSquareDisplacement" << endl;
    for(int timestep=0; timestep<numberOfTimesteps; timestep++) {
        //movie.saveState(&system); //including also the starting position in the movie, uncomment to write to movie file

        system.step(dt); //moving the particle one step

        statisticsSampler.sample(system);
        statisticsFile.saveStatistics(&system, &statisticsSampler, 0); //output done here

        if( !(timestep % 100) ) {
            // Print the timestep every 100 timesteps
            cout << system.steps() << "      " << system.time() << "      " << statisticsSampler.temperature() << "      " << statisticsSampler.kineticEnergy() << "      " << statisticsSampler.potentialEnergy() << "      " << statisticsSampler.totalEnergy() << "      " << statisticsSampler.diffusionConstant() << "      " << statisticsSampler.meanSquareDisplacement() << endl;
        }
    }

    movie.close();

    executionTimeFinish = clock();
    double totalExecutionTime = (( executionTimeFinish - executionTimeStart )/double(CLOCKS_PER_SEC ));
    cout << "Execution time: " << totalExecutionTime << endl;

    statisticsFile.saveStatistics(&system, &statisticsSampler, totalExecutionTime);
    statisticsFile.close();

    return 0;
}



