#program to run and compile cpp code

from pylab import *
import os as os

"""
------------------------------------------------------------------------------------------
"""
def read_file(filename):
    infile = open(filename, "r")
    all_lines = infile.readlines()
    infile.close()

    time = zeros(len(all_lines)-2)
    instantTemperature = zeros(len(all_lines)-2)
    kineticEnergy = zeros(len(all_lines)-2)
    potentialEnergy = zeros(len(all_lines)-2)
    totalEnergy = zeros(len(all_lines)-2)
    diffusionConstant = zeros(len(all_lines)-2)
    meanSquareDisplacement = zeros(len(all_lines)-2)

    i = 0
    while i < len(all_lines) - 2:
        time[i] = float(all_lines[i+1].split()[1])
        instantTemperature[i] = float(all_lines[i+1].split()[2])
        kineticEnergy[i] = float(all_lines[i+1].split()[3])
        potentialEnergy[i] = float(all_lines[i+1].split()[4])
    	totalEnergy[i] = float(all_lines[i+1].split()[5])
        diffusionConstant[i] = float(all_lines[i+1].split()[6])
        meanSquareDisplacement[i] = float(all_lines[i+1].split()[7])

    	i += 1
    return time, instantTemperature, kineticEnergy, potentialEnergy, totalEnergy, diffusionConstant, meanSquareDisplacement


def read_file_E(filename):
    infile = open(filename, "r")
    all_lines = infile.readlines()
    infile.close()

    totalEnergy = zeros(len(all_lines)-2)
    executionTime = float(all_lines[-1].split()[4])

    i = 0
    while i < len(all_lines) - 2:
        totalEnergy[i] = float(all_lines[i+1].split()[5])

        i += 1

    return totalEnergy, executionTime

"""
------------------------------------------------------------------------------------------
"""
#plotting sigma_E vs dt:
"""
#all in MD units!
numberOfUnitCells = 10
initialTemperature = 2.5 # measured in Kelvin
latticeConstant = 5.26 # measured in angstroms

length_of_list = 20
dt_list = linspace(0.005, 0.1, length_of_list)

#compiling once:
os.system('g++ -o main *.cpp math/*.cpp potentials/*.cpp integrators/*.cpp -I. -O3 -std=c++11')

#running cpp code
for dt in dt_list:
	os.system('./main %s %s %s %s %s' %(numberOfUnitCells, initialTemperature, latticeConstant, dt, 1))
	os.system('./main %s %s %s %s %s' %(numberOfUnitCells, initialTemperature, latticeConstant, dt, 2))

sigmaEnergyEuler = zeros(length_of_list-1)
sigmaEnergyVV = zeros(length_of_list-1)
executionTimeEuler = zeros(length_of_list-1)
executionTimeVV = zeros(length_of_list-1)

for i in range(length_of_list-1):
    totalEnergy, executionTimeEuler[i] = read_file_E('run_plot_python_output/statistics_file_NrOfCells%s_T%s_b%s_dt%s_int%s.txt' %(numberOfUnitCells, int(initialTemperature*1000), int(latticeConstant*1000), int(round(dt_list[i]*10000)), 1))
    sigmaEnergyEuler[i] = sqrt(sum(totalEnergy**2)/len(totalEnergy) - (sum(totalEnergy)/len(totalEnergy))**2)
    totalEnergy, executionTimeVV[i] = read_file_E('run_plot_python_output/statistics_file_NrOfCells%s_T%s_b%s_dt%s_int%s.txt' %(numberOfUnitCells, int(initialTemperature*1000), int(latticeConstant*1000), int(round(dt_list[i]*10000)), 2))
    sigmaEnergyVV[i] = sqrt(sum(totalEnergy**2)/len(totalEnergy) - (sum(totalEnergy)/len(totalEnergy))**2)

figure(1)
plot(dt_list[:-1], sigmaEnergyEuler, 'k-o')
hold('on')
plot(dt_list[:-1], sigmaEnergyVV, 'r-o')
xlabel('dt [MD units]', fontsize=18)
ylabel('$\sigma_E$ [MD units]', fontsize=18)
legend(['Euler', 'Velocity Verlet'], fontsize=16, loc='upper left')
title('numberOfUnitCells= %s, initialTemperature= %s\n latticeConstant= %s, numberOfTimesteps= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(totalEnergy)), fontsize=16)

#plotting execution time vs dt:
figure(3)
plot(dt_list[:-1], executionTimeEuler, 'k-o')
hold('on')
plot(dt_list[:-1], executionTimeVV, 'r-o')
xlabel('dt [MD units]', fontsize=18)
ylabel('Execution Time [seconds]', fontsize=18)
legend(['Euler', 'Velocity Verlet'], fontsize=16, loc='upper left')
title('numberOfUnitCells= %s, initialTemperature= %s\n latticeConstant= %s, numberOfTimesteps= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(totalEnergy)), fontsize=16)
"""

"""
------------------------------------------------------------------------------------------
"""
#plotting execution time vs numberOfParticles:

#all in MD units!
initialTemperature = 2.5 # measured in Kelvin
latticeConstant = 5.26 # measured in angstroms
dt = 0.05

numberOfUnitCells_list = [5, 10, 15, 20]

executionTimeEuler = zeros(len(numberOfUnitCells_list))
executionTimeVV = zeros(len(numberOfUnitCells_list))

#compiling once:
#os.system('g++ -o main *.cpp math/*.cpp potentials/*.cpp integrators/*.cpp -I. -O3 -std=c++11')

i = 0
#running cpp code
for numberOfUnitCells in numberOfUnitCells_list:
    os.system('./main %s %s %s %s %s' %(numberOfUnitCells, initialTemperature, latticeConstant, dt, 1))
    os.system('./main %s %s %s %s %s' %(numberOfUnitCells, initialTemperature, latticeConstant, dt, 2))
    totalEnergy, executionTimeEuler[i] = read_file_E('run_plot_python_output/statistics_file_NrOfCells%s_T%s_b%s_dt%s_int%s.txt' %(numberOfUnitCells, int(initialTemperature*1000), int(latticeConstant*1000), int(round(dt*10000)), 1))
    totalEnergy, executionTimeVV[i] = read_file_E('run_plot_python_output/statistics_file_NrOfCells%s_T%s_b%s_dt%s_int%s.txt' %(numberOfUnitCells, int(initialTemperature*1000), int(latticeConstant*1000), int(round(dt*10000)), 2))
    i += 1

figure(4)
plot(numberOfUnitCells_list, executionTimeEuler, 'k-o')
hold('on')
plot(numberOfUnitCells_list, executionTimeVV, 'r-o')
xlabel('Number Of Unit Cells', fontsize=18)
ylabel('Execution Time [seconds]', fontsize=18)
legend(['Euler', 'Velocity Verlet'], fontsize=16, loc='upper left')
title('dt= %s, initialTemperature= %s\n latticeConstant= %s, numberOfTimesteps= %s' %(dt, initialTemperature, latticeConstant, len(totalEnergy)), fontsize=16)


"""
------------------------------------------------------------------------------------------
""" 
#plotting statistics for one dt
"""
#all in MD units!
numberOfUnitCells = 15
initialTemperature = 2.5 # measured in MD approx 2.5=300 Kelvin
latticeConstant = 5.26 # measured in angstroms (MD)

dt = 0.0526 # Measured MD units, is about the same as 5e-15 seconds
integratorNumber = 2 # initially set to Velocity Verlet, 1 is Euler-Cromer
numberOfTimesteps = 1000

#compiling once:
#os.system('g++ -o main *.cpp math/*.cpp potentials/*.cpp integrators/*.cpp -I. -O3 -std=c++11')

#running cpp code
#os.system('./main %s %s %s %s %s %s' %(numberOfUnitCells, initialTemperature, latticeConstant, dt, integratorNumber, numberOfTimesteps))

time, instantTemperature, kineticEnergy, potentialEnergy, totalEnergy, diffusionConstant, meanSquareDisplacement = read_file('run_plot_python_output/statistics_file_NrOfCells%s_T%s_b%s_dt%s_int%s.txt' %(numberOfUnitCells, int(initialTemperature*1000), int(latticeConstant*1000), int(dt*10000), integratorNumber))

figure(2)
plot(time, meanSquareDisplacement)
xlabel('time [MD units]', fontsize=18)
ylabel(r'$\langle r^2(t) \rangle$ [MD units]', fontsize=18)
legend(['Velocity Verlet'], fontsize=16, loc='lower left')
title('numberOfUnitCells= %s, initialTemperature= %s \n latticeConstant= %s, numberOfTimesteps= %s, time= %s, dt= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(time), time[-1]), dt, fontsize=16)
"""
show()

