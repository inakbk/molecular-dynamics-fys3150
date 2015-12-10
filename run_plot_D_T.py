#program to run and compile cpp code and plot

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
#plotting statistics vs time for one dt

#all in MD units!
numberOfUnitCells = 5
latticeConstant = 5.26 # measured in angstroms (MD)

dt = 0.05 # Measured MD units, is about the same as 5e-15 seconds
integratorNumber = 2 # initially set to Velocity Verlet, 1 is Euler-Cromer

numberOfTimesteps = 10000
initialTemperature_list = linspace(0.5,18.5,55) # measured in MD approx 2.5=300 Kelvin, 8.35 = 1000 Kelvin, 16.7 = 2000 Kelvin
instantTemperatureEquilibrium = zeros(len(initialTemperature_list))
diffusionConstantEquilibrium = zeros(len(initialTemperature_list))
i = 0

#compiling once:
#os.system('g++ -o main *.cpp math/*.cpp potentials/*.cpp integrators/*.cpp -I. -O3 -std=c++11')

for initialTemperature in initialTemperature_list:
    #running cpp code
    #os.system('./main %s %s %s %s %s %s' %(numberOfUnitCells, initialTemperature, latticeConstant, dt, integratorNumber, numberOfTimesteps))
    time, instantTemperature, kineticEnergy, potentialEnergy, totalEnergy, diffusionConstant, meanSquareDisplacement = read_file('run_plot_python_output/statistics_file_NrOfCells%s_T%s_b%s_dt%s_int%s.txt' %(numberOfUnitCells, int(initialTemperature*1000), int(latticeConstant*1000), int(round(dt*10000)), integratorNumber))
    instantTemperatureEquilibrium[i] = instantTemperature[-1]
    diffusionConstantEquilibrium[i] = diffusionConstant[-1]
    i += 1

    figure(1)
    plot(time[500:], diffusionConstant[500:])
    xlabel('time [MD units]', fontsize=18)
    ylabel('$D$ [MD units]', fontsize=18)
    #legend(['Velocity Verlet'], fontsize=16, loc='lower left')
    title('numberOfUnitCells= %s, initialTemperature= [%s,%s] \n latticeConstant= %s, numberOfTimesteps= %s,dt= %s' %(numberOfUnitCells, initialTemperature_list[0], initialTemperature_list[-1], latticeConstant, len(time), dt), fontsize=16)

    figure(2)
    plot(time, meanSquareDisplacement)
    xlabel('time [MD units]', fontsize=18)
    ylabel(r'$\langle r^2(t) \rangle$ [MD units]', fontsize=18)
    #legend(['Velocity Verlet'], fontsize=16, loc='lower left')
    title('numberOfUnitCells= %s, initialTemperature= [%s,%s] \n latticeConstant= %s, numberOfTimesteps= %s, dt= %s' %(numberOfUnitCells, initialTemperature_list[0], initialTemperature_list[-1], latticeConstant, len(time), dt), fontsize=16)

    figure(3)
    plot(time[:500], instantTemperature[:500]/initialTemperature)
    xlabel('Time', fontsize=18)
    ylabel(r'$T/T_i$ $[199.735 K]$', fontsize=18)
    #legend(['Velocity Verlet'], fontsize=16, loc='lower left')
    title('numberOfUnitCells= %s, initialTemperature= [%s,%s] \n latticeConstant= %s, numberOfTimesteps= %s, dt= %s' %(numberOfUnitCells, initialTemperature_list[0], initialTemperature_list[-1], latticeConstant, len(time), dt), fontsize=16)

figure(4)
plot(instantTemperatureEquilibrium, instantTemperatureEquilibrium/initialTemperature_list, '-o')
xlabel('$T$ $[199.735 K]$', fontsize=18)
ylabel(r'$T/T_i$ $[199.735 K]$', fontsize=18)
#legend(['Velocity Verlet'], fontsize=16, loc='lower left')
title('numberOfUnitCells= %s, initialTemperature= [%s,%s] \n latticeConstant= %s, numberOfTimesteps= %s, dt= %s' %(numberOfUnitCells, initialTemperature_list[0], initialTemperature_list[-1], latticeConstant, len(time), dt), fontsize=16)

figure(5)
plot(instantTemperatureEquilibrium, diffusionConstantEquilibrium, '-o')
xlabel('$T$ $[199.735 K]$', fontsize=18)
ylabel(r'$D$ [MD units]', fontsize=18)
#legend(['Velocity Verlet'], fontsize=16, loc='lower left')
title('numberOfUnitCells= %s, initialTemperature= [%s,%s] \n latticeConstant= %s, numberOfTimesteps= %s, dt= %s' %(numberOfUnitCells, initialTemperature_list[0], initialTemperature_list[-1], latticeConstant, len(time), dt), fontsize=16)

show()
