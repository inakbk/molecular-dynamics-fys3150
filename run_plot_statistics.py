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
initialTemperature = 2.5 # measured in MD approx 2.5=300 Kelvin, 8.35 = 1000 Kelvin, 16.7 = 2000 Kelvin
latticeConstant = 5.26 # measured in angstroms (MD)

dt = 0.05 # Measured MD units, is about the same as 5e-15 seconds
integratorNumber = 2 # initially set to Velocity Verlet, 1 is Euler-Cromer
numberOfTimesteps = 10000

#compiling once:
#os.system('g++ -o main *.cpp math/*.cpp potentials/*.cpp integrators/*.cpp -I. -O3 -std=c++11')

#running cpp code
#os.system('./main %s %s %s %s %s %s' %(numberOfUnitCells, initialTemperature, latticeConstant, dt, integratorNumber, numberOfTimesteps))

time, instantTemperature, kineticEnergy, potentialEnergy, totalEnergy, diffusionConstant, meanSquareDisplacement = read_file('run_plot_python_output/statistics_file_NrOfCells%s_T%s_b%s_dt%s_int%s.txt' %(numberOfUnitCells, int(initialTemperature*1000), int(latticeConstant*1000), int(dt*10000), integratorNumber))
#time, instantTemperature, kineticEnergy, potentialEnergy, totalEnergy, diffusionConstant, meanSquareDisplacement = read_file('run_plot_python_output/statistics_file_NrOfCells%s_T%s_b%s_dt%s_int%s.txt' %(numberOfUnitCells, int(initialTemperature*1000), int(latticeConstant*1000), int(dt*10000), integratorNumber))

print len(totalEnergy)

figure(1)
plot(time, instantTemperature)
xlabel('time [MD units]', fontsize=18)
ylabel('$T_i$ [MD units]', fontsize=18)
legend(['Velocity Verlet'], fontsize=16, loc='lower left')
title('numberOfUnitCells= %s, initialTemperature= %s \n latticeConstant= %s, numberOfTimesteps= %s, dt= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(time), dt), fontsize=16)

figure(2)
plot(time, kineticEnergy)
hold('on')
#xlabel('time [MD units]', fontsize=18)
#ylabel(r'$\langle E_k \rangle$ [MD units]', fontsize=18)
#legend(['Velocity Verlet'], fontsize=16, loc='lower left')
#title('numberOfUnitCells= %s, initialTemperature= %s \n latticeConstant= %s, numberOfTimesteps= %s,dt= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(time), dt), fontsize=16)

#figure(3)
plot(time, potentialEnergy)
#xlabel('time [MD units]', fontsize=18)
#ylabel('$V$ [MD units]', fontsize=18)
#legend(['Velocity Verlet'], fontsize=16, loc='lower left')
#title('numberOfUnitCells= %s, initialTemperature= %s \n latticeConstant= %s, numberOfTimesteps= %s,dt= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(time), dt), fontsize=16)

#figure(4)
plot(time, totalEnergy)
xlabel('time [MD units]', fontsize=18)
ylabel('$E$ [MD units]', fontsize=18)
legend(['Velocity Verlet'], fontsize=16, loc='lower left')
title('numberOfUnitCells= %s, initialTemperature= %s \n latticeConstant= %s, numberOfTimesteps= %s,dt= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(time), dt), fontsize=16)

figure(5)
plot(time, diffusionConstant)
xlabel('time [MD units]', fontsize=18)
ylabel('$D$ [MD units]', fontsize=18)
legend(['Velocity Verlet'], fontsize=16, loc='lower left')
title('numberOfUnitCells= %s, initialTemperature= %s \n latticeConstant= %s, numberOfTimesteps= %s,dt= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(time), dt), fontsize=16)

figure(6)
plot(time, meanSquareDisplacement)
xlabel('time [MD units]', fontsize=18)
ylabel(r'$\langle r^2(t) \rangle$ [MD units]', fontsize=18)
legend(['Velocity Verlet'], fontsize=16, loc='lower left')
title('numberOfUnitCells= %s, initialTemperature= %s \n latticeConstant= %s, numberOfTimesteps= %s, dt= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(time), dt), fontsize=16)

show()
