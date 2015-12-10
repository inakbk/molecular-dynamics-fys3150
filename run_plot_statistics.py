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
numberOfUnitCells = [5,10]
initialTemperature = 8.35 # measured in MD approx 2.5=300 Kelvin, 8.35 = 1000 Kelvin, 16.7 = 2000 Kelvin
latticeConstant = 5.26 # measured in angstroms (MD)

dt = 0.05 # Measured MD units, is about the same as 5e-15 seconds
integratorNumber = 2 # initially set to Velocity Verlet, 1 is Euler-Cromer
numberOfTimesteps = 10000

#compiling once:
#os.system('g++ -o main *.cpp math/*.cpp potentials/*.cpp integrators/*.cpp -I. -O3 -std=c++11')

#running cpp code
#os.system('./main %s %s %s %s %s %s' %(numberOfUnitCells, initialTemperature, latticeConstant, dt, integratorNumber, numberOfTimesteps))

time, instantTemperature, kineticEnergy, potentialEnergy, totalEnergy, diffusionConstant, meanSquareDisplacement = read_file('run_plot_python_output/timesteps10000/statistics_file_NrOfCells%s_T%s_b%s_dt%s_int%s.txt' %(numberOfUnitCells[0], int(initialTemperature*1000), int(latticeConstant*1000), int(dt*10000), integratorNumber))
time2, instantTemperature2, kineticEnergy2, potentialEnergy2, totalEnergy2, diffusionConstant2, meanSquareDisplacement2 = read_file('run_plot_python_output/timesteps10000/statistics_file_NrOfCells%s_T%s_b%s_dt%s_int%s.txt' %(numberOfUnitCells[1], int(initialTemperature*1000), int(latticeConstant*1000), int(dt*10000), integratorNumber))

time = time[0:1000]
kineticEnergy = kineticEnergy[0:1000]
potentialEnergy = potentialEnergy[0:1000]
totalEnergy = totalEnergy[0:1000]
instantTemperature = instantTemperature[0:1000]

time2 = time2[0:1000]
kineticEnergy2 = kineticEnergy2[0:1000]
potentialEnergy2 = potentialEnergy2[0:1000]
totalEnergy2 = totalEnergy2[0:1000]
instantTemperature2 = instantTemperature2[0:1000]

figure(1)
subplot(2,1,1)
plot(time, instantTemperature, 'r')
ylabel('T $[199.735 K]$', fontsize=18)
legend(['numberOfUnitCells= %s' %numberOfUnitCells[0]], fontsize=16, loc='upper right')
title('numberOfUnitCells= %s, initialTemperature= %s \n latticeConstant= %s, numberOfTimesteps= %s, dt= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(time), dt), fontsize=16)

subplot(2,1,2)
plot(time2, instantTemperature2, 'b')
xlabel(r'Time $[1.00224\cdot 10^{-13} \mathrm{s}]$', fontsize=18)
ylabel(r'T $[199.735 K]$', fontsize=18)
legend(['numberOfUnitCells= %s' %numberOfUnitCells[1]], fontsize=16, loc='upper right')

figure(2)
subplot(2,1,1)
plot(time, kineticEnergy, 'b')
plot(time, potentialEnergy, 'g')
plot(time, totalEnergy, 'r')
ylabel(r'Energy $[1.651\cdot 10^{-21} \mathrm{J}]$', fontsize=18)
legend([r'$\langle E_k \rangle$', '$V$', '$E$'], fontsize=16, loc='center right')
title('numberOfUnitCells= %s, initialTemperature= %s \n latticeConstant= %s, numberOfTimesteps= %s,dt= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(time), dt), fontsize=16)

subplot(2,1,2)
plot(time2, kineticEnergy2, 'b')
plot(time2, potentialEnergy2, 'g')
plot(time2, totalEnergy2, 'r')

xlabel(r'Time $[1.00224\cdot 10^{-13} \mathrm{s}]$', fontsize=18)
ylabel(r'Energy $[1.651\cdot 10^{-21} \mathrm{J}]$', fontsize=18)
legend([r'$\langle E_k \rangle$', '$V$', '$E$'], fontsize=16, loc='center right')


"""
figure(3)
plot(time, potentialEnergy)
xlabel('time [MD units]', fontsize=18)
ylabel('$V$ [MD units]', fontsize=18)
legend(['Velocity Verlet'], fontsize=16, loc='lower left')
title('numberOfUnitCells= %s, initialTemperature= %s \n latticeConstant= %s, numberOfTimesteps= %s,dt= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(time), dt), fontsize=16)

figure(4)
plot(time, totalEnergy)
xlabel('time [MD units]', fontsize=18)
ylabel('$E$ [MD units]', fontsize=18)
legend(['Velocity Verlet'], fontsize=16, loc='lower left')
title('numberOfUnitCells= %s, initialTemperature= %s \n latticeConstant= %s, numberOfTimesteps= %s,dt= %s' %(numberOfUnitCells, initialTemperature, latticeConstant, len(time), dt), fontsize=16)

figure(5)
plot(time[500:], diffusionConstant[500:])
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
"""
show()
