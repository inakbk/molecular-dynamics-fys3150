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

    totalEnergy = zeros(len(all_lines))

    i = 1
    while i < len(all_lines) - 1:
    	totalEnergy[i] = float(all_lines[i].split()[5])

    	i += 1

    return totalEnergy

"""
------------------------------------------------------------------------------------------
"""

numberOfUnitCells = 5
initialTemperature = 300.0 # measured in Kelvin
latticeConstant = 5.26 # measured in angstroms
#dt = 0.001 # Measured MD units, is about the same as 1e-15 seconds
integratorNumber = 2 # initially set to Velocity Verlet, 1 is Euler-Cromer

#compiling once:
#os.system('g++ -o main *.cpp math/*.cpp potentials/*.cpp integrators/*.cpp -I. -O3 -std=c++11')

N = 10
dt_list = linspace(0.0001, 0.1, N)
sigmaEnergyEuler = zeros(N)
sigmaEnergyVV = zeros(N)
i = 0

for dt in dt_list:
	os.system('./main %s %s %s %s %s' %(numberOfUnitCells, initialTemperature, latticeConstant, dt, 1))
	totalEnergy = read_file('run_plot_python_output/statistics_file.txt')
	sigmaEnergyEuler[i] = sum(totalEnergy**2)/len(totalEnergy) - sum(totalEnergy)**2/len(totalEnergy)
	
	os.system('./main %s %s %s %s %s' %(numberOfUnitCells, initialTemperature, latticeConstant, dt, 2))
	totalEnergy = read_file('run_plot_python_output/statistics_file.txt')
	sigmaEnergyVV[i] = sum(totalEnergy**2)/len(totalEnergy) - sum(totalEnergy)**2/len(totalEnergy)
	
	i += 1

#print len(totalEnergy)
#print sigmaEnergyEuler #/((1.651*10**(-21))**2) what units?

plot(dt_list, sigmaEnergyEuler, 'k-o')
hold('on')
plot(dt_list, sigmaEnergyVV, 'r-o')
show()




