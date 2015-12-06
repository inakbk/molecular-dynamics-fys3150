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


    return noe

"""
------------------------------------------------------------------------------------------
"""



#compiling once:
os.system('g++ -o main ../*.cpp -larmadillo -llapack -lblas -L/usr/local/lib -I/usr/local/include -O3 -std=c++11')
os.system('./main')


