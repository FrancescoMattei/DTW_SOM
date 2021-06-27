# DTW_SOM
Dynamic Time Warping Self-Organizing Map algorithm

The software implements a Self-Organizing Map (SOM) algorithm that uses a Dynamic Time Warping distance measure (DTW). The program is written in Fortran 90 and has been used to cluster time series of the same length. The aim of the clustering was to reduce the complexity of a satellite chlorophyll a time series dataset which comprised over 37000 series. The dataset complexity reduction allowed the spatio-temporal analysis of the chlorophyll a dynamics in the Mediterranean Sea.
The program allows the customization of several parameters through a configuration file (.cfg). The configuration file let the user define:
1)	Input file name (a .txt comma separated file in which the data are stored)
2)	Best Matching Unit (BMU) output file name (Contains the BMU for each data object)
3)	Final codebooks file name (Contains the final codebooks)
4)	Training info file name (Summary of the training configuration)
5)	Number of objects in the data file
6)	Length of the objects (Variables number)
7)	Number of SOM rows
8)	Number of SOM columns
9)	Neighborhood type: gaussian (0) or bubble (1)
10)	Neighborhood shrinking: exponential (0) or linear (1)
11)	Learning rate decay: exponential (0) or linear (1)
12)	Training epochs (In which radius and learning rate decay from user defined values)
13)	Fine tuning epochs (In which radius and learning rate are fixed to a user defined value and 0.01 respectively)
14)	Radius during fine tuning epochs
15)	Radius minimum value
16)	Radius maximum value
17)	Learning rate minimum value
18)	Learning rate maximum value
19)	Window parameter (Define the maximum allowed distance from the diagonal for the DTW computation Sakoe-Chiba (1978)).

An example of configuration file can be found in the repository.
The Fortran 90 file (DTW_SOM .f90 file) can be compiled with Gfortran compiler and run from the Command Prompt tying the application name followed by the configuration file name (e.g. DTW_SOM configuration_file_example.cfg).
The configuration file can be useful to set up several batch procedure to perform a large number of tests.
