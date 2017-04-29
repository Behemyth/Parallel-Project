# Parallel-Project
Parallel Computing Final Project 2017 - RPI

Compile: mpic++ project.cpp -std=c++11 -o project.exe
Run: mpiexec -n (processor count)  ./project.exe (particle count) (ticks) (recursion level) (nearest neighbors)

Program arguments:
Recursion level - determines how smooth the final displayed sphere will be, must be greater than 0
					6 creates a very smooth sphere