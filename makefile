all: 
	g++ main.cpp OdeSolver.cpp -I/usr/include/python3.8 -lpython3.8 -fopenmp -g -Wall -O2 -o main 
	./main 
	python3 plot_odes.py