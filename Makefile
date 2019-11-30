SHELL = /bin/sh
CC=g++
CFLAGS=-std=c++11 -c -Wall
MAIN=main.cpp
SOLVER=solver.cpp 
MATRIX=matrix.cpp
MAIN=main.cpp
OBJECTS=main.o solver.o matrix.o
EXECUTABLE=Solver

.PHONY: all clean


all: $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXECUTABLE)

solver.o: $(SOLVER)
	$(CC) $(CFLAGS) $(SOLVER) -o $@

matrix.o: $(MATRIX)
	$(CC) $(CFLAGS) $(MATRIX) -o $@

main.o: $(MAIN)
	$(CC) $(CFLAGS) $(MAIN) -o $@

clean:
	rm ./*.o $(EXECUTABLE)



