


all: part1_serial part1_omp part1_mpi part2_serial part2_omp part2_mpi

part1_serial: hwk2part1_serial.cpp
	g++ hwk2part1_serial.cpp -o p1s

part1_omp: hwk2part1_omp.cpp
	g++ -fopenmp hwk2part1_omp.cpp -o p1o

part1_mpi: hwk2part1_mpi.cpp
	mpic++ hwk2part1_mpi.cpp -o p1mpi

part2_serial: hwk2part2_serial.cpp
	g++ hwk2part2_serial.cpp -o p2s


part2_omp: hwk2part2_omp.cpp
	g++ -fopenmp hwk2part2_omp.cpp -o p2o

part2_mpi: hwk2part2_omp.cpp
	 mpic++ hwk2part2_mpi.cpp -o p2mpi
