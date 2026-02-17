# Makefile using g++ for pardy.cpp

CXX=mpicxx
CXXFLAGS=-std=c++17 -fopenmp -O3 -Wall -Wfatal-errors -g -flto -march=native -ffast-math -I./rarray -DNDEBUG
LDFLAGS=-g -fopenmp -O3 -flto
LDLIBS=

APPNAME=pardy

OBJFILES=$(APPNAME).o lcg.o lattice.o inifile.o atom.o system.o global.o estimates.o parallelwork.o forces.o cells.o mpicommunication.o

TESTOUTFILES=testpardy*.out

all: $(APPNAME)

.PHONY: clean test

$(APPNAME): $(APPNAME).o lcg.o lattice.o atom.o system.o inifile.o global.o estimates.o parallelwork.o forces.o cells.o mpicommunication.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(APPNAME).o: $(APPNAME).cpp lattice.h atom.h system.h inifile.h global.h estimates.h parallelwork.h forces.h cells.h mpicommunication.h

atom.o: atom.cpp atom.h
cells.o: cells.cpp cells.h atom.h system.h forces.h global.h parallelwork.h debug.h
estimates.o: estimates.cpp global.h estimates.h
forces.o: forces.cpp forces.h global.h atom.h parallelwork.h 
global.o: global.cpp global.h
inifile.o: inifile.cpp inifile.h
lattice.o: lattice.cpp lattice.h
lcg.o: lcg.cpp lcg.h
mpicommunication.o: mpicommunication.cpp mpicommunication.h atom.h system.h parallelwork.h cells.h forces.h global.h debug.h estimates.h
parallelwork.o: parallelwork.cpp parallelwork.h atom.h
pardy.o: pardy.cpp lcg.h lattice.h debug.h system.h atom.h inifile.h global.h estimates.h parallelwork.h forces.h cells.h mpicommunication.h
system.o: system.cpp system.h

clean:
	\rm -f $(OBJFILES) $(TESTOUTFILES)

test: test_pardy_27 test_pardy_503_nt1_np1 test_pardy_503_nt2_np1 test_pardy_503_nt4_np1 test_pardy_503_nt1_np2 test_pardy_503_nt2_np2 test_pardy_503_nt1_np3 test_pardy_503_nt1_np4 test_pardy_503_nt1_np5 test_pardy_503_nt1_np8 test_pardy_512_nt1_np1 test_pardy_512_nt1_np2 test_pardy_512_nt1_np4 test_pardy_512_nt1_np8 test_pardy_32768_nt1_np1 test_pardy_32768_nt1_np2 test_pardy_32768_nt1_np3 test_pardy_32768_nt1_np4 test_pardy_32768_nt1_np5 test_pardy_32768_nt1_np8 test_pardy_262144_nt2_np2

test_pardy_27: $(APPNAME) lj27.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 1 ./$(APPNAME) lj27.ini > testpardy27.out

test_pardy_503_nt1_np1: $(APPNAME) lj503.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 1 ./$(APPNAME) lj503.ini > testpardy503_1_1.out

test_pardy_503_nt2_np1: $(APPNAME) lj503.ini
	OMP_NUM_THREADS=2 time -p mpirun -np 1 ./$(APPNAME) lj503.ini > testpardy503_2_1.out

test_pardy_503_nt2_np2: $(APPNAME) lj503.ini
	OMP_NUM_THREADS=2 time -p mpirun -np 2 ./$(APPNAME) lj503.ini > testpardy503_2_2.out

test_pardy_503_nt4_np1: $(APPNAME) lj503.ini
	OMP_NUM_THREADS=4 time -p mpirun -np 1 ./$(APPNAME) lj503.ini > testpardy503_4_1.out

test_pardy_503_nt1_np2: $(APPNAME) lj503.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 2 ./$(APPNAME) lj503.ini > testpardy503_1_2.out

test_pardy_503_nt1_np3: $(APPNAME) lj503.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 3 ./$(APPNAME) lj503.ini > testpardy503_1_3.out

test_pardy_503_nt1_np4: $(APPNAME) lj503.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 4 ./$(APPNAME) lj503.ini > testpardy503_1_4.out

test_pardy_503_nt1_np5: $(APPNAME) lj503.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 5 --use-hwthread-cpus ./$(APPNAME) lj503.ini > testpardy503_1_5.out

test_pardy_503_nt1_np8: $(APPNAME) lj503.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 8 --use-hwthread-cpus ./$(APPNAME) lj503.ini > testpardy503_1_8.out

test_pardy_512_nt1_np1: $(APPNAME) lj512.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 1 ./$(APPNAME) lj512.ini > testpardy512_1_1.out

test_pardy_512_nt1_np2: $(APPNAME) lj512.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 2 ./$(APPNAME) lj512.ini > testpardy512_1_2.out

test_pardy_512_nt1_np4: $(APPNAME) lj512.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 4 ./$(APPNAME) lj512.ini > testpardy512_1_4.out

test_pardy_512_nt1_np8: $(APPNAME) lj512.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 8 --use-hwthread-cpus ./$(APPNAME) lj512.ini > testpardy512_1_8.out

test_pardy_32768_nt1_np1: $(APPNAME) lj32768.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 1 ./$(APPNAME) lj32768.ini > testpardy32768_1_1.out

test_pardy_32768_nt1_np2: $(APPNAME) lj32768.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 2 ./$(APPNAME) lj32768.ini > testpardy32768_1_2.out

test_pardy_32768_nt1_np3: $(APPNAME) lj32768.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 3 ./$(APPNAME) lj32768.ini > testpardy32768_1_3.out

test_pardy_32768_nt1_np4: $(APPNAME) lj32768.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 4 ./$(APPNAME) lj32768.ini > testpardy32768_1_4.out

test_pardy_32768_nt1_np5: $(APPNAME) lj32768.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 5 --use-hwthread-cpus ./$(APPNAME) lj32768.ini > testpardy32768_1_5.out

test_pardy_32768_nt1_np8: $(APPNAME) lj32768.ini
	OMP_NUM_THREADS=1 time -p mpirun -np 8 --use-hwthread-cpus ./$(APPNAME) lj32768.ini > testpardy32768_1_8.out

test_pardy_262144_nt2_np2: $(APPNAME) lj262144.ini
	OMP_NUM_THREADS=2 time -p mpirun -np 2 ./$(APPNAME) lj262144.ini > testpardy262144_1_2.out
