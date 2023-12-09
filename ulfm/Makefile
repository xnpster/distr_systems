MPI_DIR=$(HOME)/distr_systems/improve_task/mpi-src/install/bin/
MPIF_NAME=mpif90
MPIRUN_NAME=mpirun

NUM_PROC=4


MPIF=$(MPI_DIR)/$(MPIF_NAME)
MPIRUN=$(MPI_DIR)/$(MPIRUN_NAME)
MPIFFLAGS= -O3 -fopenmp
MPIRUN_FLAGS= -n $(NUM_PROC) --with-ft ulfm --map-by :OVERSUBSCRIBE

all:
	$(MPIF) $(MPIFFLAGS) -c iter.f90 -o iter.o
	$(MPIF) $(MPIFFLAGS) -c print_ranks.f90 -o print_ranks.o
	
	$(MPIF) $(MPIFFLAGS) -c ex.f90 -o ex.o
	$(MPIF) $(MPIFFLAGS) -c slave.f90 -o slave.o
	
	$(MPIF) $(MPIFFLAGS) ex.o print_ranks.o iter.o -o ex
	$(MPIF) $(MPIFFLAGS) slave.o print_ranks.o iter.o -o slave

	$(MPIRUN) $(MPIRUN_FLAGS) ./ex

clean:
	rm -rf ex slave *.o *.mod backups/*