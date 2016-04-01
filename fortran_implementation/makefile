#makefile for qo code, using mac Accelerate framework
RM = rm -f

OBJS = t2g_hamiltonian_class.o linalg.o get_command_line_input.o solver_class.o main_recursive.o 

COMP = ifort

#LDFLAGS=-framework Accelerate
LDFLAGS = -L/opt/intel/compilers_and_libraries_2016.2.146/mac/mkl/lib/intel64/ -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -qopenmp -lmkl_core -lpthread

qo: $(OBJS)
	$(COMP) $(FLAGS) $(LDFLAGS) -o qo  $(OBJS) 

main_recursive.o: main_recursive.f90 t2g_hamiltonian_class.mod solver_class.mod get_command_line_input.mod
	$(COMP) $(FLAGS) $(LDFLAGS) -c $<


t2g_hamiltonian_class.mod t2g_hamiltonian_class.o: t2g_hamiltonian_class.f90
	$(COMP) $(FLAGS) $(LDFLAGS) -c $<

linalg.mod linalg.o: linalg.f90
	$(COMP) $(FLAGS) $(LDFLAGS) -c $<

get_command_line_input.mod get_command_line_input.o: get_command_line_input.f90
	$(COMP) $(FLAGS) $(LDFLAGS) -c $<

solver_class.mod solver_class.o: solver_class.f90 t2g_hamiltonian_class.mod linalg.mod
	$(COMP) $(FLAGS) $(LDFLAGS) -c $<

clean:
	$(RM) *.o *.mod Thubb2

debug: $(OBJS)
	$(COMP) $(OBJS) $(LDFLAGS) -g  -o qo_debug