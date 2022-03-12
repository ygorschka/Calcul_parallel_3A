MPI=mpif90
FLAGS=-fbounds-check
EXE=exe
all:$(EXE) clean

$(EXE): param.o communications.o fonctions.o charge.o matrice.o gradient.o main.o
	$(MPI) -o $@ $^
%.o: %.f90
	$(MPI) $(FLAGS) -c $<
clean:
	rm -f *.o *.mod
