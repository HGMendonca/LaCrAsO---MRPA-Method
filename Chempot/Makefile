FC = gfortran
#FC = g95

#FFLAGS = -O3 -fopenmp
FFLAGS = -mcmodel=medium -O3
#FFLAGS = -O0 -ggdb

#LFLAGS = -O3 -fopenmp
LFLAGS = -mcmodel=medium -O3
#LFLAGS = -O0 -ggdb

#LIBS = -lgomp /usr/lib64/liblapack.a /usr/lib64/libblas.a
#LIBS = /usr/lib64/liblapack.a /usr/lib64/libblas.a
LIBS = /usr/lib/x86_64-linux-gnu/liblapack.a /usr/lib/x86_64-linux-gnu/libblas.a

#DEBUG = -fsanitize=address
DEBUG = 

OBJECTS = interface.o chempot.o tightb.o 

MODULES = interface.mod

DATA = 
#chempot.dat saida.dat

.PHONY: clean

chempot.dat: chempot.exe
	./chempot.exe > saida.dat
#gdb baresus.exe

chempot.exe: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) $(DEBUG) -o chempot.exe $(LIBS)

%.mod : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(MODULES) $(DATA) $(FIGURES) chempot.exe

help:
	@echo "Valid targets:"
	@echo "  chempot.exe"
	@echo "  chempot.o"
	@echo "  clean:  removes .o and .exe files"
