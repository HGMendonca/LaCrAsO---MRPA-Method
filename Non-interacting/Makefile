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

OBJECTS = interface.o amu.o baresus.o tightb.o 

MODULES = interface.mod

DATA = baresubacr.dat saida.dat

.PHONY: clean

baresubacr.dat: baresus.exe
	./baresus.exe > saida.dat
#gdb baresus.exe

baresus.exe: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) $(DEBUG) -o baresus.exe $(LIBS)

%.mod : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(MODULES) $(DATA) baresus.exe

help:
	@echo "Valid targets:"
	@echo "  baresus.exe"
	@echo "  baresus.o"
	@echo "  clean:  removes .o and .exe files"
