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

OBJECTS = interface.o rpasus.o sus.o amu.o tightb.o

MODULES = interface.mod

DATA = rpasus.dat saida.dat

.PHONY: clean

rpasus.dat: rpasus.exe
	./rpasus.exe > saida.dat
#gdb rpasus.exe

rpasus.exe: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) $(DEBUG) -o rpasus.exe $(LIBS)

%.mod : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(MODULES) $(DATA) rpasus.exe

help:
	@echo "Valid targets:"
	@echo "  rpasus.exe"
	@echo "  rpasus.o"
	@echo "  clean:  removes .o and .exe files"
