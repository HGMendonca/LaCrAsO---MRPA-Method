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

OBJECTS = interface.o path.o tightb.o 

MODULES = interface.mod

DATA = 
#path.dat saida.dat

.PHONY: clean

path.dat: path.exe
	./path.exe > saida.dat
#gdb path.exe

path.exe: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) $(DEBUG) -o path.exe $(LIBS)

%.mod : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(MODULES) $(DATA) $(FIGURES) path.exe

help:
	@echo "Valid targets:"
	@echo "  cpath.exe"
	@echo "  path.o"
	@echo "  clean:  removes .o and .exe files"
