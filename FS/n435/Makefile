FC = gfortran
#FC = g95

#FFLAGS = -O3 -fopenmp
FFLAGS = -O3 -mcmodel=medium
#FFLAGS = -O0 -ggdb

#LFLAGS = -O3 -fopenmp
LFLAGS = -O3 -mcmodel=medium
#LFLAGS = -O0 -ggdb

#LIBS = -lgomp /usr/lib64/liblapack.a /usr/lib64/libblas.a
LIBS = /usr/lib/x86_64-linux-gnu/liblapack.a /usr/lib/x86_64-linux-gnu/libblas.a

DEBUG = -fsanitize=address

OBJECTS = interface.o amu.o band.o tightb.o fermi-n4.25-7p.o

MODULES = interface.mod

FIGURES = 

DATA = pockets-all.dat saida.dat

.PHONY: clean

pockets-all.dat: fermi-n4.25-7p.exe
	./fermi-n4.25-7p.exe > saida.dat
#gdb fermi-n4.25-7p.exe

fermi-n4.25-7p.exe: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) -o fermi-n4.25-7p.exe $(LIBS)

%.mod : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(MODULES) $(DATA) $(FIGURES) fermi-n4.25-7p.exe

help:
	@echo "Valid targets:"
	@echo "  fermi-n4.25-7p.exe"
	@echo "  fermi-n4.25-7p.o"
	@echo "  clean:  removes .o and .exe files"
